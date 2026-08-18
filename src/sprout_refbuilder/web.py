from __future__ import annotations

import os
import threading
import uuid
from pathlib import Path

from fastapi import BackgroundTasks, FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel, Field

from .builder import ReferenceBuilder
from .kew import KewClient
from .panel import HierarchicalPanelBuilder

app = FastAPI(title="SPrOUT Reference Builder", version="0.1.0")
_jobs: dict[str, dict] = {}
_lock = threading.Lock()


class BuildRequest(BaseModel):
    taxon: str = Field(min_length=3, max_length=120, examples=["Abatia rugosa"])
    minimum_score: float = Field(default=0.50, ge=0.0, le=1.0)
    ml_min_score: float = Field(default=-0.05, ge=-1.0, le=1.0)


class PanelRequest(BaseModel):
    level: str = Field(pattern="^(order|family|genus)$")
    parent_taxa: list[str] = Field(default_factory=list)
    include_species: list[str] = Field(default_factory=list)
    representatives: int = Field(default=2, ge=1, le=20)


def settings() -> tuple[Path, Path, Path, Path | None]:
    cache = Path(os.environ.get("SPROUT_REF_CACHE", Path.home() / ".cache/sprout-ref"))
    backbone = Path(os.environ.get("SPROUT_REF_BACKBONE", "data/backbone"))
    work = Path(os.environ.get("SPROUT_REF_WORK_DIR", "runs"))
    model_value = os.environ.get("SPROUT_REF_ML_MODEL")
    return cache, backbone, work, Path(model_value) if model_value else None


@app.get("/", response_class=HTMLResponse)
def home() -> str:
    return _INDEX


@app.get("/api/health")
def health() -> dict:
    cache, backbone, work, model = settings()
    return {
        "status": "ok",
        "cache": str(cache.resolve()),
        "backbone": str(backbone.resolve()),
        "backbone_ready": backbone.is_dir() and any(backbone.glob("*.fasta")),
        "work_dir": str(work.resolve()),
        "ml_model": str(model.resolve()) if model else None,
    }


@app.get("/api/taxa")
def search_taxa(q: str, limit: int = 12) -> list[dict]:
    cache, _, _, _ = settings()
    with KewClient(cache) as client:
        return [record.to_dict() for record in client.search(q, min(max(limit, 1), 50))]


@app.post("/api/build", status_code=202)
def start_build(request: BuildRequest, tasks: BackgroundTasks) -> dict:
    _, backbone, _, _ = settings()
    if not backbone.is_dir() or not any(backbone.glob("*.fasta")):
        raise HTTPException(
            503,
            "Backbone not configured. Set SPROUT_REF_BACKBONE to the alignment directory.",
        )
    job_id = uuid.uuid4().hex
    with _lock:
        _jobs[job_id] = {"id": job_id, "status": "queued", "taxon": request.taxon}
    tasks.add_task(_run_build, job_id, request)
    return _jobs[job_id]


@app.post("/api/panel", status_code=202)
def start_panel(request: PanelRequest, tasks: BackgroundTasks) -> dict:
    _, backbone, _, _ = settings()
    summary = Path(
        os.environ.get("SPROUT_REF_SPECIES_SUMMARY", backbone / "species_summary.csv")
    )
    if not backbone.is_dir() or not any(backbone.glob("*.fasta")):
        raise HTTPException(503, "Backbone is not configured")
    if not summary.is_file():
        raise HTTPException(503, "Species summary is not configured")
    job_id = uuid.uuid4().hex
    with _lock:
        _jobs[job_id] = {
            "id": job_id,
            "status": "queued",
            "type": "panel",
            "level": request.level,
        }
    tasks.add_task(_run_panel, job_id, request, summary)
    return _jobs[job_id]


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    if job_id not in _jobs:
        raise HTTPException(404, "Unknown job")
    return _jobs[job_id]


@app.get("/api/jobs/{job_id}/download")
def download(job_id: str) -> FileResponse:
    job = _jobs.get(job_id)
    if not job or job.get("status") != "complete":
        raise HTTPException(409, "Job is not complete")
    path = Path(job["bundle"])
    return FileResponse(path, filename=path.name, media_type="application/zip")


def _run_build(job_id: str, request: BuildRequest) -> None:
    cache, backbone, work, model = settings()
    output = work / job_id
    with _lock:
        _jobs[job_id]["status"] = "running"
    try:
        with KewClient(cache) as client:
            report = ReferenceBuilder(
                backbone,
                output,
                minimum_score=request.minimum_score,
                ml_model=model,
                ml_min_score=request.ml_min_score,
            ).build_from_kew(request.taxon, client)
        bundle = Path(str(output) + "_bundle.zip")
        with _lock:
            _jobs[job_id].update(status="complete", report=report.to_dict(), bundle=str(bundle))
    except Exception as exc:  # noqa: BLE001 - surfaced to the local user through job state
        with _lock:
            _jobs[job_id].update(status="failed", error=str(exc))


def _run_panel(job_id: str, request: PanelRequest, summary: Path) -> None:
    cache, backbone, work, _ = settings()
    output = work / job_id
    with _lock:
        _jobs[job_id]["status"] = "running"
    try:
        with KewClient(cache) as client:
            report = HierarchicalPanelBuilder(
                backbone,
                summary,
                output,
                representatives_per_taxon=request.representatives,
            ).build(
                request.level,
                client,
                parent_taxa=request.parent_taxa,
                include_species=request.include_species,
            )
        bundle = Path(str(output) + "_bundle.zip")
        with _lock:
            _jobs[job_id].update(status="complete", report=report, bundle=str(bundle))
    except Exception as exc:  # noqa: BLE001 - surfaced to the local user through job state
        with _lock:
            _jobs[job_id].update(status="failed", error=str(exc))


def main() -> None:
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8000)


_INDEX = """<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SPrOUT Reference Builder</title><style>
:root{color-scheme:dark;--bg:#081512;--panel:#10231d;--ink:#e9f5ef;--muted:#9db9ad;--accent:#6ee7a8;--line:#27443a}
*{box-sizing:border-box}body{margin:0;background:radial-gradient(circle at 15% 5%,#174333 0,transparent 35%),var(--bg);font:16px/1.5 system-ui;color:var(--ink)}
main{max-width:860px;margin:8vh auto;padding:24px}.brand{letter-spacing:.16em;text-transform:uppercase;color:var(--accent);font-weight:800}h1{font-size:clamp(2.2rem,7vw,4.6rem);line-height:.95;margin:.5rem 0 1.25rem;max-width:760px}p{color:var(--muted);max-width:680px}
.panel{margin-top:32px;padding:24px;border:1px solid var(--line);border-radius:18px;background:color-mix(in srgb,var(--panel) 90%,transparent);box-shadow:0 24px 80px #0008}
label{display:block;margin:12px 0 8px;color:var(--muted)}.row{display:flex;gap:10px}input,select{flex:1;min-width:0;padding:14px 16px;border:1px solid var(--line);border-radius:10px;background:#091813;color:var(--ink);font-size:1rem}button,a.button{padding:14px 18px;border:0;border-radius:10px;background:var(--accent);color:#062016;font-weight:800;cursor:pointer;text-decoration:none}pre{white-space:pre-wrap;word-break:break-word;color:#cfe7dc;background:#091813;padding:16px;border-radius:10px;min-height:64px}.fine{font-size:.85rem}.grid{display:grid;grid-template-columns:1fr 2fr;gap:10px}
@media(max-width:600px){main{margin:3vh auto;padding:18px}.row{flex-direction:column}}
</style></head><body><main><div class="brand">Angiosperms353 · PAFTOL · SPrOUT</div>
<h1>构建可信 reference</h1><p>从 Kew 精确物种 recovery 或 order → family → genus 分层骨架生成 SPrOUT <code>-r</code> 可用的 alignment bundle；全程保留自动 QC 与 provenance。</p>
<section class="panel"><form id="form"><label for="taxon">目标物种</label><div class="row"><input id="taxon" placeholder="例如 Abatia rugosa" required><button>构建</button></div></form><p class="fine">系统不会用近缘种 consensus 冒充未测序物种。无精确数据时会停止并返回候选名称。</p><pre id="status">等待输入。</pre><a class="button" id="download" hidden>下载 reference bundle</a></section>
<section class="panel"><form id="panel-form"><div class="grid"><div><label for="level">层级</label><select id="level"><option>order</option><option>family</option><option>genus</option></select></div><div><label for="parents">上一级 calls（逗号分隔）</label><input id="parents" placeholder="Asparagales, Brassicales, Rosales"></div></div><label for="species">可选：整合 Kew 精确物种（二名法，逗号分隔）</label><div class="row"><input id="species" placeholder="Allium sativum, Brassica oleracea"><button>构建 panel</button></div></form><p class="fine">order 可留空；family 输入 orders；genus 输入 families。每个 taxon 默认自动选择两个高覆盖且尽量分散的骨架代表。</p></section>
<script>
const form=document.querySelector('#form'), status=document.querySelector('#status'), dl=document.querySelector('#download');
form.addEventListener('submit',async e=>{e.preventDefault();dl.hidden=true;status.textContent='提交任务…';const r=await fetch('/api/build',{method:'POST',headers:{'content-type':'application/json'},body:JSON.stringify({taxon:document.querySelector('#taxon').value})});const j=await r.json();if(!r.ok){status.textContent=j.detail||JSON.stringify(j);return}poll(j.id)});
document.querySelector('#panel-form').addEventListener('submit',async e=>{e.preventDefault();dl.hidden=true;status.textContent='提交 panel 任务…';const split=id=>document.querySelector(id).value.split(',').map(x=>x.trim()).filter(Boolean);const body={level:document.querySelector('#level').value,parent_taxa:split('#parents'),include_species:split('#species')};const r=await fetch('/api/panel',{method:'POST',headers:{'content-type':'application/json'},body:JSON.stringify(body)});const j=await r.json();if(!r.ok){status.textContent=j.detail||JSON.stringify(j);return}poll(j.id)});
async function poll(id){const r=await fetch('/api/jobs/'+id),j=await r.json();status.textContent=JSON.stringify(j,null,2);if(j.status==='queued'||j.status==='running')setTimeout(()=>poll(id),1500);if(j.status==='complete'){dl.href='/api/jobs/'+id+'/download';dl.hidden=false}}
</script></main></body></html>"""
