from fastapi.testclient import TestClient

from sprout_refbuilder.web import app


def test_home_and_unconfigured_build(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("SPROUT_REF_BACKBONE", str(tmp_path / "missing"))
    client = TestClient(app)
    assert client.get("/").status_code == 200
    health = client.get("/api/health").json()
    assert health["backbone_ready"] is False
    response = client.post("/api/build", json={"taxon": "Abatia rugosa"})
    assert response.status_code == 503
    panel = client.post("/api/panel", json={"level": "family", "parent_taxa": ["Rosales"]})
    assert panel.status_code == 503
