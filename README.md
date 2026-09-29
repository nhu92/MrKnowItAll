# SPrOUT Reference Builder

自动从 [Kew Tree of Life / PAFTOL](https://treeoflife.kew.org/) 当前公开 release 获取
Angiosperms353 assembled recoveries，用经过校对的多物种 alignment 骨架做自动 QC，并生成
[SPrOUT](https://github.com/nhu92/SPrOUT) `02_exon_trees.py -r` 可直接使用的 reference 目录。

> 当前是一个可运行的 research prototype。它不会为未测序物种“猜”出一条物种级序列；只有目标
> 物种存在公开 recovery，或用户提供该物种自己的 target-sequencing reads/assembly 时，才会给
> 新序列加上该物种名称。

## 核心区别

- `target file`：HybPiper 从 reads 捕获/组装基因时使用（SPrOUT `01_exons_assembly.py -m`）。
- `reference alignments`：每个 gene 一个多物种 aligned FASTA，供 MAFFT 加入 query exon、建树和
  SPrOUT 分类使用（SPrOUT `02_exon_trees.py -r`）。
- 本项目生成第二种；在 FASTQ 模式下只把 target file 临时交给 HybPiper。

## 安装

```bash
conda env create -f environment.yml
conda activate sprout-refbuilder
pip install -e .
```

用原始 ZIP 安装你提供的完整骨架；程序会检查 alignment 等宽和 353-locus 不变量，并生成
SHA-256 清单。实际审计结果是 188,167 条 locus sequences、821 个 specimen/source entries、
809 个唯一 binomials（不是工作描述中的 871 个唯一物种）。

```bash
sprout-ref install-backbone \
  --archive /path/to/angiosperms_353_v2_interim_targetfile_gene_alignments.zip \
  --destination data/backbone
sprout-ref inspect-backbone data/backbone
sprout-ref sync
sprout-ref search "Abatia rugosa"
sprout-ref build \
  --taxon "Abatia rugosa" \
  --backbone data/backbone \
  --output runs/abatia-rugosa
```

安装时默认应用 `corrections.xlsx` 转录出的 48 条 taxonomy/name replacements，并从全部 loci
排除 9 个系统发育位置异常、疑似 mislabelled 的 specimen。异常记录不会按照其错误落点被强行
改名。规则随软件版本保存在 `resources/taxonomy_corrections.tsv`；命中数、未命中规则和规则文件
SHA-256 写入 `backbone_install.json`。可用 `--taxonomy-corrections custom.tsv` 替换规则，或用
`--no-taxonomy-corrections` 明确安装原始 header。

### 分层 reference panel

完整骨架安装后可按上一级 SPrOUT call 构建更窄的 reference：

```bash
sprout-ref build-panel --level order \
  --backbone data/backbone --summary data/backbone/species_summary.csv \
  --output runs/order

sprout-ref build-panel --level family \
  --parent-taxon Asparagales --parent-taxon Brassicales --parent-taxon Rosales \
  --backbone data/backbone --summary data/backbone/species_summary.csv \
  --output runs/family
```

mix7 受控 demo 见 [examples/mix7](examples/mix7/README.md)，或运行：

```bash
sprout-ref demo-mix7 --backbone data/backbone \
  --summary data/backbone/species_summary.csv --output runs/mix7-demo
```

### 用 SPrOUT call 从 KEW 逐级细化

`build-kew-panel` 直接读取 SPrOUT 输出的候选 TXT 或 prediction CSV，从 Kew Release 4.0
下载 assembled Angiosperms353 recoveries，并构造数量平衡的下一层 reference。每个 family/genus
使用相同数目的代表，避免 SPrOUT 当前按 reference 累加分数时偏向 Kew 采样更多的类群。

```bash
# order_candidates.txt 含 Rosales：构建 Rosales 内各 family 的平衡 reference
sprout-ref build-kew-panel --level family \
  --parent-file sprout_test.order_candidates.txt \
  --backbone data/backbone --output runs/sprout_test-family-kew \
  --representatives 4 --threads 16

# family 结果出来后，为候选 families 构建 genus reference
sprout-ref build-kew-panel --level genus \
  --parent-file sprout_test.family_candidates.txt \
  --backbone data/backbone --output runs/sprout_test-genus-kew \
  --representatives 2 --threads 16

# genus 结果出来后，纳入候选 genera 中每个有 recovery 的物种
sprout-ref build-kew-panel --level species \
  --parent-file sprout_test.genus_candidates.txt \
  --backbone data/backbone --output runs/sprout_test-species-kew \
  --representatives 1 --threads 16
```

每一步的 `ref/` 传给 SPrOUT `02_exon_trees.py -r`，`gene.list.txt` 传给 `-g`。程序先以完整
校对骨架执行自动 QC 和 alignment projection，最终 panel 只保留本层 KEW references；缺少精确
recovery 的物种不会被近缘 consensus 冒名代替。`selection.tsv`、`qc.json` 和
`panel_report.json` 保存下载 accession、接受位点数及完整 provenance。
QC 和 alignment 均按 locus 多进程并行；HPC 上将 `--threads` 设为 Slurm 分配的 CPU 数。每个
MAFFT 子进程固定使用 1 核，避免“进程数 × MAFFT 线程数”的嵌套超配。

### 任意混合 reads 一键鉴定到 genus

TTU/Nocona 上可使用 `examples/hpc/ttu_mix_to_genus.sbatch`。它接收任意 paired
FASTQ/FASTQ.GZ，依次运行 assembly、order、Kew family 和 Kew genus，并输出：

```text
FINAL_GENUS_CANDIDATES.txt
FINAL_GENUS_RANKING.csv
FINAL_GENUS_EVIDENCE.csv
```

层级候选不再只由全局 z-score 决定。程序同时读取每个 exon tree 的距离矩阵，对每条 query
exon 投最近 reference 票，并要求至少两个独立 locus 支持；候选是 locus-vote、z-score 和总体
第一名的并集。这能保留混合样品中被强组分压低全局 z-score 的低丰度类群。

UMD Zaratan 上对原始 SPrOUT 50-gene 示例做隔离 smoke test 的 Slurm 脚本见
[examples/hpc](examples/hpc/README.md)。它使用独立输出目录，并限制嵌套 MAFFT 并发以避免超配。

可选的自动异常模型：`sprout-ref train-qc --backbone data/backbone --output data/qc.joblib`，
然后在 `build` 增加 `--ml-model data/qc.joblib`，或为 Web 设置 `SPROUT_REF_ML_MODEL`。默认只将
anomaly score `< -0.05` 的强异常作为额外 rejection，可用 `--ml-min-score` 在独立验证后调整。

输出：

```text
runs/abatia-rugosa/
├── ref/                         # SPrOUT -r 指向这里；每 gene 一个 aligned FASTA
├── gene.list.txt
├── requested_species_sequences.fasta
├── qc.json
└── report.json
runs/abatia-rugosa_bundle.zip
```

在 SPrOUT 中使用：

```bash
python 02_exon_trees.py \
  -e 02_exon_extracted \
  -r runs/abatia-rugosa/ref \
  -g runs/abatia-rugosa/gene.list.txt \
  -p sample_name -t 16
```

## Web app

```bash
export SPROUT_REF_BACKBONE=/absolute/path/to/353-alignments
export SPROUT_REF_SPECIES_SUMMARY=/absolute/path/to/species_summary.csv
export SPROUT_REF_WORK_DIR=/absolute/path/to/runs
sprout-ref web --host 127.0.0.1 --port 8000
```

打开 <http://127.0.0.1:8000>。Web UI/API 支持精确物种构建、order/family/genus panel、Kew
物种整合、后台状态轮询和 ZIP 下载。

## 本地 reads 模式

```bash
sprout-ref build-reads \
  --taxon "Genus species" --order Order --family Family \
  --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz \
  --target-file Angiosperms353_targetSequences.fasta \
  --backbone data/backbone --output runs/genus-species
```

这条路径需要 `hybpiper`、`bwa`；推荐同时安装 `mafft`。若没有 MAFFT，程序使用 Biopython
pairwise projection 保持骨架列宽，适合开发和小规模验证。

详细原理、QC、机器学习、数据边界与生产部署见
[中文技术设计](docs/technical-design.zh-CN.md)。

## 验证

```bash
pytest
ruff check src tests
```

完整骨架上的 mix7 受控 demo 实测生成三个 353-locus 等宽 reference panel。Kew Release 4.0
精确序列自动 QC 的接受数为：Allium sativum 281、Asparagus officinalis 347、Brassica oleracea
350、Artocarpus heterophyllus 349。该 demo 没有 mix7 FASTQ，因此 order/family calls 是按已知正确
上游结果模拟的；它验证的是分层建库与 Kew integration，不应被表述成盲样分类结果。

## 数据与引用

软件代码为 MIT。Kew、SPrOUT、Angiosperms353 target 和你的校对骨架各自保留原始许可与引用
要求；构建结果的 `report.json` 保存 Kew release、accession 和处理方法。请引用：

- [Kew Tree of Life Explorer public releases](https://sftp.kew.org/pub/treeoflife/)
- [Baker et al. 2022, PAFTOL platform](https://doi.org/10.1093/sysbio/syab035)
- [Johnson et al. 2019, Angiosperms353](https://doi.org/10.1093/sysbio/syy086)
- [Hu et al. 2026, SPrOUT](https://doi.org/10.64898/2026.02.20.707031)
