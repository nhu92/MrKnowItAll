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

把完整的 871 物种、353 gene 校对 alignment 放到 `data/backbone/`。公开 SPrOUT 仓库目前只带
50 genes 的演示面板，并不等于完整 871 骨架。

```bash
sprout-ref inspect-backbone data/backbone
sprout-ref sync
sprout-ref search "Abatia rugosa"
sprout-ref build \
  --taxon "Abatia rugosa" \
  --backbone data/backbone \
  --output runs/abatia-rugosa
```

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
export SPROUT_REF_BACKBONE=/absolute/path/to/871-alignments
export SPROUT_REF_WORK_DIR=/absolute/path/to/runs
sprout-ref web --host 127.0.0.1 --port 8000
```

打开 <http://127.0.0.1:8000>。Web API 包括 taxon 搜索、后台构建、状态轮询和 ZIP 下载。

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

## 数据与引用

软件代码为 MIT。Kew、SPrOUT、Angiosperms353 target 和你的 871 骨架各自保留原始许可与引用
要求；构建结果的 `report.json` 保存 Kew release、accession 和处理方法。请引用：

- [Kew Tree of Life Explorer public releases](https://sftp.kew.org/pub/treeoflife/)
- [Baker et al. 2022, PAFTOL platform](https://doi.org/10.1093/sysbio/syab035)
- [Johnson et al. 2019, Angiosperms353](https://doi.org/10.1093/sysbio/syy086)
- [Hu et al. 2026, SPrOUT](https://doi.org/10.64898/2026.02.20.707031)
