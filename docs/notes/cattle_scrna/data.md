# Cattle scRNA-seq Data Processing

---

# **1. 原始数据获取与预处理（Raw Data Acquisition & Preprocessing）**

## **1.1 批量下载 SRA（支持断点续传、多线程与完整性检查）**

- 使用 `wget -c` 对所有 AWS 链接进行断点续传下载

- 多线程并发下载（常用 48 线程）

- 对空文件、下载失败文件自动重新尝试

- 可启用 `vdb-validate` 进行 SRA 文件一致性校验

示例流程（简化版）：

```bash
# 批量下载 SRR（伪代码）
while read url; do
    fname=$(basename "$url")
    wget -c "$url" -O "data/raw/${fname}"
done < aws_links.txt

# 可选：完整性校验
vdb-validate data/raw/SRRxxxx.sra
```

输出：`data/raw/*.sra`

---

## **1.2 SRA → FASTQ 转换（fasterq-dump + pigz 并行压缩）**

- 使用 fasterq-dump 解压 `.sra`

- 使用多线程 pigz 压缩，显著提升速度

- 每个样本使用独立临时目录，处理更稳健

- 自动跳过已有转换结果

示例流程：

```bash
# 创建临时目录
mkdir -p tmp/${sample}

# 解压为 fastq
fasterq-dump data/raw/${sample}.sra \
    -O tmp/${sample} \
    --split-files \
    --skip-technical \
    -e 12

# 多线程压缩
pigz -p 12 tmp/${sample}/${sample}_1.fastq
pigz -p 12 tmp/${sample}/${sample}_2.fastq

# 移动至输出目录
mv tmp/${sample}/*.fastq.gz data/fastq/
```

输出：  
`data/fastq/SRRxxxx_1.fastq.gz`  
`data/fastq/SRRxxxx_2.fastq.gz`

---

## **1.3 测序质量评估：FastQC + MultiQC**

- 为每个样本执行 FastQC

- 自动识别 CPU 数量并合理分配线程

- 使用 MultiQC 汇总所有样本质量报告

示例：

```bash
# FastQC
fastqc -t 24 data/fastq/*.fastq.gz \
    -o data/qc_fastqc

# 汇总
multiqc data/qc_fastqc \
    -n multiqc_report.html \
    -o data/qc_fastqc
```

输出：

- 每个样本的 `*_fastqc.html`

- 汇总的 `multiqc_report.html`

---

## **1.4 Cell Ranger v9 计数（对每个样本自动调用）**

- 自动扫描 `data/fastq/cellranger_links/`

- 20 CPU + 90GB RAM 分配给单个样本

- 同时运行最多 5 个样本

- 自动跳过已有输出

示例：

```bash
cellranger count \
  --id="SRRxxxx" \
  --transcriptome="refdata-gex-bos_taurus-2020-A" \
  --fastqs="data/fastq/cellranger_links/SRRxxxx" \
  --sample="SRRxxxx" \
  --localcores=20 \
  --localmem=90
```

输出目录：  
`data/cellranger/SRRxxxx/outs/filtered_feature_bc_matrix/`

---

# **2. 表达矩阵合并与质量控制（Merging & Quality Control）**

## **2.1 合并 Cell Ranger 输出为一个 AnnData**

- 自动识别每个样本的 matrix 目录

- 加入 sample_id 与前缀化的 obs_names

- 使用 scanpy 的 `anndata.concat()` 合并

示例：

```python
import scanpy as sc
from pathlib import Path

paths = sorted(Path("data/cellranger").glob("SRR*/outs/filtered_feature_bc_matrix"))

adatas = []
for p in paths:
    sample = p.parts[-3]    # 例如 SRRxxxx
    ad = sc.read_10x_mtx(p, var_names="gene_symbols")
    ad.var_names_make_unique()

    # 添加样本信息
    ad.obs["sample_id"] = sample
    ad.obs_names = [f"{sample}_" + bc for bc in ad.obs_names]

    adatas.append(ad)

adata = adatas[0].concatenate(adatas[1:], join="outer", fill_value=0)
```

输出：`01_raw_merged.h5ad`

---

## **2.2 计算 QC 指标**

- per-cell：基因数、UMI 数

- per-gene：被检测的细胞数

- 计算线粒体比例

示例：

```python
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    inplace=True
)
```

---

## **2.3 过滤低质量细胞**

常用标准：

- `n_genes_by_counts >= 200`

- `n_genes_by_counts <= 6000`

- `total_counts >= 1000`

- `pct_counts_mt < 10`

- 基因需在 ≥3 个细胞表达

示例：

```python
adata = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["n_genes_by_counts"] <= 6000) &
    (adata.obs["pct_counts_mt"] < 10) &
    (adata.obs["total_counts"] >= 1000),
    :
].copy()

sc.pp.filter_genes(adata, min_cells=3)
```

---

## **2.4 Scrublet 双细胞检测（按 sample 分组）**

- 为每个 sample 分别运行 Scrublet

- 添加 doublet_score / predicted_doublet

- 最终过滤掉预测双细胞

示例：

```python
import scrublet as scr
import numpy as np

adata.obs["predicted_doublet"] = False

for sid in adata.obs["sample_id"].unique():
    sub = adata[adata.obs["sample_id"] == sid]

    s = scr.Scrublet(sub.X, expected_doublet_rate=0.06)
    scores, preds = s.scrub_doublets()

    # 写回全局对象
    adata.obs.loc[sub.obs_names, "doublet_score"] = scores
    adata.obs.loc[sub.obs_names, "predicted_doublet"] = preds

# 过滤双细胞
adata = adata[~adata.obs["predicted_doublet"]].copy()
```

输出：`02_qc_filtered.h5ad`

---

## **2.5 添加组织注释（metadata 映射）**

示例：

```python
import pandas as pd

meta = pd.read_csv("metadata/SraRunTable.csv")
map_dict = dict(zip(meta["Run"], meta["tissue"]))

adata.obs["tissue"] = adata.obs["sample_id"].map(map_dict)
```

输出：`02_qc_filtered_annotated.h5ad`

---

## **2.6 组织构成统计（按组织计数）**

```python
tissue_counts = (
    adata.obs
    .groupby("tissue")
    .size()
    .sort_values(ascending=False)
)
```

可输出 CSV 或绘制条形图。

---

# **3. 高变基因选择（Per-Tissue HVG Selection）**

目标：消除组织数量差异造成的偏差。

详细流程：

```python
import numpy as np

all_hvg = set()

for tissue in adata.obs["tissue"].unique():
    sub = adata[adata.obs["tissue"] == tissue].copy()

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    sc.pp.highly_variable_genes(sub, n_top_genes=2000)
    hvg = sub.var[sub.var["highly_variable"]].index

    all_hvg |= set(hvg)

# 截取 HVG 并合并
adata = adata[:, list(all_hvg)].copy()

# 全局标准化
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

输出：`03_preprocessed.h5ad`

---

# **4. BBKNN 整合 + UMAP + Leiden 聚类**

## **4.1 PCA（消除 UMI、线粒体影响）**

```python
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
```

---

## **4.2 BBKNN 整合**

```python
import scanpy.external as sce

sce.pp.bbknn(
    adata,
    batch_key="sample_id",
    n_pcs=50,
    neighbors_within_batch=3
)
```

---

## **4.3 UMAP 与 Leiden**

```python
sc.tl.umap(
    adata,
    min_dist=0.25,
    spread=1.2
)

sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="leiden"
)
```

输出：`04_bbknn_full.h5ad`

---

# **5. 细胞大类注释（Marker-driven Cluster Annotation）**

示例流程：

```python
markers = pd.read_csv("data/markers.csv")
type_to_genes = markers.groupby("large_type")["marker_gene"].apply(list)
cluster_to_type = {}

for cl in adata.obs["leiden"].unique():
    sub = adata[adata.obs["leiden"] == cl]

    scores = {
        t: sub[:, genes].X.sum()
        for t, genes in type_to_genes.items()
        if all(g in adata.var_names for g in genes)
    }

    best_type = max(scores, key=scores.get)
    cluster_to_type[cl] = best_type

adata.obs["large_type"] = adata.obs["leiden"].map(cluster_to_type)
```

---

# **6. DOCK9 表达分析（组织层面 + UMAP 空间 + 细胞类型）**

## **6.1 raw（QC 层面）按组织统计 DOCK9**

```python
expr = adata[:, "DOCK9"].X.toarray().flatten()

adata.obs["DOCK9"] = expr
adata.obs["DOCK9_binary"] = (expr > 0).astype(int)

summary = (
    adata.obs
    .groupby("tissue", observed=False)
    .agg(
        mean_expr=("DOCK9", "mean"),
        median_expr=("DOCK9", "median"),
        detect_fraction=("DOCK9_binary", "mean"),
    )
)
```

---

## **6.2 DOCK9 在整合 UMAP 空间中的表达**

```python
adata.obs["DOCK9_expr"] = expr_df.loc[adata.obs_names, "DOCK9_expr"]

sc.pl.umap(adata, color="DOCK9_expr", cmap="viridis")
sc.pl.umap(adata, color="DOCK9_binary", cmap="Set1")
```

---

## **6.3 按大类细胞类型统计 DOCK9**

```python
dock9_stats = (
    adata.obs
    .groupby("large_type")
    .agg(
        mean_expr=("DOCK9_expr", "mean"),
        detect_fraction=("DOCK9_binary", "mean")
    )
)
```

---


