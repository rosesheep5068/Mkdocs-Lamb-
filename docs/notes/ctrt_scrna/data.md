# CTCR Project — Data Processing Notes

---

# **1. 质量控制（Quality Control, QC）**

在对 **320 万细胞**进行后续整合前，需要对 IGT / PanCancer 原始 h5ad 执行 QC，包括计算质控指标、过滤低质量细胞与保存干净的对象。

## **1.1 加载原始数据并计算 QC 指标**

```python
import scanpy as sc

# 载入原始 count 矩阵（含 raw counts）
adata = sc.read_h5ad("igt_s9_fine_counts.h5ad")

# 计算 per-cell QC 指标
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["MT"],   # 若已在 var 中标记线粒体基因
    percent_top=None, # 可加入 n_top_counts
    inplace=True
)
```

得到的关键 QC 指标：

- `total_counts` —— 总 UMI

- `n_genes_by_counts` —— 检测到的基因数

- `pct_counts_mt` —— 线粒体占比

这些指标用于后续过滤与回归校正。

## **1.2 过滤低质量细胞**

```python
# 过滤高线粒体比例细胞
adata = adata[adata.obs["pct_counts_mt"] < 20, :]

# 过滤低基因数量细胞
adata = adata[adata.obs["n_genes_by_counts"] > 200, :]

# （可选）过滤极端高 UMI 或 potential doublets
adata = adata[adata.obs["total_counts"] < 5e5, :]
```

最终保存为：  
`qc_igt.h5ad`、`qc_pancancer.h5ad`

---

# **2. 归一化（Normalize）与高变基因筛选（HVG）**

QC 后的数据执行 NormalizeTotal → Log1p → HVG 筛选，为降维做准备。

## **2.1 NormalizeTotal + Log1p**

```python
# 每细胞归一化到 1e4 counts（Cell Ranger 默认）
sc.pp.normalize_total(
    adata,
    target_sum=1e4
)

# log1p transform 稳定表达分布
sc.pp.log1p(adata)
```

该步骤消除文库大小差异，使不同细胞在表达层面可比较。

## **2.2 筛选高变基因（HVG）**

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,     # 经典数量
    flavor="cell_ranger", # 大规模数据常用
    subset=False          # 不自动 subset，手动 subset 更安全
)

# 仅保留 HVG 用于降维
adata = adata[:, adata.var["highly_variable"]].copy()
```

输出：  
`norm_hvg_igt.h5ad`、`norm_hvg_pancancer.h5ad`

---

# **3. 未批次校正的降维与聚类（Neighbors）**

这一步用于**观察原始批次效应是否明显**，同时作为对照与 BBKNN 比较。

## **3.1 回归技术变量**

```python
# 回归掉文库大小与线粒体比例，减少技术噪声
sc.pp.regress_out(
    adata,
    keys=["total_counts", "pct_counts_mt"]
)

# 标准化，使所有 HVG 的方差一致
sc.pp.scale(adata, max_value=10)
```

## **3.2 PCA + Neighbors + UMAP**

```python
# PCA 降维
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# 构建 KNN 图（未批次校正）
sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=30
)

# 低维可视化
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
```

## **3.3 Leiden 聚类**

```python
sc.tl.leiden(
    adata,
    resolution=0.5,    # 可视数据分辨率
    flavor="igraph"
)
```

输出：  
`preproc_igt_neighbors.h5ad`

---

# **4. BBKNN 批次校正与整合（适用于百万级细胞）**

由于全量数据 **>3.2M 细胞**，直接运行 BBKNN 会爆内存；因此采用“下采样 + ingest full”的高效策略。

## **4.1 进行下采样（如 100 万细胞）**

```python
import numpy as np

# 假设 IGT 全量为 3.2M，这里取 1M 支撑 BBKNN
keep_idx = np.random.choice(
    adata.n_obs,
    size=1_000_000,
    replace=False
)
adata_sub = adata[keep_idx].copy()
```

下采样必须保持**不同组织的细胞类型代表性**。

## **4.2 回归变量 → PCA → BBKNN 整合**

```python
# 回归 & scale
sc.pp.regress_out(adata_sub, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata_sub, max_value=10)

# PCA
sc.tl.pca(adata_sub, n_comps=50)

# BBKNN 批次整合（batch_key 对应 sampleID）
import scanpy.external as sce

sce.pp.bbknn(
    adata_sub,
    batch_key="sampleID",
    n_pcs=50,
    neighbors_within_batch=3,  # 控制跨批次混合程度
    use_faiss=True             # GPU 优化（若可用）
)
```

## **4.3 UMAP + Leiden 聚类**

```python
sc.tl.umap(adata_sub, min_dist=0.25, spread=1.2)

sc.tl.leiden(
    adata_sub,
    resolution=0.4,
    flavor="igraph"
)
```

输出：  
`preproc_igt_bbknn_subsample.h5ad`

## **4.4 将全量数据 ingest 到 subsample embedding 中**

```python
adata_full = adata.copy()

# 与 subsample 保持相同的 PCA 设置
sc.pp.regress_out(adata_full, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata_full, max_value=10)
sc.tl.pca(adata_full, n_comps=50)

# ingested UMAP（映射全量细胞）
sc.tl.ingest(adata_full, adata_sub, obs="leiden")
```

输出：  
`preproc_igt_bbknn_full.h5ad`

---

# **5. BBKNN / Neighbors 可视化与 QC 评估**

通过比较两个 embedding 的 UMAP，可以判断：

- 原始批次效应是否显著

- BBKNN 是否成功整合组织差异

- 细胞是否按预期的生物类别聚集

示例：

```python
# 查看批次是否被纠正
sc.pl.umap(adata_neighbors, color="sampleID")

# 查看 QC 是否影响结构
sc.pl.umap(adata_bbknn, color=["total_counts", "pct_counts_mt"])

# 查看 cluster 结构差异
sc.pl.umap(adata_bbknn, color="leiden", legend_loc="right margin")
```

输出到：  
`figs_neighbors/`、`figs_bbknn/`



---

# **6. DOCK9 在不同 embedding 中的表达对比**

为了确认 DOCK9 是否受 **批次效应（batch effect）** 影响，我们将其表达分别映射到：

- **Neighbors（未校正）embedding**

- **BBKNN（批次校正）embedding**

并观察表达模式是否一致。

```python
# Neighbors 版本（未批次校正）
sc.pl.umap(
    adata_neighbors,
    color="DOCK9",
    cmap="viridis",
    frameon=False,
    title="DOCK9 expression (Neighbors)",
)

# BBKNN 版本（批次校正）
sc.pl.umap(
    adata_bbknn,
    color="DOCK9",
    cmap="viridis",
    frameon=False,
    title="DOCK9 expression (BBKNN)",
)
```

---

# **7. DOCK9 表达提取、映射与 cluster 统计**

用于后续：

- TF × celltype 相关性分析

- 生成连续 / 二值化 DOCK9 图

- 生成 DOCK9 每 cluster 的检测比例

## **7.1 从 QC 数据中提取 DOCK9 表达**

```python
expr = adata_qc[:, "DOCK9"].X
if hasattr(expr, "toarray"):  # 稀疏矩阵处理
    expr = expr.toarray().flatten()
```

## **7.2 将表达映射回 neighbors / BBKNN**

```python
adata.obs["DOCK9_expr"] = expr
adata.obs["DOCK9_binary"] = (expr > 0).astype(int)
```

## **7.3 按 cluster 统计 DOCK9 表达与检出率**

```python
cluster_stats = (
    adata.obs.groupby("leiden")["DOCK9_expr"]
    .agg(
        mean="mean",
        median="median",
        detect_fraction=lambda x: (x > 0).mean(),
        n_cells="count"
    )
)
```

输出：

- per-cell DOCK9（csv）

- per-cluster DOCK9 summary（csv）

- 连续/二值化 DOCK9 表达图

---

# **8. 按组织（tissue / system）绘制 UMAP**

用于检查跨组织拓扑结构是否合理，并用于筛选 DOCK9 相关的细胞大类（上皮、间充质、内皮）。

```python
sc.pl.umap(
    adata,
    color="tissue",
    legend_loc="right margin",
    frameon=False,
)

sc.pl.umap(
    adata,
    color="system",
    legend_loc="right margin",
    frameon=False,
)
```

---

# **9. FOXC2 表达提取与分析（TF 候选验证）**

方法与 DOCK9 相同，只是将基因换成 FOXC2，用于验证其作为候选 TF 的可能性。

```python
expr = adata_qc[:, "FOXC2"].X
adata.obs["FOXC2_expr"] = np.array(expr).flatten()
adata.obs["FOXC2_binary"] = (adata.obs["FOXC2_expr"] > 0).astype(int)
```

FOXC2 的结果用于后续生成 **TF × celltype** 的平均表达矩阵（对应 meanu1.png）。

---

# **10. 细胞类型注释（八大类）**

对 Neighbors + BBKNN embedding 进行八大类注释（B cell、T cell、Myeloid、Epithelial、Stromal、Endothelial 等）。

## **10.1 构建 marker 字典**

```python
df = pd.read_csv("canonical_marker_8types.csv")
marker_dict = {
    cat: group["Marker"].unique().tolist()
    for cat, group in df.groupby("Cell_category")
}
```

## **10.2 逐细胞计算 marker score 并预测 celltype**

```python
scores = {
    cat: adata[:, [g for g in genes if g in adata.var_names]].X.mean(axis=1)
    for cat, genes in marker_dict.items()
}

score_df = pd.DataFrame(scores, index=adata.obs_names)
adata.obs["celltype_pred"] = score_df.idxmax(axis=1)
```

## **10.3 cluster-level consensus**

```python
adata.obs["celltype_majority"] = (
    adata.obs.groupby("leiden")["celltype_pred"]
    .transform(lambda x: x.mode()[0])
)
```

## **10.4 绘制注释后的 UMAP**

```python
sc.pl.umap(adata, color="celltype_pred", legend_loc="right margin")
sc.pl.umap(adata, color="celltype_majority", legend_loc="right margin")
```

输出：

- `annotated_neighbors.h5ad`

- `annotated_bbknn_subsample.h5ad`

---

# **11. 注释 UMAP 批量绘制**

为每个标签生成独立注释图（celltype_pred / majority）。

```python
for label in ["celltype_pred", "celltype_majority"]:
    sc.pl.umap(
        adata,
        color=label,
        legend_loc="right margin",
        frameon=False
    )
```

输出文件：  
`umap_annotation_neighbors_*`、`umap_annotation_bbknn_subsample_*`

---

# **12. 基于 marker 重新筛选三类关键细胞**

根据 DOCK9 表达分布特征与项目重点，我们筛选：  
✔ 上皮 epithelial  
✔ 间充质 stromal  
✔ 内皮 endothelial

## **12.1 per-cell marker scoring**

```python
adata, score_df = per_cell_marker_score(adata, marker_dict)
```

## **12.2 筛选三大类**

```python
target = ["stromal", "epithelial", "endothelial"]
adata_sub = adata[adata.obs["marker_pred"].isin(target)]
```

剩余约 **110 万细胞**。

输出：  
`norm_hvg_igt_stromal_epi_endo_marker.h5ad`

---

# **13. 三类细胞的 BBKNN 整合**

用于生成最终 U2U1 / EX2U1 图。

## **13.1 内存优化 + 再下采样（如 >150 万）**

```python
if adata.n_obs > MAX_CELLS:
    keep_idx = np.random.choice(adata.n_obs, MAX_CELLS, replace=False)
    adata_sub = adata[keep_idx].copy()
```

## **13.2 回归 + PCA + BBKNN**

```python
sc.pp.regress_out(adata_sub, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata_sub)
sc.tl.pca(adata_sub, n_comps=50)

sce.pp.bbknn(
    adata_sub,
    batch_key="sampleID",
    n_pcs=50,
    neighbors_within_batch=3,
)
```

## **13.3 再聚类**

```python
sc.tl.umap(adata_sub)
sc.tl.leiden(adata_sub, resolution=0.5)
```

输出：

- `preproc_stromal_epi_endo_bbknn_subsample.h5ad`

- `preproc_stromal_epi_endo_bbknn_full.h5ad`

---

# **14. 三类细胞整合结果的 UMAP 可视化**

绘制 Leiden / sampleID / QC 指标。

```python
sc.pl.umap(adata, color="leiden")
sc.pl.umap(adata, color="sampleID")
sc.pl.umap(adata, color="pct_counts_mt")
```

输出：  
`figs_stromal_epi_endo_bbknn/`

---

# **15. DOCK9 Topic Analysis（表达 + TF 面板 + UMAP）**

## **15.1 从 QC 读取 DOCK9 + TF 表达**

```python
expr_df = adata_qc[:, GENE_PANEL].to_df()
expr_df["cellID"] = adata_qc.obs_names
expr_df.to_csv("TF_expression_matrix.csv")
```

GENE_PANEL 包括：

- DOCK9

- FOXC2、PBX1/2、MEIS1/2、ONECUT3、HNF1A/B … 共 15 个 TF

## **15.2 映射表达到 subset**

```python
adata_sub.obs["DOCK9"] = expr_df.loc[adata_sub.obs_names, "DOCK9"]
adata_sub.obs["DOCK9_binary"] = (adata_sub.obs["DOCK9"] > 0)
```

## **15.3 cluster 层面统计**

```python
summary = (
    adata_sub.obs.groupby("leiden")["DOCK9"]
    .agg(["mean", "median", lambda x: (x>0).mean(), "count"])
)
```

## **15.4 绘图（UMAP / Violin / DotPlot）**

包括：

- DOCK9 连续表达图

- DOCK9 detected 二值图

- cluster-level violin

- celltype × DOCK9 dotplot

---

# **16. 三大类细胞的精细化注释**

基于更细的 marker（如腺体、基底、成纤维细胞、毛细血管等）。

## **16.1 加载专用 marker panel**

```python
marker_dict = load_markers("marker_stromal_epi_endo.csv")
```

## **16.2 per-cell 注释**

```python
adata.obs["celltype_pred"] = score_df.idxmax(axis=1)
```

## **16.3 cluster 共识注释**

```python
adata.obs["celltype_majority"] = adata.obs["leiden"].map(cluster_label_dict)
```

## **16.4 输出与可视化**

- `annotated_stromal_epi_endo.h5ad`

- `umap_annotation_stromal_epi_endo_celltype_pred.png`

- `umap_annotation_stromal_epi_endo_celltype_majority.png`

---

# **17. 构建 TF × cell metadata 表（用于 meanu1.png）**

```python
merged_df = pd.merge(meta_df, tf_df, on="cellID")
merged_df.to_csv("TF_expression_percell.csv")
```

用于：

- TF 平均表达矩阵

- TF × celltype 相关性图

---

# **18. 样本数量统计（3.2M 细胞结构）**

```python
summary = {
    "Total_cells": adata.n_obs,
    "Unique_samples": obs['sampleID'].nunique(),
    "Unique_donors": obs['donorID'].nunique(),
    "Unique_datasets": obs['datasetID'].nunique(),
}
```

输出：

- sample_summary.csv

- sample_cell_counts.csv

---

# **19. DOCK9 按组织层面表达总结**

```python
expr = adata[:, "DOCK9"].X
adata.obs["DOCK9"] = np.array(expr).flatten()
```

按组织统计：

```python
summary = adata.obs.groupby("tissue").agg(...)
```

绘制条形图作为 Supplementary。

---

# **20. 按细胞类型总结 DOCK9 表达**

## **20.1 合并注释 + DOCK9 表达**

```python
meta = adata.obs.merge(dock9_df, on="cellID", how="inner")
```

## **20.2 按 celltype 聚合**

```python
df = meta.groupby("celltype").agg(
    mean_expr=("DOCK9", "mean"),
    pct_detected=("DOCK9", lambda x: (x>0).mean()),
    n_cells=("DOCK9", "count")
)
```

## **20.3 柱状图可视化**

- mean DOCK9 expression

- detection fraction

用于展示不同细胞类型的 DOCK9 功能偏好性。

---
