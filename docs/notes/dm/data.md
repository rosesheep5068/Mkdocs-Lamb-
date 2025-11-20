# DM Proteomics Data Processing

---

# 1. 数据来源（Label-free Proteomics）

本项目使用：

- **Ctrl（对照）**

- **UVB（紫外损伤）**

- **DM（处理组）**

并基于标准化后的蛋白表达矩阵进行分析。

---

# 2. 数据导入与准备

## **2.1 读取差异分析结果与表达矩阵**

```r
UVB_ctrl_data <- read.csv("UVBvsctrl_deg_all.csv")
DM_UVB_data   <- read.csv("DMvsUVB_deg_all.csv")
all_exp       <- read.csv("all_compare.csv")
```

---

# 3. 差异蛋白筛选（UVB vs Ctrl & DM vs UVB）

## **3.1 设置差异筛选阈值（log2FC + p-value）**

```r
UVB_ctrl <- UVB_ctrl_data %>%
  mutate(regulate = case_when(
    log2FoldChange > 0.263 & pvalue < 0.05 ~ "Up",
    log2FoldChange < -0.263 & pvalue < 0.05 ~ "Down",
    TRUE ~ "Normal"
  ))

DM_UVB <- DM_UVB_data %>%
  mutate(regulate = case_when(
    log2FoldChange > 0.263 & pvalue < 0.05 ~ "Up",
    log2FoldChange < -0.263 & pvalue < 0.05 ~ "Down",
    TRUE ~ "Normal"
  ))
```

## **3.2 提取差异蛋白并保存**

```r
UVB_ctrl_diff <- UVB_ctrl[UVB_ctrl$regulate %in% c("Up", "Down"), ]
DM_UVB_diff   <- DM_UVB[DM_UVB$regulate %in% c("Up", "Down"), ]

write.csv(UVB_ctrl_diff, "result/DM_result3.0/UVB_ctrl.csv", row.names = FALSE)
write.csv(DM_UVB_diff, "result/DM_result3.0/DM_UVB.csv", row.names = FALSE)
```

---

# 4. 差异蛋白集合交集分析

## **4.1 提取基因名集合**

```r
UVB_ctrl_diff_genes <- UVB_ctrl_diff$gene_name
DM_UVB_diff_genes   <- DM_UVB_diff$gene_name
```

## **4.2 绘制 Venn 图（集合交集展示）**

```r
gene_list <- list(
  UVB_ctrl_diff_genes = UVB_ctrl_diff_genes,
  DM_UVB_diff_genes   = DM_UVB_diff_genes
)

ggvenn(gene_list,
       fill_color = c("#4682B4", "#32CD32"),
       show_percentage = FALSE,
       set_name_size = 15,
       text_size = 10)
```

---

# 5. 目标基因集比对（Intersection with CSIG genes）

## **5.1 定义目标基因集**

```r
CSIG_genes <- c("PAI-1","MMP3","TNF-alpha","HMGB1", ... )
```

## **5.2 找出交集基因**

```r
UVB_ctrl_target_overlap <- intersect(UVB_ctrl_diff_genes, CSIG_genes)
DM_UVB_target_overlap   <- intersect(DM_UVB_diff_genes, CSIG_genes)
common_target_overlap   <- intersect(UVB_ctrl_target_overlap, DM_UVB_target_overlap)
```

---

# 6. 提取目标基因表达并绘制热图

## **6.1 从表达矩阵提取目标基因**

```r
meta_data <- read.csv("all_compare.csv")
target_genes <- common_target_overlap
subset_data <- meta_data[meta_data$gene_name %in% target_genes, ]
```

## **6.2 Z-score 标准化 + 分组注释**

```r
fpkm_data <- subset_data[, c("con1_fpkm","con2_fpkm","con3_fpkm",
                             "UVB1_fpkm","UVB2_fpkm","UVB3_fpkm",
                             "DM1_fpkm","DM2_fpkm","DM3_fpkm")]

fpkm_data_zscore <- t(scale(t(fpkm_data)))

annotation_col <- data.frame(group = rep(c("Ctrl","UVB","DM"), each = 3))
rownames(annotation_col) <- colnames(fpkm_data_zscore)
```

## **6.3 绘制热图**

```r
pheatmap(
  fpkm_data_zscore,
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE
)
```

---

# 7. GO 富集分析（BP / CC / MF）

## **7.1 进行 ID 转换**

```r
combined_go <- bind_rows(UVB_ctrl_diff, DM_UVB_diff)
combined_go$gene_id <- gsub("_\\d+$", "", combined_go$gene_id)

combined_go_entrez <- bitr(combined_go$gene_id,
                           fromType = "ENSEMBL",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)
```

## **7.2 三类 GO（BP/CC/MF）富集**

```r
combined_go_BP <- enrichGO(combined_go_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP")
combined_go_CC <- enrichGO(combined_go_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC")
combined_go_MF <- enrichGO(combined_go_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF")
```

## **7.3 合并并绘制 GO Top20 气泡图**

```r
combined_GO <- rbind(combined_go_BP@result,
                     combined_go_CC@result,
                     combined_go_MF@result)

top_GO <- combined_GO %>% arrange(qvalue) %>% head(20)

ggplot(top_GO, aes(RichFactor, Description, size = Count, color = qvalue)) +
  geom_point()
```

---

# 8. GO 分类柱状图（BP / CC / MF）

```r
dt <- read.csv("top_go.csv")
df <- dt %>% arrange(Class, desc(number_of_gene))

ggplot(df, aes(x = GO_Term, y = number_of_gene, fill = Class)) +
  geom_col() +
  theme_classic()
```

---

# 9. KEGG 富集分析（Top20 Pathways）

## **9.1 KEGG 通路富集**

```r
entrez_genes <- mapIds(org.Hs.eg.db, keys = combined_go$gene_id,
                       column = "ENTREZID", keytype = "ENSEMBL")
entrez_genes <- na.omit(entrez_genes)

kegg_enrichment <- enrichKEGG(gene = entrez_genes, organism = "hsa")
```

## **9.2 绘制 KEGG 气泡图**

```r
kegg_data <- as.data.frame(kegg_enrichment) %>% head(20)

ggplot(kegg_data, aes(RichFactor, Description,
                      size = GeneRatio_value, color = pvalue)) +
  geom_point(alpha = 0.7)
```

---

# 10. GSEA 分析（GO / KEGG / Reactome）

## **10.1 构建基因排序表**

```r
genelist <- DM_UVB_gsea$log2FoldChange
names(genelist) <- DM_UVB_gsea$entrez_id
genelist <- sort(genelist, decreasing = TRUE)
```

## **10.2 GSEA 分析**

```r
Go_gseresult    <- gseGO(genelist, OrgDb = org.Hs.eg.db)
KEGG_gseresult  <- gseKEGG(geneList = genelist, organism = "hsa")
Reactome_result <- gsePathway(geneList = genelist, organism = "human")
```

---

# 11. Reactome 通路富集与可视化

```r
reactome_result <- enrichPathway(gene = reactome_genes$ENTREZID,
                                 organism = "human")

dotplot(reactome_result)
```

---

# 12. 特殊筛选（示例：某单个通路前后 ±10 条）

项目中包含特定展示需求（如某条 GO/KEGG 通路前后随机选取 20 条），可在下方示例代码中实现：

```r
go_to_replace <- combined_GO[combined_GO$ID == "GO:0090398", ]
```

（此部分根据项目需求使用）

---
