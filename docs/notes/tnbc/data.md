# TNBC Data Processing

---

# 1. 数据来源（TCGA + METABRIC）

本项目使用：

- **TCGA-BRCA**：用于差异分析 & 表达量比较

- **METABRIC**：用于生存分析验证

---

# 2. TCGA 数据获取与准备

## **2.1 下载 TCGA-BRCA STAR-counts 数据**

```r
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)
```

## **2.2 提取表达矩阵与样本信息**

```r
counts_data <- assay(data)
counts_data <- log2(counts_data + 1)
sample_info <- colData(data)
```

## **2.3 提取 TNBC（Basal）与 Normal（No）样本**

```r
t_needed = c("vital_status", "paper_BRCA_Subtype_PAM50")
meta = sample_info[, t_needed]
meta$paper_BRCA_Subtype_PAM50[is.na(meta$paper_BRCA_Subtype_PAM50)] <- "No"

basal_normal_samples <- meta$paper_BRCA_Subtype_PAM50 %in% c("Basal", "No")
counts_basal_normal <- counts_data[, basal_normal_samples]
meta_basal_normal <- meta[basal_normal_samples, ]
```

---

# 3. 基因注释与矩阵清洗

## **3.1 Ensembl → Gene Symbol 映射**

```r
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(counts_basal_normal)

mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = sub("\\..*", "", ensembl_ids),
                 mart = ensembl)
```

## **3.2 替换矩阵行名**

```r
rownames(counts_basal_normal) <- mapping$external_gene_name[
  match(sub("\\..*", "", rownames(counts_basal_normal)), mapping$ensembl_gene_id)
]
```

---

# 4. 差异表达分析（Basal vs Normal）

## **4.1 过滤低表达基因**

```r
expr_matrix <- as.matrix(counts_basal_normal)
keep <- rowMeans(expr_matrix > 1) > 0.5
expr_matrix <- expr_matrix[keep, ]
```

## **4.2 limma 差异分析**

```r
library(limma)
group <- factor(meta_basal_normal$paper_BRCA_Subtype_PAM50, levels = c("No", "Basal"))
design <- model.matrix(~ group)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "groupBasal", number = Inf)
```

---

# 5. 铁代谢基因交集筛选

```r
given_genes <- c("ABCB6","GLRX3","SLC46A1", ..., "FTMT","ATP6V1G1")
significant_gene_names <- rownames(subset(results, adj.P.Val < 0.05 & abs(logFC) > 1))
common_genes <- intersect(given_genes, significant_gene_names)
```

---

# 6. 关键基因表达提取与可视化

## **6.1 提取三个核心基因表达（FTH1 / NCOA4 / SUV39H2）**

```r
selected <- c("FTH1", "NCOA4", "SUV39H2")
exp_data <- counts_data_renamed[rownames(counts_data_renamed) %in% selected, ]
```

## **6.2 转换为长格式**

```r
exp_data_long <- exp_data %>% as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = ifelse(Sample %in% TNBC_samples, "TNBC", "Normal"))
```

## **6.3 小提琴图示例（以 NCOA4 为例）**

```r
ggplot(exp_data_long %>% filter(Gene == "NCOA4"),
       aes(Group, Expression, colour = Group)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.5) +
  geom_jitter(alpha = 0.3) +
  geom_signif(comparisons = list(c("TNBC","Normal")),
              test = "wilcox.test", map_signif_level = TRUE) +
  theme_bw()
```

---

# 7. 基因相关性分析（TNBC 样本）

```r
SN_exp_data <- counts_data_renamed[c("SUV39H2","NCOA4"), TNBC_samples]
cor_value <- cor(SN_exp_data[1,], SN_exp_data[2,], method = "pearson")

ggplot(data.frame(SUV39H2 = SN_exp_data[1,], NCOA4 = SN_exp_data[2,]),
       aes(SUV39H2, NCOA4)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson")
```

---

# 8. METABRIC 生存分析（验证）

## **8.1 TNBC 样本筛选**

```r
cli_data <- merge(cli_sample_data, cli_patient_data, by = "PATIENT_ID")
cli_info <- subset(cli_data, ER_STATUS=="Negative" & HER2_STATUS=="Negative" & PR_STATUS=="Negative")
```

## **8.2 表达矩阵整理与合并**

```r
exp_info <- t(exp_data)
colnames(exp_info) <- exp_info[1, ]
exp_info <- exp_info[-c(1,2), ]
meta_data <- merge(exp_info, cli_info, by = "SAMPLE_ID")
meta_data$OS_YEARS <- meta_data$OS_MONTHS / 12
```

---

# 9. 单基因生存分析（NCOA4 示例）

## **9.1 中位数分组**

```r
NCOA4_merged <- meta_data[, c("SAMPLE_ID","NCOA4","OS_YEARS","OS_STATUS")]
NCOA4_merged$group <- ifelse(NCOA4_merged$NCOA4 > median(NCOA4_merged$NCOA4), "High","Low")
fit <- survfit(Surv(OS_YEARS, OS_STATUS) ~ group, data = NCOA4_merged)
```

## **9.2 绘制 Kaplan–Meier 曲线**

```r
ggsurvplot(fit, data = NCOA4_merged,
           pval = TRUE, conf.int = TRUE,
           risk.table = "absolute",
           surv.median.line = "hv")
```

---

# 10. 双基因组合生存分析（示例：SUV39H2 × NCOA4）

## **10.1 基于 SUV39H2 中位数分组**

```r
SN_data$SUV39H2_group <- ifelse(SN_data$SUV39H2 >= median(SN_data$SUV39H2), "high", "low")
```

## **10.2 对每个子组再次按 NCOA4 中位数分组**

```r
median_NCOA4_high <- median(SN_data$NOCA4[SN_data$SUV39H2_group=="high"])
SN_data$NCOA4_group_high <- ifelse(SN_data$NOCA4 >= median_NCOA4_high, "high","low")
```

## **10.3 KM 曲线绘制**

```r
fit2 <- survfit(Surv(time, status) ~ group, data = SN_data)

ggsurvplot(fit2, data = SN_data,
           pval = TRUE, conf.int = TRUE,
           risk.table = "absolute",
           legend.title = "Strata")
```

---


