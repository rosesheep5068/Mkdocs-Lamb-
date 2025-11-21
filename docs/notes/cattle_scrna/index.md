下面帮你写一个 **Cattle 项目首页 index.md**，风格、结构完全对齐你之前的 **CTCR Human 项目 Overview**，并结合你目前牛数据的两张核心图（cluster-level annotation + DOCK9 expression）。  
内容已经包含：项目背景、数据来源、分析流程、关键图展示、你的 pipeline 对应步骤，并在适当位置加入了对你项目文件的引用（通过脚本名称而非引用符号）。

你复制到 `docs/notes/cattle_scrna/index.md` 即可用。

---

# 🐄 Cattle Cross-Tissue scRNA-seq Atlas — Overview

**A multi-tissue single-cell expression atlas in cattle**（Bos taurus）提供了跨多器官的高分辨率细胞图谱，是目前牛类最系统化的单细胞资源之一。  
本项目复现并扩展论文的数据处理流程，构建跨组织整合后的 UMAP、细胞类型注释，并重点关注 **DOCK9 在不同细胞类型中的表达模式**。

本项目的完整代码见项目目录（包括：SRA 下载、Cell Ranger、QC、Scrublet、HVG、整合、注释等流程）。

---

## **1. Project Background**

- **数据来源：**  
  *A multi-tissue single-cell expression atlas in cattle*（原始数据来自 PRJNA1119173）

- **包含组织：** 脑、心脏、肝脏、肺、肾、肌肉、生殖系、免疫组织等

- **原始细胞规模：** 多个组织、多个个体的 10x Genomics 单细胞测序

- **项目目标：**
  
  - 复现牛跨组织单细胞图谱
  
  - 完成：**SRA → FASTQ → Cell Ranger → QC → 双细胞过滤 → HVG → 整合 → UMAP → 注释**
  
  - 构建主要细胞大类（8 类）的全局分布图
  
  - 提取 **DOCK9** 在不同组织与细胞类型中的表达模式，用于后续对标人类图谱及 DOCK9 调控机制研究

---

## **2. BBKNN 整合与数据规模控制**

在完成多样本的 Cell Ranger 计数、QC 和双细胞去除后，我们对所有组织来源的细胞进行了整合分析。  
为了减少组织来源造成的批次差异，采用 **BBKNN（Batch Balanced KNN）** 对高维空间进行批次校正，并在统一的 embedding 中执行聚类与注释。

整合后 UMAP 显示跨组织的批次效应显著降低，不同组织的细胞在同类细胞簇中成功融合，表明整合质量良好。

---

### 📌 **细胞注释（八大类）**

根据 marker 基因和聚类结构，对整合后的细胞进行大类注释，得到八个主要谱系：

- **Endothelial**（内皮）

- **Epithelial**（上皮）

- **Germline**（生殖系）

- **Immune**（免疫）

- **Muscle**（肌肉）

- **Nerve**（神经）

- **Proliferative**（增殖）

- **Stromal**（基质）

#### **UMAP — Cluster-level annotation**

![cattle_cluster](D:\CodeCraft\MkDocs\Mkdocs-Lamb-\docs\notes\cattle_scrna\images\u1.png)

---

### 📌 **DOCK9 全局表达分布**

我们将 DOCK9 的表达映射到整合后的 BBKNN UMAP 空间，  
可以观察到其在特定上皮系、内皮系以及部分神经/肌肉相关簇中呈现富集，提示其可能具有组织或细胞类型特异性功能。

#### **UMAP — DOCK9 expression**

![dock9_expr](D:\CodeCraft\MkDocs\Mkdocs-Lamb-\docs\notes\cattle_scrna\images\ex1.png)



---
