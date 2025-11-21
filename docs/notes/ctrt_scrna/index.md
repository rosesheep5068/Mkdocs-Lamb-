# CTCR Project — Overview

Human Multi-Organ scRNA-seq Atlas (3.2M cells)

---

## 1. Project Background

数据来源Cross-tissue multicellular coordination and its rewiring in cancer 一文，为一个跨多组织的大规模人类单细胞转录组数据集，总计约 **320 万个细胞**。  
项目目标包括：

- 构建跨器官的高分辨率细胞图谱  
- 完成 QC → 标准化 → 批次整合 → 聚类 → 注释 的完整流程  
- 结合 DOCK9 的表达特征进一步筛选相关细胞群  
- 针对上皮、间充质、内皮细胞进行深入亚群分析  
- 探索 DOCK9 及其候选调控因子（FOXC2、PBX1/2、MEIS1/2 等）在不同细胞类型中的表达模式与潜在调控关系  

---

## 2. BBKNN 整合与数据规模控制

CTCR 全量包含 **3.2M 细胞**，在执行 BBKNN 时遇到显著内存瓶颈。  
我们通过**多次向下采样**进行可行性测试，最终选择：

- **100 万细胞规模（1M）** 作为 BBKNN 批次校正的最佳折中点  
- 在保证不同组织细胞类型完整性的前提下，实现计算可行性  
- 得到高质量的跨组织 UMAP 与聚类结构  

### 📌 细胞注释（八大类）

![u1](./images/u1.png)

### 📌 DOCK9 全局表达（3.2M → 1M）

![ex1u1](./images/ex1.png)

-----

## 3. 细胞类型注释与 DOCK9 相关细胞群筛选

在 QC → 标准化（Normalize + HVG）→ 降维 → 聚类后，我们基于 **八大类 canonical markers** 完成初步注释，包括：

- Epithelial  
- Stromal / Mesenchymal  
- Endothelial  
- Myeloid  
- B cells  
- T cells (CD4/CD8)  
- Others  

由于 DOCK9 在多数免疫细胞中表达极低，我们进一步**根据 DOCK9 的低表达与空间分布特征对细胞类型进行筛选**，最终选择：

### ✔ 上皮细胞（Epithelial）

### ✔ 间充质细胞（Stromal / Mesenchymal）

### ✔ 内皮细胞（Endothelial）

三类共约 **110 万细胞（1.1M）** 作为后续深度整合与分析对象。

---

## 4. 上皮、间质、内皮细胞的精细亚群注释

我们将 1.1M 细胞重新进行 BBKNN 批次整合、聚类并建立新的 UMAP embedding，对细胞亚群进行细化注释，包括：

- Arterial / Venous / Capillary endothelial  
- Basal / Luminal / Ciliated epithelial  
- Fibroblast, Myofibroblast  
- Goblet cell, Pericyte, RPE-like cells 等  

### 📌 选定三大类细胞（Epi / Stromal / Endo）的重新注释

![u2u1](./images/u2.png)

### 📌 DOCK9 在精细亚群中的表达分布

![ex2u1](./images/ex2.png)

---

## 5. DOCK9–TF–Celltype 关联分析

基于 rs1408718 及其连锁 SNP，我们使用 **JASPAR TFBS Extraction Tool** 预测了潜在调控 DOCK9 的转录因子，包括：

- **FOXC2**
- **PBX1 / PBX2**
- **MEIS1 / MEIS2**
- **HOXD13, HOXB13**
- **HNF1A/B**
- **ONECUT3 等**

我们将这些 TF 的平均表达量与各亚群的细胞类型预测进行关联分析，用以推测：

- 哪些亚群具备潜在的 DOCK9 调控能力  
- 哪些 TF 与特定细胞类型在表达模式上呈现共变关系  

### 📌 TF × Celltype 平均表达矩阵（候选 DOCK9 调控因子）

![mean](./images/mean.png)

---
