# DM Proteomics Analysis Overview

本页面简要介绍 DM 蛋白质组学分析的生物信息学流程、数据来源以及主要可视化结果。

## 📊 数据来源

- Label-free 定量蛋白质组数据（Ctrl / UVB / DM 三组）
- 使用统一标准化处理后的蛋白表达矩阵进行后续分析

## 🧬 分析流程（Bioinformatics Workflow）

- 数据导入与缺失值过滤  
- 表达矩阵标准化（log2 转换、归一化）  
- PCA 分析评估样本间整体差异  
- 差异蛋白筛选（UVB vs Ctrl；DM vs UVB）  
- 差异可视化（火山图、热图、Venn 图）  
- KEGG 与 GO 功能富集分析  

## 📌 主要结果展示

以下为 DM 项目差异蛋白分析与功能注释的可视化示例：

### • PCA、热图与火山图

![DM Proteomics Result](./images/DM.png)
