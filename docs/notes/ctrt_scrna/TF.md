# 🧰 JASPAR TFBS Extraction Tool 使用指南

## 📌 环境说明

- 已配置好 conda 环境：
  
  ```bash
  conda activate clgeno
  ```

- 工具目录：
  
  ```
  /home/jsnu_bioinfo/JASPAR_TFBS/bin/
  ```
  
  内含 `extract_TFBSs_JASPAR.sh` 等脚本。

- 依赖（已可用）：  
  `parallel`、`bedtools`、`bigBedToBed`

---

## 📂 数据位置

- **bigBed 数据库**（已下载好）：
  
  ```
  /cache/workspace/jsnu_bioinfo/data/DEHCD/data/JASPAR2026_hg38.bb
  ```

- **结果目录**（可写入结果）：
  
  ```
  /cache/workspace/jsnu_bioinfo/proj_face/DOCK9/DEHCD/results/04_dock9_fimo/
  ```

---

## 📝 创建输入 BED 文件

工具需要一个 **Tab 分隔**的 BED 文件（至少 3 列：chr, start, end）。

例如：提取 **rs1408718 ±20bp** 区域：

```bash
cat <<EOF > /cache/workspace/jsnu_bioinfo/data/DEHCD/data/rs1408718_ucsc.bed
chr13    98973111    98973152    rs1408718
EOF
```

⚠️ 注意：

- bigBed 文件是 **hg38 UCSC 风格染色体** → 用 `chr13`

- 

---

## 🚀 运行命令

调用 `extract_TFBSs_JASPAR.sh`：

```bash
/home/jsnu_bioinfo/JASPAR_TFBS/bin/extract_TFBSs_JASPAR.sh \
  -i /cache/workspace/jsnu_bioinfo/data/DEHCD/data/rs1408718_ucsc.bed \
  -b /cache/workspace/jsnu_bioinfo/data/DEHCD/data/JASPAR2026_hg38.bb \
  -o /cache/workspace/jsnu_bioinfo/proj_face/DOCK9/DEHCD/results/04_dock9_fimo/rs1408718_TFBS.tsv \
  -p 16
```

### 参数说明

| 参数   | 含义                     |
| ---- | ---------------------- |
| `-i` | 输入 BED 文件              |
| `-b` | JASPAR bigBed 数据库文件    |
| `-o` | 输出文件路径                 |
| `-p` | 并行核数（建议 8–32，根据服务器情况）  |
| `-t` | （可选）指定 TF 列表，只提取目标 TF  |
| `-m` | （可选）指定 JASPAR 矩阵 ID 列表 |
| `-s` | （可选）分数阈值，过滤低分 hits     |

---

## 📊 查看结果

输出文件：

```
/cache/workspace/jsnu_bioinfo/proj_face/DOCK9/DEHCD/results/04_dock9_fimo/rs1408718_TFBS.tsv
```

查看前几行：

```bash
head -n 20 rs1408718_TFBS.tsv
```

统计总数：

```bash
wc -l rs1408718_TFBS.tsv
```

转为 CSV：

```bash
awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6}' rs1408718_TFBS.tsv \
> rs1408718_TFBS.csv
```

---

## 🎯 用途

- 快速获得 **rs1408718 区域**的 TFBS 注释

- 验证 **FOXC2 及相关 FOX 家族**是否命中该区间

- 输出可直接用于 **绘图 / downstream 分析**

---
