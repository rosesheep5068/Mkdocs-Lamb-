---
hide:
  - footer
  - feedback
comments: false
---

<!-- ===== 动画标题 ===== -->

<center>
  <font class="custom-font ml3">Yang Li Research Hub</font>
</center>

<script src="https://cdn.statically.io/libs/animejs/2.0.2/anime.min.js"></script>

<style>
.custom-font {
  font-size: 40px;
  font-weight: 600;
  color: #444;
  letter-spacing: 2px;
}
@media (max-width: 768px) {
  .custom-font { font-size: 28px; }
}
</style>

---

## 🔍 Research Overview

- **Facial Genomics** — Enhancer–TF interactions (DOCK9 / FOXC2) shaping human craniofacial features.  
- **Single-Cell Multi-Omics** — Developmental trajectories across tissues and species (Human & Cattle).  
- **Medical Imaging & AI** — CT-based radiomics and deep learning for renal prognosis (ΔeGFR).

---

<div class="grid cards" markdown>

- :material-flask-outline:{ .lg .middle } __DOCK9–FOXC2 Regulatory Mechanism__
  ---
  
  SNP–enhancer functional validation and transcriptional regulation.  
  [:octicons-arrow-right-24: View Project](projects/eye_snp.md)

- :material-brain:{ .lg .middle } __MedIR Kidney Prognosis (ΔeGFR)__
  ---
  
  Multi-modal deep learning framework integrating CT and clinical data.  
  [:octicons-arrow-right-24: Explore](projects/ctrc_kidney_dl.md)

- :material-cow:{ .lg .middle } __Cattle Cross-Tissue Single-Cell Atlas__
  ---
  
  Transcriptomic integration and cellular landscape across tissues.  
  [:octicons-arrow-right-24: Open Notes](notes/cattle_scrna/index.md)

</div>

---

## 🧰 Data & Code Resources

- [scRNA-seq Integration (Scanpy + Harmony)](notes/ctrt_scrna/integration.md)  
- [ATAC-seq Pseudotime Analysis (SnapATAC2)](notes/cattle_scrna/integration.md)  
- [CT Radiomics Feature Selection & LASSO](notes/medir_ct/model.md)

---

## 👤 About

**Li Yang (Lambé)**  
College of Life Sciences, Jiangsu Normal University  
Research collaborator at Prof. Li Xichuan’s lab, Tianjin Normal University  

📧 [506837558@qq.com](mailto:506837558@qq.com)  
🌐 [GitHub](https://github.com/rosesheep5068)

---

<!-- ===== 背景网格效果 ===== -->

<style>
body {
  position: relative;
}
body::before {
  --size: 35px;
  --line: color-mix(in hsl, canvasText, transparent 85%);
  content: '';
  height: 100vh;
  width: 100%;
  position: absolute;
  background: 
    linear-gradient(90deg, var(--line) 1px, transparent 1px var(--size)) 50% 50% / var(--size) var(--size),
    linear-gradient(var(--line) 1px, transparent 1px var(--size)) 50% 50% / var(--size) var(--size);
  -webkit-mask: linear-gradient(-20deg, transparent 50%, white);
  mask: linear-gradient(-20deg, transparent 50%, white);
  top: 0;
  pointer-events: none;
  z-index: -1;
}
@media (max-width: 768px) {
  body::before { display: none; }
}
.md-grid { max-width: 1000px; }
</style>