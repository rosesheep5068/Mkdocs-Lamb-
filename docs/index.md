---
hide:
  - footer
  - feedback
comments: false
---

<!-- ===== Title with animation ===== -->

<center>
  <font class="custom-font ml3">Yang Li Research Hub</font>
</center>

<script src="https://cdn.statically.io/libs/animejs/2.0.2/anime.min.js"></script>

<style>
.custom-font {
  font-size: 42px;
  font-weight: 600;
  color: #333;
  letter-spacing: 2px;
}
@media (max-width: 768px) {
  .custom-font { font-size: 28px; }
}
</style>

---

## 🔬 Research Overview

- **Facial Genomics** — SNP → Enhancer → TF → Gene Expression → Phenotype.
- **Single-Cell Multi-Omics** — Human / Cattle / Mouse developmental mapping.
- **Medical Imaging & AI** — CT radiomics + deep learning for ΔeGFR prediction.

---

## 🚀 Featured Projects

<div class="grid cards" markdown>

- :material-face-recognition:{ .lg .middle } __Craniofacial Regulatory Mechanisms__
  
  SNP–enhancer perturbations, TF binding, and gene regulatory logic.  
  [:octicons-arrow-right-24: View Project](projects/eye_snp.md)

- :material-brain:{ .lg .middle } __MedIR Kidney Prognosis (ΔeGFR)__
  
  Multi-modal 3D CT + clinical fusion deep learning framework.  
  [:octicons-arrow-right-24: Explore](projects/ctrc_kidney_dl.md)

- :material-dna:{ .lg .middle } __Cross-Tissue Single-Cell Atlas__
  
  Comparative regulatory landscapes across species.  
  [:octicons-arrow-right-24: Open](notes/ctrt_scrna/index.md)

</div>

---

## 📊 Data & Code Resources

- [scRNA-seq Integration Pipeline](notes/ctrt_scrna/integration.md)  
- [ATAC-seq Pseudotime (SnapATAC2)](notes/cattle_scrna/integration.md)  
- [Radiomics Feature Selection & LASSO](notes/medir_ct/model.md)

---

## 👤 About Me

**Li Yang (Lambé)**  
Graduate Student, College of Life Sciences, JSNU  
Collaborator at Prof. Li Xichuan’s Lab, TJNU  

📧 506837558@qq.com  
🌐 [GitHub](https://github.com/rosesheep5068)

---

<!-- ===== Background grid ===== -->

<style>
body {
  position: relative;
}
body::before {
  --size: 36px;
  --line: color-mix(in hsl, canvasText, transparent 88%);
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