# Proinflammatory and cytotoxic CD38+HLA-DR+ effector memory CD8+ T cells are peripherally expanded in human cardiac allograft vasculopathy 

Cardiac allograft vasculopathy (CAV) is the leading cause of late graft failure post-heart transplantation (HTx). While interferon-gamma (IFNG) is considered to play a central role in chronic allogeneic immune responses mediating CAV, circulating lymphocytes participating in the IFNG-axis remain largely unknown in human CAV. In this study, we sought to characterize circulating lymphocytes in ISHLT grade 2 or 3 CAV to identify potentially pathogenic populations participating in the IFNG-axis. 

PBMCs collected from high grade CAV (n=6) and normal HTx (n=12) patients were analyzed using CITE-seq and VDJ-seq. CITE-seq and VDJ-seq were performed using the BD Rhapsody Single-Cell Analysis System. Cell surface epitopes were labeled with 56 antibody-derived tags (ADT). Targeted mRNA amplification was performed for the BD Human Immune Response Panel (398 genes) and a custom panel of ~100 genes. One CAV sample and two normal HTx samples were pooled by hash-tagging and processed together (6 batches total). The sequencing data also contains a non-specific graft dysfunction (NGD) group, which was removed from the analysis before data integration. 

1. CITE-seq data was analyzed using the Seurat (v5.2.0) R package.
2. Batch effects were corrected based on the canonical correlation analysis.
3. Major cell types were identified using selected ADT markers
4. Unsupervised clustering (WNN) was performed on CD8T, CD4T, B, and NK cells to identify subpopulations.
5. Key subpopulations were identified by cell compositional and differential gene expression analysis.
6. Clonal expansion was compared between clusters and also between the disease groups. 

In this study, we found CD8+ T cells expressed IFNG most highly. Among the CD8+ T cell clusters, the CD38+HLA-DR+ CD8+ effector memory T (Tem) cell cluster was significantly increased in CAV compared to normal HTx. This cluster showed clonal expansion, activated IFNG signaling as well as enhanced cytotoxicity with granzyme B (GZMB) and perforin 1 (PRF1) overexpression. Circulating CD38+HLA-DR+ CD8+ Tem cells may contribute to the pathogenic IFNG-axis of human CAV.

Data are available at NIH GEO repositories:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291290

Citation: 
