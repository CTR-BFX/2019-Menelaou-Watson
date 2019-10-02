# Blastocyst transfer in mice alters placental transcriptome and growth

**Katerina Menelaou<sup>1,2</sup>, Malwina Prater<sup>2</sup>, Simon J. Tunster<sup>1,2,</sup> Georgina E.T. Blake<sup>1,2</sup>, Colleen Geary-Joo<sup>3</sup>, James C. Cross<sup>4</sup>, Russell S. Hamilton<sup>2,5</sup>, and Erica D. Watson<sup>1,2,*</sup>*.**


<sup>1</sup>Department of Physiology, Development, and Neuroscience, University of Cambridge, Cambridge, UK

<sup>2</sup>Centre for Trophoblast Research, University of Cambridge, Cambridge, UK

<sup>3</sup>Transgenic Services, Clara Christie Centre for Mouse Genomics, University of Calgary, Calgary, Canada
<sup>4</sup>Department of Comparative Biology and Experimental Medicine, University of Calgary, Calgary, Canada
<sup>5</sup>Department of Genetics, University of Cambridge, Cambridge, UK


*Corresponding author: E.D.W., Department of Physiology, Development, and Neuroscience, University of Cambridge, Physiological Laboratory, Downing Street, Cambridge, CB2 3EG, UK; email: edw23@cam.ac.uk 



## Abstract

Assisted reproduction technologies (ART) are becoming increasingly common. Therefore, how these procedures influence gene regulation and feto-placental development are important to explore. Here, we assess the effects of blastocyst transfer on mouse placental growth and transcriptome. C57Bl/6 blastocysts were transferred into uteri of B6D2F1 pseudopregnant females and dissected at embryonic day 10.5 for analysis. Compared to non-transferred controls, placentas from transferred conceptuses weighed less even though the embryos were larger on average. This suggested a compensatory increase in placental efficiency. RNA-sequencing of whole male placentas revealed 543 differentially expressed genes (DEGs) after blastocyst transfer: 188 and 355 genes were down-regulated and up-regulated, respectively. DEGs were independently validated in male and female placentas. Bioinformatic analyses revealed that DEGs represented expression in all major placental cell types and included genes that are critical for placenta development and/or function. Furthermore, the direction of transcriptional change in response to blastocyst transfer implied an adaptive response to improve placental function to maintain fetal growth. Some DEGs appeared in genomic clusters inferring shared common cis-regulatory mechanisms. Our analysis suggested that altered DNA methylation is unlikely primary cause transcriptional disruption in transferred placentas since significantly fewer DEGs were associated with CpG islands compared to the genome at large. Furthermore, CpG methylation was unchanged in two DEGs as determined by bisulfite pyrosequencing in female placentas of transferred conceptuses. Whether blastocyst transfer disrupts other epigenetic mechanisms (e.g., histone modifications) should be explored. We conclude that embryo transfer, a protocol required for ART, significantly impacts the placental transcriptome and function.



## Data Processing

Data can be accessed from **ArrayExpress** with accession number [[**E-MTAB-8036**]](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8036). 

Fastq files (located in subdirectories of /storage/CTR-Projects/CTR_edw23/CTR_edw23_0002/) were QC'ed (FastQC and fastq_screen), trimmed with Trim Galore! and aligned to GRCm38 mouse genome using STAR aligner. Alignments and QC were processed using custom ClusterFlow (v0.5dev) pipelines and assessed using MultiQC (0.9.dev0). Gene quantification was determined with HTSeq-Counts (v0.6.1p1). Since Sequencing was performed in duplicate to provide at least 18 million reads per sample, HTSeq counts were combined before downstream analysis (script: `combine_htseq_counts.pl`). Additional quality control was performed with rRNA and mtRNA counts script, feature counts (v 1.5.0-p2) and qualimap (v2.2). 

 Principle component analysis was performed on the rlog transformed count data for top 5000 most variable genes. Differential gene expression was performed with DESeq2 package (v1.16.1, R v3.4.0) and with the same package read counts were normalised on the estimated size factors. Heatmaps were generated with 'pheatmap' R package. Karyoplots were generated with karyoploteR  -abs log2FC > 1(v1.8.8). Pathway and GO analyses were done and visualised with R packages: 'pathview', 'gage', 'DOSE', 'clusterProfiler', 'ReactomePA' and 'GOplot'. Heatmaps were generated with “ComplexHeatmap” R package (v 1.20.0). 






## Scripts 

Scripts are provided in directory `Scripts` (R Script to reproduce paper figures and Python script to combine HTseq counts).



## Figures 


|  Figure    |  Output Filename                              |   Description     |
| :-------:  |      :----                                  |         :---    |
|     2A     |    Fig_2A_PCA_top5000MV.png           |  PCA plot         |
|     2B     |    Fig_2B_MAplot.png                  |  MA plot          |
|     2C     |    Fig_2C_ComplexHeatmap.png           |  Heatmap of changes in genes expression between Control vs Transferred |
|     5A     |    Fig_5A_l2fc1_SigGenes_karyoplot.png  |    Karyoplot showing distribution of DEGs on chromosomes        |
|     5B     |    Fig_5B_EPU6b.png                |   Example of  Enhancer-promoter unit (EPU) with DEGs.     |


![Fig.2A.](https://github.com/nmalwinka/2019-Menelaou/blob/master/Figures/Fig_2A_PCA_top5000MV.png)

Fig.2A.

![Fig.2B.](https://github.com/nmalwinka/2019-Menelaou/blob/master/Figures/Fig_2B_MAplot.png)

Fig.2B.

![Fig.2C.](https://github.com/nmalwinka/2019-Menelaou/blob/master/Figures/Fig_2C_ComplexHeatmap.png)

Fig.2C.

![Fig.5A.](https://github.com/nmalwinka/2019-Menelaou/blob/master/Figures/Fig_5A_l2fc1_SigGenes_karyoplot.png)

Fig.5A.

![Fig.5B.](https://github.com/nmalwinka/2019-Menelaou/blob/master/Figures/Fig_5B_EPU6b.png)

Fig.5B.











## Contact
Contact Malwina Prater (mn367@cam.ac.uk) or Russell Hamilton (rsh46@cam.ac.uk) for bioinformatics related queries.
