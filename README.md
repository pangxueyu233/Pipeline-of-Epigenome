



# Pipeline of Epigenome data analysis

## Overview

This workshop recorded the whole processing steps of Epigenome data analysis, including the CUT&Tag-seq (Chip-seq) and ATAC-seq data analysis, in CC-LY Lab written by Shawn (Xiangyu) Pan and Xuelan Chen. This page would be helpful and easy to be read and operated, especially for the bioinformatic new-hand. We will try to keep updating of ` Pipeline-of-Epigenome`. And this pipeline is flexible, you could broaden more analysis steps and tools integrated into this page, such as TF enrichment, bulk ATAC-seq data deconvolution and et al. We also expected you could add comments and provide request to improve this page. Hope you could had a good grip of the basic Epigenome data analysis rapidly and smoothly

## **The analysis pipeline included**

## 	1. CUT&Tag-seq (Chip-seq)

1. [x] The introduction of the tools on `Linux` system
2. [x] Quality Control of raw data
3. [x] Alignment
4. [x] PCR duplicates removing
5. [x] Peak calling
   1. [x] MACS2/3
   2. [x] SEACR
6. [x] The `.bw` and `.bedgraph` files generation
7. [x] Visualization of global distribution
8. [ ] Differential Peaks analysis
   1. [ ] Quantify peaks in all samples
   2. [ ] PCA analysis of all samples.
   3. [ ] hierarchical clustering of all samples
   4. [ ] data normalization
   5. [ ] Differential Peaks identification
   6. [ ] Volcano map of Differential Peaks 
   7. [ ] gene annotation
9. [ ] Motif/TF identification
10. [ ] et al.

## 	2. ATAC-seq

1. [ ] The introduction of the tools on `Linux` system
2. [ ] Quality Control of raw data
3. [ ] Alignment
4. [ ] PCR duplicates removing
5. [ ] Peak calling
   1. [ ] MACS2/3
6. [ ] The `.bw` and `.bedgraph` files generation
7. [ ] Visualization of global distribution
8. [ ] Differential Peaks analysis
9. [ ] Visualization of Differential Peaks
10. [ ] Motif/TF identification
11. [ ] CNV and Mutation
12. [ ] et al.

---

## For  CUT&Tag-seq data (ChIP-seq)

You could learn the standard analysis pipeline from the official document by clicking [here](https://yezhengstat.github.io/CUTTag_tutorial/). However, this pipeline is designed based on the single-index 25x25 PE Illumina sequencing data. In our lab, we prefer use the `Novaseq 6000` platform to sequence the library with PE150. Hence, we modified the analysis pipeline to be better in our data processing. 

Besides, you could get more skills and knowledge by reading this [post](https://github.com/crazyhottommy/ChIP-seq-analysis), which comprehensively collected the algorithms and tools in CUT&Tag-seq (Chip-seq) analysis. 

In this page, `fastp` was used to make a quality control of raw fastq files.  The`bowtie2` was used to align the raw data with references. After the PCR duplicates removing and the peaks calling in each sample, the `chromVAR` and `GenomicRanges` were used to quantify the counts of CUT&Tag-seq data. In latest version, `DESeq2` normalized data , which was much better to reduce the effect of peak size and library size, were used to identify the differential peaks/genes. 

## 1.  The introduction of the tools on `Linux` system

Before the learning of this page, you should install the following tools firstly on your `Liunx` system:

~~~shell
# The quality control processing of raw data with multiple threads. 用于多线程数据质控
# https://github.com/OpenGene/fastp
> fastp

# Alignment tools for CUT&Tag-seq analysis. 用于数据比对
# https://github.com/BenLangmead/bowtie2
> bowtie2

# For processing the .bam files, including indexing, alignment summary and extraction. 用于bam文件处理：文件索引,比对统计和过滤操作
# https://github.com/samtools/samtools
> samtools

# For PCR duplicates removing.用于PCR去重
# https://github.com/broadinstitute/gatk/releases
> gatk

# Conversion form .bam files to .bedgraph. 用于转换bam文件为bedgraph
# https://github.com/arq5x/bedtools2
> bedtools

# The peaks calling in CUT&Tag-seq data. 用于检测CUT&Tag样本的peaks
# https://github.com/FredHutch/SEACR
> SEACR_1.3.sh

# The peaks calling in CUT&Tag-seq data. 用于检测CUT&Tag样本的peaks
# https://github.com/macs3-project/MACS
> macs2/macs3

# Visulization and peaks distribution summary 用于可视化bam文件和峰值分布量化，包括bw文件生成
# https://github.com/deeptools/deepTools
> deeptools
~~~

After you have installed the softwares on your `ubuntu` system, you could begin to learn the processing steps of analysis. 

*The codes of CUT&Tag-seq (ChIP-seq) data analysis were divided into two parts, including pre-processing on `Linux` system and statistic calculation on `R` environment*.

You could visit them by clicking:

## [Part1. the pre-processing on `Linux` system in each CUT&Tag-seq sample](CUTTAG_pre.md)

## [Part2.  the statistic calculation on `R` environment in each CUT&Tag-seq data (updating)]() 

## [ ]. Keep updating















