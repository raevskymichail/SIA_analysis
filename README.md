<h1 align="center">
  Analysis of spinal intramedullary astrocytoma (SIA)
  <br>
</h1>

[![](https://img.shields.io/github/languages/code-size/raevskymichail/SIA_analysis)](https://img.shields.io/github/languages/code-size/raevskymichail/SIA_analysis)
[![](https://img.shields.io/github/languages/top/raevskymichail/SIA_analysis)](https://img.shields.io/github/languages/top/raevskymichail/SIA_analysis)
[![](https://img.shields.io/github/issues/raevskymichail/SIA_analysis)](https://img.shields.io/github/issues/raevskymichail/SIA_analysis)
[![](https://img.shields.io/github/license/raevskymichail/SIA_analysis)](https://img.shields.io/github/license/raevskymichail/SIA_analysis)

This repository contains primary source code for *"Transcriptomic portraits and molecular pathway activation features of spinal intramedullary astrocytomas"* manuscript. 

## Introduction

Spinal intramedullary astrocytoma (SIA) is a rare subtype of glioma comprising about **2–4%** of all primary central nervous system (CNS) neoplasms and approximately 6-8% of gliomas occurring in the spinal cord. However, clinical data on prognostic biomarkers as well as tumor molecular data associated with treatment outcomes are needed for patients with spinal astrocytoma due to a particularly low frequency of these tumors and lack of successful therapeutic regimens. In addition, diagnosis and treatment of these neoplasms is often challenging given their ambiguous manifestations such as back pain, limb weakness, paresthesia, and bowel and bladder dysfunction

In this study we report **31 spinal intramedullary astrocytoma (SIA) RNA sequencing profiles** for 25 adult patients with documented clinical annotations. To our knowledge, this is the first clinically annotated RNA sequencing dataset of spinal astrocytomas derived from the intradural intramedullary compartment. We compared these tumor profiles with the previous healthy CNS RNAseq data for spinal cord and brain and identified SIA-specific gene sets and molecular pathways.

Our findings suggest a trend for SIA-upregulated pathways governing interactions with the immune cells, and downregulated pathways for the neuronal functioning in the context of normal CNS activity. In two patient tumor biosamples, we identified diagnostic *KIAA1549-BRAF* fusion oncogenes, and we also found sixteen new SIA-associated fusion transcripts. In addition, we bioinformatically simulated activities of targeted cancer drugs in SIA samples and predicted that several tyrosine kinase inhibitory drugs and thalidomide analogs could be potentially effective as the second-line treatment agents to aid in the prevention of SIA recurrence and progression.

## 📝 Requirements

Main dependencies are:
* *DESeq2* >= 1.34.0
* *EnhancedVolcano* >= 1.12.0
* *pca2d* >= 3.6.0
* *preprocessCore* >= 1.50.0
* *biomaRt* >= 2.44.4
* *pacman* >= 0.5.1

Other (minor) dependecies will be automatically installed if they are missing by `pacman` package within an execution of source scripts.

## 🚀 Quick start

Source scripts that can be used for a reproduction is located at `./scripts`.
Annotation as well as supplementary data files used across analyses can be found at `./annotations`.
Callculations for drugs scores and molecular pathway activation levels were done using Oncobox platform (open.oncobox.com/)

## 🆘 Help

Please feel free to contact Mikhail Raevskiy (raevskii.mm@phystech.edu) if you have any questions about the software.

## 📃 License

This project is [Apache 2.0](https://github.com/raevskymichail/SIA_analysis/blob/main/LICENSE) licensed.