# Master Thesis Repository – Gonçalo Alves

## Overview

This repository contains the supplementary material associated with my Master's thesis on the characterization of cellular populations and intercellular communication following spinal cord injury (SCI) using single-nucleus RNA sequencing (snRNA-seq).

The project investigates how different cellular populations—including immune, vascular, and ependymal cells—respond to SCI and how their communication networks change during injury progression. By combining cell-type annotation, differential expression analysis, and ligand–receptor interaction inference using CellPhoneDB, the study identifies dynamic changes in the spinal cord microenvironment that may contribute to tissue repair or chronic pathology.

The repository includes all scripts, figures, and supplementary tables used throughout the thesis.

> **Important**
>
> This repository is **not intended to be downloaded and executed as a standalone project**. It does not contain the complete computational environment, raw datasets, or software dependencies required to reproduce the analyses from scratch.
>
> Instead, it serves as a companion repository to the written thesis, providing transparency and access to the analysis scripts, figures, and supplementary material referenced throughout the dissertation.

---

# Repository Structure

```text
.
├── Figures/
├── Scripts/
├── Supplementary_Tables/
└── README.md
```

| Folder | Description |
|---------|-------------|
| `Figures/` | All figures included in the thesis. |
| `Scripts/` | Analysis scripts organized by analysis type. |
| `Supplementary_Tables/` | Supplementary tables referenced in the thesis. |

---

# Figures

The **Figures** directory contains every figure included in the thesis.

Figures are organized according to the thesis chapters:

```text
Figures/
├── Introduction/
├── Materials_Methods/
└── Results/
```

Each figure follows the naming convention:

```text
x.1
x.2
...
x.n
```

where:

- **1.x** → Introduction
- **2.x** → Materials and Methods
- **3.x** → Results

This numbering directly matches the figure numbering used throughout the thesis.

---

# Supplementary Tables

The **Supplementary_Tables** folder contains all supplementary tables referenced in the dissertation.

Currently, it includes the two supplementary tables accompanying the thesis.

---

# Scripts

The **Scripts** directory contains all analysis scripts developed during the project.

The scripts are organized according to their purpose within the analysis workflow.

## CellphoneDB

Contains all scripts related to the CellPhoneDB analysis used to infer ligand–receptor interactions and investigate intercellular communication between cell populations.

---

## General

General-purpose scripts used throughout the project.

These include scripts for:

- Clustering
- Cell lineage annotation
- Data processing
- Reusable analysis functions
- Utility workflows shared across multiple datasets

---

## Immune

Scripts used exclusively for the annotation and analysis of the immune cell dataset.

---

## Meningeal

Scripts used exclusively for the annotation and analysis of the meningeal cell dataset.

---

## Neuron

Scripts used exclusively for the annotation and analysis of neuronal cell populations.

---

## codex

Contains the Bash (`.sh`) scripts used to organize and execute the computational pipeline.

The naming convention is as follows:

- `_immune` – Immune cell annotation pipeline.
- `_meningeal` – Meningeal cell annotation pipeline.
- `_neu` – Neuronal cell annotation pipeline.
- `_cellphonedb` – CellPhoneDB and cell–cell communication analyses.

These scripts document the order in which analyses were performed and were primarily used to manage computational workflows on a high-performance computing (HPC) cluster.

---

# Repository Purpose

This repository was created to accompany my Master's thesis by providing:

- All figures included in the dissertation.
- Supplementary tables.
- Analysis scripts used throughout the project.
- Documentation of the computational workflow.

The repository is intended as a reference resource for readers of the thesis and to promote transparency of the computational analyses performed.

It should **not** be considered a fully reproducible analysis pipeline.

---

# Thesis Summary

Spinal cord injury (SCI) causes permanent neurological impairment due to the limited regenerative capacity of the central nervous system and the complex inflammatory response that follows injury. This project uses single-nucleus RNA sequencing (snRNA-seq) data from a mouse model of SCI to characterize cellular populations across different stages of injury and investigate how these cells communicate through ligand–receptor interactions.

By integrating clustering, differential gene expression, cell-type annotation, and CellPhoneDB analyses, the study identifies dynamic changes in immune, vascular, and ependymal cell populations. The results reveal increased intercellular communication during chronic injury, with fibroblasts emerging as major signaling hubs involved in fibrosis, vascular remodeling, and immune regulation. These findings contribute to a better understanding of the SCI microenvironment and highlight potential therapeutic targets to promote tissue repair and regeneration.

---

# Citation

If you use any material from this repository, please cite the corresponding Master's thesis.

---

**Author:** Gonçalo Alves

**Master's Thesis Repository**
