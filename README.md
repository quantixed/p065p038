# p065p038

Code and data for Downie et al. manuscript `#p065p038`

**Non-disruptive inducible labeling of ER-membrane contact sites using the Lamin B Receptor**

Laura Downie, Nuria Ferrandiz, Elizabeth Courthold, Megan Jones, & Stephen J. Royle

Preprint: *bioRxiv* 2024. DOI: [https://doi.org/10.1101/2024.05.31.596797](https://doi.org/10.1101/2024.05.31.596797)

--

# R Code and data for recreation of plots

### Main Figures

- **Fig 1D,E,G** LBR cluster analysis (plasma membrane, 3D) - `SJR233` analysis of LD317
- **Fig 1I** LBR cluster analysis (plasma membrane, 2D movie) - `SJR217` analysis of LD295
- **Fig 2B,C** MAPPER + LBR cluster analysis - `SJR239` analysis of LD347 (6 repeats)
- **Fig 3C** Individual line profiles with image outputs - `SJR248`
- **Fig 3D** Line profile comparison - `SJR266` analysis of LD237, LD239, LD352, LD365, LD360
- **Fig 4C,D** Mitochondria-ER contact analysis from SBF-SEM - `SJR242`
- **Fig 5B** Line profile comparison - `SJR266` analysis of LD237, LD239, LD352, LD365, LD360
- **Fig 7B** Golgi coloc - `SJR267` analysis of LD376

### Supplementary Figures

- **Fig S2B** Flow cytometry - `NF229`
- **Fig S5B** Thapsigargin experiment - `SJR265` analysis of LD446
- **Fig S6B** Line profile comparison (cell lines) - `profile_CellType`
- **Fig S8B-D** LBR cluster sizes after 2h or 4h (3D) - `LBR_LongTerm`
- **Fig S9** Individual line profiles with image outputs - `SJR248`
- **Fig S11** Line profile comparison (constructs) - `profile_Regions`
- **Fig S12D** ER expression of FKBP-GFP-tagged constructs under different promoters - `ERExpr`
- **Fig S12E,F** Sec61 and LBR cluster analysis (3D) under different promoters - `Sec61Expr`


### Notes

- The R code can be executed with data placed into the `Data` folder of each RProject.
- For RProject folders where data is absent from the `Data` folder, it is available seprately at zenodo ([10.5281/zenodo.15263507](https://doi.org/10.5281/zenodo.15263507)), due to GitHub space constraints. Download the data, add it to the `Data` folder and the analysis can proceed.
- Fiji scripts for image analysis are included in the `Script` folder of each RProject.
- Plots used in the paper can be found in `Output/Plots` or regenerated using the scripts provided.
- The values plotted in each Figure are exported to `Output/Data` 
- For convenience, all of these files are collated in `plot_data_files`

**The folder `plot_data_files` contains the numerical data underlying all plots in the manuscript.**


### Additional data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15263507.svg)](https://doi.org/10.5281/zenodo.15263507)

### Archived version of the code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15582238.svg)](https://doi.org/10.5281/zenodo.15582238)