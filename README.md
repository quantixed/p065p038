# p065p038

Code and data for Downie et al. manuscript `#p065p038`

**Non-disruptive inducible labeling of ER-membrane contact sites using the Lamin B Receptor**

Laura Downie, Nuria Ferrandiz, Elizabeth Courthold, Megan Jones, & Stephen J. Royle

Preprint: *bioRxiv* 2024. DOI: [https://doi.org/10.1101/2024.05.31.596797](https://doi.org/10.1101/2024.05.31.596797)

--

# R Code and data for recreation of plots

### Main Figures

- LBR cluster analysis (plasma membrane, 3D) - `SJR233` analysis of LD317
- LBR cluster analysis (plasma membrane, 2D movie) - `SJR217` analysis of LD295
- MAPPER + LBR cluster analysis - `SJR239` analysis of LD347 (6 repeats)
- Mitochondria-ER contact analysis from SBF-SEM - `SJR242`
- Individual line profiles with image outputs - `SJR248`
- Line profile comparison - `SJR266` analysis of LD237, LD239, LD352, LD365, LD360
- Thapsigargin experiment - `SJR265` analysis of LD446
- Golgi coloc - `SJR267` analysis of LD376
- ER expression of FKBP-GFP-tagged constructs under different promoters - `ERExpr`
- Sec61 and LBR cluster analysis (3D) under different promoters - `Sec61Expr`
- LBR cluster sizes after 2h or 4h (3D) - `LBR_LongTerm`


### Notes

- The R code can be executed with data placed into the `Data` folder of the RProject.
- For project folders where data is absent from this repo, it is available seprately at zenodo ([10.5281/zenodo.14011426](https://doi.org/10.5281/zenodo.14011426)), due to GitHub space constraints.
- Fiji scripts for image analysis are included in the `Script` folder of the respective project.