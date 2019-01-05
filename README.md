# iTALK
`iTALK` is an R toolkit for characterizing and illustrating intercellular communication, developed and maintained by [Linghua Wang lab](https://www.mdanderson.org/research/departments-labs-institutes/labs/linghua-wang-laboratory.html) at the University of Texas MD Anderson Cancer Center. `iTALK` can be used to visualize the complexity, diversity and dynamics of cell-cell communication in a wide range of biological processes. For more information, please refer to [our manuscript](https://www.biorxiv.org/content/early/2019/01/04/507871).

# Installation
To install the developmental version from GitHub:

```R
if(!require(devtools)) install.packages("devtools");
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
```
To load the installed `iTALK` in R:
```R
library(iTALK)
```
# Citation
This package is intended for research use only. For any bugs, enhancement requests and other issues, please use the [`iTALK` GitHub issues tracker](https://github.com/Coolgenome/iTALk/issues) or email [Yuanxin Wang](mailto:ywang65@mdanderson.org). If you find iTALK useful and use iTALK in your publication, please cite the paper: [iTALK: an R Package to Characterize and Illustrate Intercellular Communication](https://www.biorxiv.org/content/early/2019/01/04/507871)
