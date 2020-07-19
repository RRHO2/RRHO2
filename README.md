# RRHO2
Implementation of improved RRHO2, for which all regions in RRHO plots are meaningful.


## Install This Package from github
First you need R `devtools` package installed.

* In R console

```R
library(devtools)
install_github("RRHO2/RRHO2", build_vignettes = TRUE)
```


## Citation

* Cahill, K. M., Huo, Z., Tseng, G. C., Logan, R. W., & Seney, M. L. (2018). Improved identification of concordant and discordant gene expression signatures using an updated rank-rank hypergeometric overlap approach. *Scientific reports*, 8(1), 1-11.

* Plaisier, S. B., Taschereau, R., Wong, J. A., & Graeber, T. G. (2010). Rank–rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures. *Nucleic acids research*, 38(17), e169-e169.


## Full tutorial

* http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html

## Short tutorial for circadian pattern detection

* Simulate data

```
set.seed(15213)
nGenes <- 2000
nDE <- 200
Genes <- paste0("Genes",1:nGenes)

list1_pvalue_1_200 <- runif(nDE,0,0.05)
list1_pvalue_201_400 <- runif(nDE,0,0.05) 
list1_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1)
list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), 
	-log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))

gene_list1 <- data.frame(Genes=Genes,DDE = list1_DDE, stringsAsFactors = FALSE)

list2_pvalue_1_200 <- runif(nDE,0,0.05)
list2_pvalue_201_400 <- runif(nDE,0,0.05) 
list2_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1)
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), 
	-log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))

gene_list2 <- data.frame(Genes=Genes,DDE = list2_DDE, stringsAsFactors = FALSE)
```

* Create the RRHO2 object
```
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("list1", "list2"), log10.ind=TRUE)
```

* Visualize the heatmap
```
RRHO2_heatmap(RRHO_obj)
```

* Visualize the Venn Diagram
```
RRHO2_vennDiagram(RRHO_obj, type="dd")
```
