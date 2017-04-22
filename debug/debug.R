if(F){
  library(devtools)
  install_github("Caleb-Huo/RRHO2")
}

library(RRHO2)
?RRHO2

setwd('~/Desktop/')
plotFolder <- 'plot'
system(paste('mkdir -p', plotFolder))
list.length <- 10000
list.names <- paste('Gene',1:list.length, sep='')
set.seed(15213)
gene.list1<- data.frame(list.names, sample(list.length)*sample(c(1,-1),list.length,replace=TRUE))
gene.list2<- data.frame(list.names, sample(list.length)*sample(c(1,-1),list.length,replace=TRUE))
# Enrichment alternative
RRHO.example <-  RRHO2(gene.list1, gene.list2, 
                       labels=c('x','y'), plots=TRUE, outputdir=plotFolder, BY=TRUE, log10.ind=TRUE,res=150)


list1 <- gene.list1
list2 <- gene.list2
labels=c('x','y')
plots=TRUE
outputdir=plotFolder
BY=TRUE
log10.ind=TRUE
res=150
boundary = 0.1
stepsize = defaultStepSize(list1, list2)

