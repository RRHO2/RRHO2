
## Compute the overlaps between two *numeric* lists:
numericListOverlap<- function(sample1, sample2, stepsize, method="hyper"){
  n<- length(sample1)
  
  overlap_hyper <- function(a,b) {
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))    
    log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)         
    signs<- 1L
    
    return(c(counts=count, 
             log.pval=as.numeric(log.pval),
             signs=as.integer(signs)
    ))    
  }
  
  overlap_fisher <- function(a,b) {
    s1 <- sample1[1:a]
    s2 <- sample2[1:b]
    commonSample <- intersect(s1, s2)
    lenA <- length(commonSample)
    lenB <- length(s1)
    lenC <- length(s2)
    lenD <- length(sample1)
    fisherTable <- matrix(c(lenA, lenB-lenA, lenC-lenA, lenD*2 - lenB - lenC + lenA),2,2)
    fisherResult <- fisher.test(fisherTable,alternative = "greater")
    logOdds <- log(fisherResult$estimate)
    logOdds[!is.finite(logOdds)] <- sign(logOdds[!is.finite(logOdds)]) * 100
    signs<- 1L
    
    return(c(counts=lenA, 
             log.pval=as.numeric(logOdds),
             signs=as.integer(signs)
    ))    
  }
  
  indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize))
  if(method=="hyper"){
    overlaps<- apply(indexes, 1, function(x) overlap_hyper(x['i'], x['j']))
    
  } else if(method=="fisher"){
    overlaps<- apply(indexes, 1, function(x) overlap_fisher(x['i'], x['j']))
  }
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows)  
  matrix.signs<- matrix(overlaps['signs',], ncol=nrows)  
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals))
}
### Testing:
# n<- 112
# sample1<- sample(n)
# sample2<- sample(n)  
# .test<- RRHO:::numericListOverlap(sample1, sample2, 10, method="hyper")
# .test<- RRHO:::numericListOverlap(sample1, sample2, 10, method="fisher")
# dim(.test$log.oval)
# library(lattice)
# levelplot(.test$counts)
# levelplot(.test$log.pval)
# table(is.na(.test$log.pval))
