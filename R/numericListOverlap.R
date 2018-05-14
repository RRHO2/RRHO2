
## Compute the overlaps between two *numeric* lists:
numericListOverlap<- function(sample1, sample2, stepsize, method="hyper", alternative, tol = 0.5, maximum){
  n<- length(sample1)
  
  overlap_hyper <- function(a,b) {
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))    
    signs<- 1L
    switch(alternative,
           enrichment={
             log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)         
             signs<- 1L
           },
           two.sided={
             the.mean<- a*b/n
             signs<- sign(count - the.mean)
             if(signs < 0){
               lower<- count 
               upper<- 2*the.mean - count 
             } else{
               lower<- 2*the.mean - count 
               upper<- count 
             }
             log.pval<- -log(phyper(q=lower+tol, m=a, n=n-a+1, k=b, lower.tail=TRUE) +
                 phyper(q= upper-tol, m=a, n=n-a+1, k=b, lower.tail=FALSE))   
             #max<-log.pval[is.finite(log.pval)==TRUE]
            log.pval[!is.finite(log.pval)]<- maximum
             },
          split={
              log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)    
              log.pval[is.na(log.pval)]<-0
              signs<- 1L})

    return(c(counts=count, 
             log.pval=as.numeric(log.pval),
             signs=as.integer(signs)
    ))    
  }
  
  overlap_fisher <- function(a,b) {
    s1 <- sample1[1:a]
    s2 <- sample2[1:b]
    lenA <- as.integer(sum(as.numeric(s1 %in% s2))) 
    lenB <- length(s1)
    lenC <- length(s2)
    
    Odds <- lenA/(lenB-lenA)/(lenC-lenA)*(n - lenB - lenC + lenA)
    if(Odds == 0){
      Odds <- 1
      }
    logOdds <- log(abs(Odds))*sign(Odds)
    #logOdds[Odds == 0]<- maximum 
    #logOdds[logOdds<0]<- -maximum
    #logOdds[!is.finite(logOdds)] <- sign(logOdds[!is.finite(logOdds)]) * 100
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
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals, signs = matrix.signs))
}
### Testing:
# n<- 112
# sample1<- sample(n)
# sample2<- sample(n)  
# .test<- numericListOverlap(sample1, sample2, 10, method="hyper")
# .test<- numericListOverlap(sample1, sample2, 10, method="fisher")
# dim(.test$log.oval)
# library(lattice)
# levelplot(.test$counts)
# levelplot(.test$log.pval)
# table(is.na(.test$log.pval))
