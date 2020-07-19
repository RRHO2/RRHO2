##' An improved version for RRHO, which aims to correct the intepretation for top left region (up in x and down in y) nad bottom right region.
##'
##' We improved the algorithm such that all four regions of RRHO plot are meaningful
##' @title RRHO2
##' @param list1 data.frame. First column is the element (possibly gene) identifier, and the second is its value on which to sort. For differential gene expression, values are often -log10(P-value) * sign(effect).
##' @param list2 data.frame. Same as list1.
##' @param stepsize Controls the resolution of the test: how many items between any two overlap tests.
##' @param labels A two element vector indicating the label of list1 and list2.
##' @param log10.ind Logical. Should pvalues be reported and plotted in -log10 scale and not -log scale?
##' @param boundary Size of the white strip. 0.1 indicates 10% of the heatmap size.
##' @param method method for odds ratio or pvalue representation "fisher" used odds ratio and "hyper" uses p-value 
##' @return list of result
##' \item{hypermat}{Matrix of -log(pvals) of the test for the first i,j elements of the lists##' the overlapping test between the first ith element of the sorted list1, 
##' using the Stratified representation.}
##' \item{genelist_uu}{Genes are are up-regulated in list 1 and up-regulated list 2, at the most significant pixel of the uu quadrant}
##' \item{genelist_dd}{Genes are are down-regulated in list 1 and down-regulated list 2, at the most significant pixel of the dd quadrant}
##' \item{genelist_du}{Genes are are down-regulated in list 1 and up-regulated list 2, at the most significant pixel of the du quadrant}
##' \item{genelist_ud}{Genes are are up-regulated in list 1 and down-regulated list 2, at the most significant pixel of the ud quadrant}
##' @author Kelly and Caleb
##' @export
##' @examples
##' 
##' ## A total of 2000 genes in both list 1 and list 2.
##' ## In list 1, genes 1-200 are up-regulated; genes 201-400 are down-regulated; the rest of the 1600 genes are noise genes. 
##' ## In list 2, genes 1-200 are up-regulated; genes 201-400 are down-regulated; the rest of the 1600 genes are noise genes.  
##' 
##' set.seed(15213)
##' nGenes <- 2000
##' nDE <- 200
##' Genes <- paste0("Genes",1:nGenes)
##' 
##' list1_pvalue_1_200 <- runif(nDE,0,0.05)
##' list1_pvalue_201_400 <- runif(nDE,0,0.05) 
##' list1_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1)
##' list1_DDE <- c(-log10(list1_pvalue_1_200), -log10(list1_pvalue_201_400) * (-1), -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))
##' 
##' gene_list1 <- data.frame(Genes=Genes,DDE = list1_DDE, stringsAsFactors = FALSE)
##' 
##' list2_pvalue_1_200 <- runif(nDE,0,0.05)
##' list2_pvalue_201_400 <- runif(nDE,0,0.05) 
##' list2_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1)
##' list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))
##' 
##' gene_list2 <- data.frame(Genes=Genes,DDE = list2_DDE, stringsAsFactors = FALSE)
##' 
##' RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("list1", "list2"), log10.ind=TRUE)
##' 


RRHO2_initialize <- function (list1, list2, stepsize = defaultStepSize(list1, list2),
                              labels = NULL, log10.ind = FALSE,
                              boundary = 0.1, method="hyper")
{
  if (any(duplicated(list1[,1])))
    stop("Non-unique gene identifier found in list1")
  if (any(duplicated(list1[,2])))
    stop("Non-unique gene identifier found in list2")
  if(!is.null(labels)){
    stopifnot(length(labels) == 2)
  }
  
  list1 <- list1[order(list1[, 2], decreasing = TRUE), ]
  list2 <- list2[order(list2[, 2], decreasing = TRUE), ]
  
  nlist1 <- length(list1[, 1])
  nlist2 <- length(list2[, 1])
  
  N <- max(nlist1, nlist2)
  ####Add options for old method#####
  #####Return to stratified method###
  .hypermat_normal<- numericListOverlap(list1[, 1], list2[, 1], stepsize, method=method)
  hypermat_normal<- .hypermat_normal$log.pval
  
  .hypermat_flipX <- numericListOverlap(rev(list1[, 1]), list2[, 1], stepsize, method=method)
  hypermat_flipX <- .hypermat_flipX$log.pval
  hypermat_flipX2 <- hypermat_flipX[nrow(hypermat_flipX):1,]
  
  stepList1 <- seq(1, nlist1, stepsize)
  stepList2 <- seq(1, nlist2, stepsize)
  
  len1 <- length(stepList1)
  len2 <- length(stepList2)
  
  lenStrip1 <- round(len1*boundary)
  lenStrip2 <- round(len2*boundary)
  
  boundary1 <- sum(list1[stepList1,2] > 0)
  boundary2 <- sum(list2[stepList2,2] > 0)
  
  hypermat <- matrix(NA,nrow=nrow(hypermat_normal) + lenStrip1,ncol=ncol(hypermat_normal) + lenStrip2)
  hypermat[1:boundary1,1:boundary2] <- hypermat_normal[1:boundary1,1:boundary2] ## u1u2, quadrant III
  hypermat[lenStrip1 + (boundary1+1):len1,lenStrip2 + (boundary2+1):len2] <- hypermat_normal[(boundary1+1):len1,(boundary2+1):len2] ## d1d2, quadrant I
  hypermat[1:boundary1,lenStrip2 + (boundary2+1):len2] <- hypermat_flipX[len1:(len1 - boundary1 + 1),(boundary2+1):len2] ## u1d2, quadrant II
  hypermat[lenStrip1 + (boundary1+1):len1,1:boundary2] <- hypermat_flipX[(len1 - boundary1):1,1:boundary2] ## u1d2, quadrant IV
  
  if (log10.ind){
    hypermat <- hypermat * log10(exp(1))
  }
  
  #### dd: down in 1 and down in 2
  maxind.dd <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  maxind.dd <- maxind.dd[maxind.dd[,1]>=lenStrip1 + (boundary1+1) & maxind.dd[,1]<=lenStrip1 +len1 & 
                           maxind.dd[,2]>=lenStrip2 + (boundary2+1) & maxind.dd[,2]<=lenStrip2 + len2,]
  
  indlist1.dd <- seq(1, nlist1, stepsize)[maxind.dd[1] - lenStrip1]
  indlist2.dd <- seq(1, nlist2, stepsize)[maxind.dd[2] - lenStrip2]
  gene_list1_dd <- list1[1:indlist1.dd, 1]
  gene_list2_dd <- list2[1:indlist2.dd, 1]
  gene_list_overlap_dd <- intersect(gene_list1_dd,
                                    gene_list2_dd)
  genelist_dd <- list(gene_list1_dd=gene_list1_dd, 
                      gene_list2_dd=gene_list2_dd,
                      gene_list_overlap_dd=gene_list_overlap_dd
                      )
  #### uu: up in 1 and up in 2
  maxind.uu <- which(max(hypermat[1:boundary1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  maxind.uu <- maxind.uu[maxind.uu[,1]>=1 & maxind.uu[,1]<=boundary1 & maxind.uu[,2]>=1 & maxind.uu[,2]<=boundary2,]
  
  indlist1.uu <- seq(1, nlist1, stepsize)[maxind.uu[1]]
  indlist2.uu <- seq(1, nlist2, stepsize)[maxind.uu[2]]
  gene_list1_uu <- list1[1:indlist1.uu, 1]
  gene_list2_uu <- list2[1:indlist2.uu, 1]
  gene_list_overlap_uu <- intersect(gene_list1_uu,
                                    gene_list2_uu)
  genelist_uu <- list(gene_list1_uu=gene_list1_uu, 
                      gene_list2_uu=gene_list2_uu,
                      gene_list_overlap_uu=gene_list_overlap_uu
  )
  #### ud: up in 1 and down in 2
  maxind.ud <- which(max(hypermat[1:boundary1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.ud <- maxind.ud[maxind.ud[,1]>=1 & maxind.ud[,1]<=boundary1 & maxind.ud[,2]>= lenStrip2 + (boundary2+1) & maxind.ud[,2]<=lenStrip2 + len2,]
  
  indlist1.ud <- seq(1, nlist1, stepsize)[maxind.ud[1]]
  indlist2.ud <- seq(1, nlist2, stepsize)[maxind.ud[2] - lenStrip2]
  gene_list1_ud <- list1[1:indlist1.ud, 1]
  gene_list2_ud <- list2[1:indlist2.ud, 1]
  gene_list_overlap_ud <- intersect(gene_list1_ud,
                                    gene_list2_ud)
  genelist_ud <- list(gene_list1_ud=gene_list1_ud, 
                      gene_list2_ud=gene_list2_ud,
                      gene_list_overlap_ud=gene_list_overlap_ud
  )
  
  #### du: down in 1 and up in 2
  maxind.du <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.du <- maxind.du[maxind.du[,1]>=lenStrip1 + (boundary1+1) & maxind.du[,1]<=lenStrip1 + len1 & maxind.du[,2]>=1 & maxind.du[,2]<=boundary2,]
  indlist1.du <- seq(1, nlist1, stepsize)[maxind.du[1] - lenStrip1]
  indlist2.du <- seq(1, nlist2, stepsize)[maxind.du[2]]
  gene_list1_du <- list1[1:indlist1.du, 1]
  gene_list2_du <- list2[1:indlist2.du, 1]
  gene_list_overlap_du <- intersect(gene_list1_du,
                                    gene_list2_du)
  genelist_du <- list(gene_list1_du=gene_list1_du, 
                      gene_list2_du=gene_list2_du,
                      gene_list_overlap_du=gene_list_overlap_du
  )
  
  
  result <- list(hypermat=hypermat, 
                 method=method,
                 labels=labels,
                 log10.ind=log10.ind,
                 genelist_uu=genelist_uu,
                 genelist_dd=genelist_dd,
                 genelist_ud=genelist_ud, 
                 genelist_du=genelist_du
                 )
  return(result)
}
