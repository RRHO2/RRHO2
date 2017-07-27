##' An improved version for RRHO, which aims to correct the intepretation for top left region (up in x and down in y) nad bottom right region.
##'
##' We improved the algorithm such that all four regions of RRHO plot are meaningful
##' @title RRHO2
##' @param list1 data.frame. First column is the element (possibly gene) identifier, and the second is its value on which to sort. For differential gene expression, values are often -log10(P-value) * sign(effect).
##' @param list2 data.frame. Same as list1.
##' @param stepsize Controls the resolution of the test: how many items between any two overlap tests.
##' @param labels Character vector with two elements: the labels of the two lists.
##' @param plots Logical. Should output plots be returned?
##' @param outputdir Path name where plots ae returned.
##' @param BY Logical. Should Benjamini-Yekutieli FDR corrected pvalues be computed?
##' @param log10.ind Logical. Should pvalues be reported and plotted in -log10 scale and not -log scale?
##' @param maximum maximum value for a union scale, default is 200.
##' @param boundary boundary interval between different quadrant.
##' @param sort determines whether gene list should be sorted by p-values or effect size
##' @return list of result
##' \item{hypermat}{Matrix of -log(pvals) of the test for the first i,j elements of the lists.}
##' @author Caleb
##' @export
##' @examples
##' 
##' plotFolder <- 'plot'
##' system(paste('mkdir -p', plotFolder))
##' list.length <- 100
##' list.names <- paste('Gene',1:list.length, sep='')
##' set.seed(15213)
##' gene.list1<- data.frame(list.names, sample(100)*sample(c(1,-1),100,replace=TRUE))
##' gene.list2<- data.frame(list.names, sample(100)*sample(c(1,-1),100,replace=TRUE))
##' # Enrichment alternative
##' RRHO.example <-  RRHO2(gene.list1, gene.list2, 
##'                       labels=c('x','y'), plots=TRUE, outputdir=plotFolder, BY=TRUE, log10.ind=TRUE)
##'

RRHO2 <- function (list1, list2, stepsize = defaultStepSize(list1, list2),
          labels, plots = FALSE, outputdir = NULL, BY = FALSE,
          log10.ind = FALSE, maximum=200, boundary = 0.1, res=30, method="fisher")
{
  if (length(list1[, 1]) != length(unique(list1[, 1])))
    stop("Non-unique gene identifier found in list1")
  if (length(list2[, 1]) != length(unique(list2[, 1])))
    stop("Non-unique gene identifier found in list2")
  if (plots && (missing(outputdir) || missing(labels)))
    stop("When plots=TRUE, outputdir and labels are required.")
  result <- list(hypermat = NA, hypermat.counts = NA, hypermat.signs = NA,
                 hypermat.by = NA, n.items = nrow(list1), stepsize = stepsize,
                 log10.ind = log10.ind, call = match.call())
  
  list1 <- list1[order(list1[, 2], decreasing = TRUE), ]
  list2 <- list2[order(list2[, 2], decreasing = TRUE), ]

  nlist1 <- length(list1[, 1])
  nlist2 <- length(list2[, 1])

  N <- max(nlist1, nlist2)
  .hypermat_normal <- numericListOverlap(list1[, 1], list2[, 1], stepsize, method=method)
  hypermat_normal <- .hypermat_normal$log.pval
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
	Q3<-hypermat[1:boundary1,1:boundary2]
	max<-max(Q3[is.finite(Q3)])
	Q3[Q3 == Inf] <-max
        hypermat[1:boundary1,1:boundary2]<-Q3

  hypermat[lenStrip1 + (boundary1+1):len1,lenStrip2 + (boundary2+1):len2] <- hypermat_normal[(boundary1+1):len1,(boundary2+1):len2] ## d1d2, quadrant I
	Q1 <-  hypermat[lenStrip1 + (boundary1+1):len1,lenStrip2 + (boundary2+1):len2]
	max<-max(Q1[is.finite(Q1)])
	Q1[Q1==Inf]<-max
	hypermat[lenStrip1 + (boundary1+1):len1,lenStrip2 + (boundary2+1):len2] <- Q1
	
  hypermat[1:boundary1,lenStrip2 + (boundary2+1):len2] <- hypermat_flipX[len1:(len1 - boundary1 + 1),(boundary2+1):len2] ## u1d2, quadrant II
	Q2<- hypermat[1:boundary1,lenStrip2 + (boundary2+1):len2]
	max<-max(Q2[is.finite(Q2)])
	Q2[Q2 == Inf] <-max
        hypermat[1:boundary1,lenStrip2 + (boundary2+1):len2]<-Q2
	
 hypermat[lenStrip1 + (boundary1+1):len1,1:boundary2] <- hypermat_flipX[(len1 - boundary1):1,1:boundary2] ## u1d2, quadrant IV
	Q4<-hypermat[lenStrip1 + (boundary1+1):len1,1:boundary2]
	min<-min(Q4[is.finite(Q4)]
	Q4[Q4 == -Inf] <- min
	hypermat[lenStrip1 + (boundary1+1):len1,1:boundary2] <-Q4

	
	
  if (log10.ind){
  	hypermat <- hypermat * log10(exp(1))
  }
    
  if (BY) {
    hypermatvec <- matrix(hypermat, nrow = nrow(hypermat) *
                            ncol(hypermat), ncol = 1)
    hypermat.byvec <- p.adjust(exp(-hypermatvec), method = "BY")
    hypermat.by <- matrix(-log(hypermat.byvec), nrow = nrow(hypermat),
                          ncol = ncol(hypermat))
    if (log10.ind)
      hypermat.by <- hypermat.by * log10(exp(1))
    result$hypermat.by <- hypermat.by
  }
  
   
  if (plots) {
    try({
      color.bar <- function(lut, min, max = -min, nticks = 11,
                            ticks = seq(min, max, len = nticks), title = "") {
        scale <- (length(lut) - 1)/(max - min)
        plot(c(0, 10), c(min, max), type = "n", bty = "n",
             xaxt = "n", xlab = "", yaxt = "n", ylab = "")
        mtext(title, 2, 2.3, cex = 0.8)
        axis(2, round(ticks, 0), las = 1, cex.lab = 0.8)
        for (i in 1:(length(lut) - 1)) {
          y <- (i - 1)/scale + min
          rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
        }
      }
	  
	  

  	  
	  ################
	  ## without highest point
	  ################
	  
      .filename <- paste("RRHOMap_combined_", labels[1], "_VS_",
                         labels[2], ".tiff", sep = "")
      tiff(filename = paste(outputdir, .filename, sep = "/"),
           width = 8, height = 8, units = "in", 
           res = res)
      jet.colors <- colorRampPalette(c("#00007F", "blue",
                                       "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                       "red", "#7F0000"))
      layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
      image(hypermat, xlab = "", ylab = "", col = jet.colors(101),
            axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map")
	  segments(x0 = boundary1/len1 ,x1 = boundary1 /len1 ,y0 = -0.2,y1 = 1.2,lwd=4,col='white')
	  segments(x0 = -0.2,x1 = 1.2,y0 = boundary2/len2,y1 = boundary2/len2,lwd=4,col='white')	  
      mtext(labels[2], 2, 0.5)
      mtext(labels[1], 1, 0.5)
	  
      finite.ind <- is.finite(hypermat)
      color.bar(jet.colors(101), min = min(hypermat[finite.ind],
                                           na.rm = TRUE), max = max(hypermat[finite.ind],
                                                                    na.rm = TRUE), nticks = 6, title = "-log(P-value)")
      dev.off()   
    })
  }
  result$hypermat <- hypermat
  return(result)
}
