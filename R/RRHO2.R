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
          log10.ind = FALSE, maximum=200, boundary = 0.1, res=30, method="fisher", alternative)
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
  .hypermat<- numericListOverlap(list1[, 1], list2[, 1], stepsize, method=method, alternative = alternative)
  hypermat<- .hypermat$log.pval
  .hypermat_flipX <- numericListOverlap(rev(list1[, 1]), list2[, 1], stepsize, method=method, alternative = alternative)
	
	####Add options for old method#####
	if(alternative == "two.sided" | alternative == "enrichment"){
		
if(log10.ind) hypermat<- hypermat *log10(exp(1))  
    
  if(BY){
    hypermatvec  <- matrix(hypermat,
                           nrow=nrow(hypermat)*ncol(hypermat),ncol=1)
    hypermat.byvec  <- p.adjust(exp(-hypermatvec),method="BY")
    hypermat.by <- matrix(-log(hypermat.byvec),
                                             nrow=nrow(hypermat),ncol=ncol(hypermat))     
    
    if(log10.ind) hypermat.by<- hypermat.by *log10(exp(1))
    result$hypermat.by<- hypermat.by
  }
  
  if(plots) {
    try({
    hypermat.signed<- hypermat * .hypermat$signs 
    
    ## Function to plot color bar
    ## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
    color.bar <- function(lut, min, max=-min, 
                          nticks=11, 
                          ticks=seq(min, max, len=nticks), 
                          title='') {
      scale  <- (length(lut)-1)/(max-min)
      plot(c(0,10), c(min,max), type='n', bty='n', 
           xaxt='n', xlab='', yaxt='n', ylab='')
      mtext(title,2,2.3, cex=0.8)
      axis(2, round(ticks,0), las=1,cex.lab=0.8)
      for (i in 1:(length(lut)-1)) {
        y  <- (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
      }
    }
    
    .filename <-paste("RRHOMap", labels[1], "_VS_", labels[2], ".jpg", sep="") 
    jpeg(filename = paste(outputdir,.filename,sep="/"), 
         width=8, height=8, 
         units="in", quality=100, res=150)
    
    jet.colors  <- colorRampPalette(
      c("#00007F", "blue", "#007FFF", "cyan", 
        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE))
    
    image(hypermat.signed, xlab='', ylab='', col=jet.colors(101), 
          axes=FALSE,breaks=c(seq(0,maximum,length.out = 101),1e10), main="Rank Rank Hypergeometric Overlap Map")
      
    mtext(labels[2],2,0.5)
    mtext(labels[1],1,0.5)
    ##mtext(paste("-log(BY P-value) =",max(hypermat.by)),3,0.5,cex=0.5)
    
    finite.ind<- is.finite(hypermat.signed)
    color.bar(jet.colors(101),
              min=0,
              max=maximum,
              nticks=6,
              title="-log(P-value)")
    dev.off()
    
    ## Make a rank scatter plot
    list2ind  <- match(list1[,1],list2[,1])
    list1ind  <- 1:nlist1
    corval  <- cor(list1ind,list2ind,method="spearman")
    .filename <-paste("RankScatter",labels[1],"_VS_",labels[2],".jpg",sep="") 
    jpeg(paste(outputdir,.filename,sep="/"), width=8, 
         height=8, units="in", quality=100, res=150)
    plot(list1ind,list2ind,xlab=paste(labels[1],"(Rank)"), 
         ylab=paste(labels[2],"(Rank)"), pch=20, 
         main=paste(
           "Rank-Rank Scatter (rho = ",signif(corval,digits=3),")"
           ,sep=""), cex=0.5)
    ## TODO: Replace linear fit with LOESS
    model  <- lm(list2ind~list1ind)
    lines(predict(model),col="red",lwd=3)
    dev.off()
    
    ## Make a Venn Diagram for the most significantly associated points
    ## Upper Right Corner (Downregulated in both)
    maxind.ur  <- which(
      max(hypermat.signed[ceiling(nrow(hypermat.signed)/2):nrow(hypermat.signed),
                   ceiling(ncol(hypermat.signed)/2):ncol(hypermat.signed)],
          na.rm=TRUE) == hypermat.signed, 
      arr.ind=TRUE)
    indlist1.ur  <- seq(1,nlist1,stepsize)[maxind.ur[1]]
    indlist2.ur  <- seq(1,nlist2,stepsize)[maxind.ur[2]]
    genelist.ur  <- intersect(
      list1[indlist1.ur:nlist1,1],
      list2[indlist2.ur:nlist2,1])
    ## Lower Right corner (Upregulated in both)
    maxind.lr  <- which(
      max(hypermat.signed[1:(ceiling(nrow(hypermat.signed)/2)-1), 
                   1:(ceiling(ncol(hypermat.signed)/2)-1)],
          na.rm=TRUE) == hypermat.signed, arr.ind=TRUE)
    indlist1.lr  <- seq(1,nlist1,stepsize)[maxind.lr[1]]
    indlist2.lr  <- seq(1,nlist2,stepsize)[maxind.lr[2]]
    genelist.lr  <- intersect(
      list1[1:indlist1.lr,1], 
      list2[1:indlist2.lr,1])
    
    ## Write out the gene lists of overlapping
    .filename <- paste(
      outputdir,"/RRHO_GO_MostDownregulated",labels[1],"_VS_",labels[2],".csv",
      sep="")
    write.table(genelist.ur,.filename,row.names=F,quote=F,col.names=F)
    .filename <- paste(
      outputdir,"/RRHO_GO_MostUpregulated",labels[1],"_VS_",labels[2],".csv",
      sep="")
    write.table(genelist.lr,.filename,row.names=F,quote=F,col.names=F)
    
    .filename <- paste(
      outputdir,"/RRHO_VennMost",labels[1],"__VS__",labels[2],".jpg", 
      sep="")
    jpeg(.filename,width=8.5,height=5,units="in",quality=100,res=150)
    vp1  <- viewport(x=0.25,y=0.5,width=0.5,height=0.9)
    vp2  <- viewport(x=0.75,y=0.5,width=0.5,height=0.9)
    
    pushViewport(vp1)
    h1  <- draw.pairwise.venn(length(indlist1.ur:nlist1),
                              length(indlist2.ur:nlist2),
                              length(genelist.ur), 
                              category=c(labels[1],labels[2]),
                              scaled=TRUE,
                              lwd=c(0,0),
                              fill=c("cornflowerblue", "darkorchid1"),
                              cex=1,
                              cat.cex=1.2,
                              cat.pos=c(0,0),
                              ext.text=FALSE,
                              ind=FALSE,
                              cat.dist=0.01)
    grid.draw(h1)
    grid.text("Down Regulated",y=1)
    upViewport()
    pushViewport(vp2)
    h2  <-  draw.pairwise.venn(length(1:indlist1.lr),
                               length(1:indlist2.lr),
                               length(genelist.lr),
                               category=c(labels[1],labels[2]),
                               scaled=TRUE,
                               lwd=c(0,0),
                               fill=c("cornflowerblue", "darkorchid1"),
                               cex=1,
                               cat.cex=1.2,
                               cat.pos=c(0,0),
                               ext.text=FALSE,
                               main="Negative",
                               ind=FALSE,
                               cat.dist=0.01)
    grid.draw(h2)
    grid.text("Up Regulated",y=1)
    dev.off()
  })
  if(length(h2)==0L) message('Unable to output JPG plots.')
  }
		
  result$hypermat <- hypermat
  result$hypermat.counts <- .hypermat$counts
  result$hypermat.signs <- .hypermat$signs
  
  return(result)
		
} else {
	#####Return to split method###
 .hypermat_flipX <- numericListOverlap(rev(list1[, 1]), list2[, 1], stepsize, method=method, alternative = alternative)
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
  maxind.dd <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.dd <- maxind.dd[maxind.dd[,1]>=lenStrip1 + (boundary1+1) & maxind.dd[,1]<=lenStrip1 +len1 & 
  							maxind.dd[,2]>=lenStrip2 + (boundary2+1) & maxind.dd[,2]<=lenStrip2 + len2,]

  indlist1.dd <- seq(1, nlist1, stepsize)[maxind.dd[1] - lenStrip1]
  indlist2.dd <- seq(1, nlist2, stepsize)[maxind.dd[2] - lenStrip2]
  genelist.dd <- intersect(list1[indlist1.dd:nlist1,
                                 1], list2[indlist2.dd:nlist2, 1])
  maxind.uu <- which(max(hypermat[1:boundary1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.uu <- maxind.uu[maxind.uu[,1]>=1 & maxind.uu[,1]<=boundary1 & maxind.uu[,2]>=1 & maxind.uu[,2]<=boundary2,]

  indlist1.uu <- seq(1, nlist1, stepsize)[maxind.uu[1]]
  indlist2.uu <- seq(1, nlist2, stepsize)[maxind.uu[2]]
  genelist.uu <- intersect(list1[1:indlist1.uu, 1],
                           list2[1:indlist2.uu, 1])
  #
  maxind.ud <- which(max(hypermat[1:boundary1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.ud <- maxind.ud[maxind.ud[,1]>=1 & maxind.ud[,1]<=boundary1 & maxind.ud[,2]>= lenStrip2 + (boundary2+1) & maxind.ud[,2]<=lenStrip2 + len2,]

  indlist1.ud <- seq(1, nlist1, stepsize)[maxind.ud[1]]
  indlist2.ud <- seq(1, nlist2, stepsize)[maxind.ud[2] - lenStrip2]
  genelist.ud <- intersect(list1[1:indlist1.ud,
                                 1], list2[indlist2.ud:nlist2, 1])
  maxind.du <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.du <- maxind.du[maxind.du[,1]>=lenStrip1 + (boundary1+1) & maxind.du[,1]<=lenStrip1 + len1 & maxind.du[,2]>=1 & maxind.du[,2]<=boundary2,]

  indlist1.du <- seq(1, nlist1, stepsize)[maxind.du[1] - lenStrip1]
  indlist2.du <- seq(1, nlist2, stepsize)[maxind.du[2]]
  genelist.du <- intersect(list1[indlist1.du:nlist1, 1],
                           list2[1:indlist2.du, 1])
   
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
	  ## with highest point
	  ################
	  
      .filename <- paste("RRHOMap_markH_combined_", labels[1], "_VS_",
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
	  x1.dd <- (maxind.dd[1] - 1)/(nrow(hypermat) - 1)
	  x2.dd <- (maxind.dd[2] - 1)/(ncol(hypermat) - 1)
	  #points(x1.dd,x2.dd,pch=18,cex=4)	  	  
	  x1.uu <- (maxind.uu[1] - 1)/(nrow(hypermat) - 1)
	  x2.uu <- (maxind.uu[2] - 1)/(ncol(hypermat) - 1)
	  #points(x1.uu,x2.uu,pch=18,cex=4)	  	  
	  x1.ud <- (maxind.ud[1] - 1)/(nrow(hypermat) - 1)
	  x2.ud <- (maxind.ud[2] - 1)/(ncol(hypermat) - 1)
	  #points(x1.ud,x2.ud,pch=18,cex=4)	  	  
	  x1.du <- (maxind.du[1] - 1)/(nrow(hypermat) - 1)
	  x2.du <- (maxind.du[2] - 1)/(ncol(hypermat) - 1)
	  #points(x1.du,x2.du,pch=18,cex=4)	  	  
	  
      mtext(labels[2], 2, 0.5)
      mtext(labels[1], 1, 0.5)
	  
      finite.ind <- is.finite(hypermat)
      color.bar(jet.colors(101), min = min(hypermat[finite.ind],
                                           na.rm = TRUE), max = max(hypermat[finite.ind],
                                                                    na.rm = TRUE), nticks = 6, title = "-log(P-value)")
      dev.off()
  	  
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
	    
	  ## maximum
      .filename <- paste("RRHOMap_fixMax_combined_", labels[1], "_VS_",
                         labels[2], ".tiff", sep = "")
      tiff(filename = paste(outputdir, .filename, sep = "/"),
           width = 8, height = 8, units = "in", 
           res = res)
      jet.colors <- colorRampPalette(c("#00007F", "blue",
                                       "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                       "red", "#7F0000"))
      layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
      image(hypermat, xlab = "", ylab = "", col = jet.colors(101),breaks=c(seq(0,maximum,length.out = 101),1e10),
            axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map")
	  segments(x0 = boundary1/len1 ,x1 = boundary1 /len1 ,y0 = -0.2,y1 = 1.2,lwd=4,col='white')
	  segments(x0 = -0.2,x1 = 1.2,y0 = boundary2/len2,y1 = boundary2/len2,lwd=4,col='white')	  
      mtext(labels[2], 2, 0.5)
      mtext(labels[1], 1, 0.5)
      finite.ind <- is.finite(hypermat)
      color.bar(jet.colors(101), min = 0, max = maximum, nticks = 6, title = "-log(P-value)")
      dev.off()
	    
  	  ## maximum
      .filename <- paste("RRHOMap_fixMax_combined_", labels[1], "_VS_",
                         labels[2], ".tiff", sep = "")
      tiff(filename = paste(outputdir, .filename, sep = "/"),
           width = 8, height = 8, units = "in", 
           res = res)
      jet.colors <- colorRampPalette(c("#00007F", "blue",
                                       "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                       "red", "#7F0000"))
      layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
      image(hypermat, xlab = "", ylab = "", col = jet.colors(101),breaks=c(seq(0,maximum,length.out = 101),1e10),
            axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map")
	  segments(x0 = boundary1/len1 ,x1 = boundary1 /len1 ,y0 = -0.2,y1 = 1.2,lwd=4,col='white')
	  segments(x0 = -0.2,x1 = 1.2,y0 = boundary2/len2,y1 = boundary2/len2,lwd=4,col='white')	  
      mtext(labels[2], 2, 0.5)
      mtext(labels[1], 1, 0.5)
      finite.ind <- is.finite(hypermat)
      color.bar(jet.colors(101), min = 0, max = maximum, nticks = 6, title = "-log(P-value)")
      dev.off()
  
	
	      .filename <- paste(outputdir, "/RRHO_down_",
                         labels[1], "_VS_down_", labels[2], ".csv", sep = "")
      write.table(genelist.dd, .filename, row.names = F,
                  quote = F, col.names = F)
      .filename <- paste(outputdir, "/RRHO_up_",
                         labels[1], "_VS_up_", labels[2], ".csv", sep = "")
      write.table(genelist.uu, .filename, row.names = F,
                  quote = F, col.names = F)
      .filename <- paste(outputdir, "/RRHO_VennCon", labels[1],
                         "_VS_", labels[2], ".tiff", sep = "")
      tiff(.filename, width = 8.5, height = 5, units = "in",
            res = res)
      vp1 <- viewport(x = 0.25, y = 0.5, width = 0.5, height = 0.9)
      vp2 <- viewport(x = 0.75, y = 0.5, width = 0.5, height = 0.9)
      pushViewport(vp1)
      h1 <- draw.pairwise.venn(length(indlist1.dd:nlist1),
                               length(indlist2.dd:nlist2), length(genelist.dd),
                               category = c(labels[1], labels[2]), scaled = TRUE,
                               lwd = c(0, 0), fill = c("cornflowerblue", "darkorchid1"),
                               cex = 1, cat.cex = 1.2, cat.pos = c(0, 0), ext.text = FALSE,
                               ind = FALSE, cat.dist = 0.01)
      grid.draw(h1)
      grid.text(paste("Down",labels[1],"Down",labels[2]), y = 1)
      upViewport()
      pushViewport(vp2)
      h2 <- draw.pairwise.venn(length(1:indlist1.uu), length(1:indlist2.uu),
                               length(genelist.uu), category = c(labels[1],
                                                                 labels[2]), scaled = TRUE, lwd = c(0, 0), fill = c("cornflowerblue",
                                                                                                                    "darkorchid1"), cex = 1, cat.cex = 1.2, cat.pos = c(0,
                                                                                                                                                                        0), ext.text = FALSE, main = "Negative", ind = FALSE,
                               cat.dist = 0.01)
      grid.draw(h2)
      grid.text(paste("Up",labels[1],"Up",labels[2]), y = 1)
      dev.off()

      .filename <- paste(outputdir, "/RRHO_down_",
                         labels[1], "_VS_up_", labels[2], ".csv", sep = "")
      write.table(genelist.du, .filename, row.names = F,
                  quote = F, col.names = F)
      .filename <- paste(outputdir, "/RRHO_up_",
                         labels[1], "_VS_down_", labels[2], ".csv", sep = "")
      write.table(genelist.ud, .filename, row.names = F,
                  quote = F, col.names = F)
	  #
      .filename <- paste(outputdir, "/RRHO_VennDis", labels[1],
                         "_VS_", labels[2], ".tiff", sep = "")
      tiff(.filename, width = 8.5, height = 5, units = "in",
            res = res)
      vp1 <- viewport(x = 0.25, y = 0.5, width = 0.5, height = 0.9)
      vp2 <- viewport(x = 0.75, y = 0.5, width = 0.5, height = 0.9)
      pushViewport(vp1)
      h1 <- draw.pairwise.venn(length(indlist1.du:nlist1),
                               length(1:indlist2.du), length(genelist.du),
                               category = c(labels[1], labels[2]), scaled = TRUE,
                               lwd = c(0, 0), fill = c("cornflowerblue", "darkorchid1"),
                               cex = 1, cat.cex = 1.2, cat.pos = c(0, 0), ext.text = FALSE,
                               ind = FALSE, cat.dist = 0.01)
      grid.draw(h1)
      grid.text(paste("Down",labels[1],"Up",labels[2]), y = 1)
      upViewport()
      pushViewport(vp2)
      h2 <- draw.pairwise.venn(length(1:indlist1.ud), length(indlist2.ud:nlist2),
                               length(genelist.ud), category = c(labels[1], labels[2]), scaled = TRUE,
							   lwd = c(0, 0), fill = c("cornflowerblue", "darkorchid1"), cex = 1, cat.cex = 1.2, cat.pos = c(0, 0), ext.text = FALSE,
							   main = "Negative", ind = FALSE,
                               cat.dist = 0.01)
      grid.draw(h2)
      grid.text(paste("Up",labels[1],"Down",labels[2]), y = 1)
      dev.off()
    
	       

    })
  }
  result$hypermat <- hypermat
  return(result)
}
}
