##' An improved version for RRHO, which aims to correct the intepretation for top left region (up in x and down in y) nad bottom right region.
##'
##' We improved the algorithm such that all four regions of RRHO plot are meaningful
##' @title RRHO2
##' @param RRHO_obj RRHO object. See RRHO2_initialize for details.
##' @param maximum Maximum value of the heatmap.
##' @param minimum Maximum value of the heatmap.
##' @param colorGradient A vector of gradient colors. Default NULL is the rainbow color.
##' @param ... other parameter for the figure control.
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
##' RRHO2_heatmap(RRHO_obj)
##' 

RRHO2_heatmap <- function(RRHO_obj, maximum=NULL,minimum=NULL, colorGradient = NULL, ...)
{
  
  hypermat <- RRHO_obj$hypermat
  labels <- RRHO_obj$labels
  method <- RRHO_obj$method
  
  if(!is.null(maximum)){
    hypermat[hypermat>maximum] <- maximum
  } else {
    maximum <- max(hypermat,na.rm=TRUE)
  }
    
  if(!is.null(minimum)){
    hypermat[hypermat<minimum] <- minimum
  } else {
    minimum <- min(hypermat,na.rm=TRUE)
  }
  
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
  
	if(is.null(colorGradient)){
	  jet.colors  <- colorRampPalette(
	    c("#00007F", "blue", "#007FFF", "cyan", 
	      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		colorGradient <- jet.colors(101)
	}
  layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  
  breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
  image(hypermat, col = colorGradient,breaks=breaks,
        axes = FALSE, ...)
  
  if(!is.null(labels)){
    mtext(labels[2],2,0.5)
    mtext(labels[1],1,0.5)
  }
  
  if(method == "hyper"){
    atitle <- ifelse(RRHO_obj$log10.ind, "-log10(P-value)", "-log(P-value)")
    color.bar(colorGradient, min = min(0,minimum), max = maximum, nticks = 6, title = atitle)
  } else if (method == "fisher"){
    atitle <- "log Odds"
    color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
  } else {
    stop("internal error (1), please report this error to https://github.com/RRHO2/RRHO2/issues")
  }
  invisible(hypermat)
}
