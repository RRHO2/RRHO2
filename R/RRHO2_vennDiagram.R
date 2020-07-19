##' An improved version for RRHO, which aims to correct the intepretation for top left region (up in x and down in y) nad bottom right region.
##'
##' We improved the algorithm such that all four regions of RRHO plot are meaningful
##' @title RRHO2
##' @param RRHO_obj RRHO object. See RRHO2_initialize for details.
##' @param type Must be one of the following. 1, dd; 2, uu; 3, ud; 4, du.
##' The Venn diagram is always counting from the highest point within its quadrant (the pixel such that the p-value is the most significant).
##' The first letter represents the direction in list 1 (d: down; u: up). For example
##' dd represent down regulation in list 1 and down regulation in list 2.
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
##' RRHO2_vennDiagram(RRHO_obj, type="dd")
##' 


RRHO2_vennDiagram <- function(RRHO_obj, type = NULL)
{
  if(is.null(RRHO_obj$labels)){
    labels <- c("list1", "list2")
  } else {
    labels <- RRHO_obj$labels
  }
  
  if(is.null(type)){
    stop("need to specify type: dd, uu, ud, du, see R help file for details.")
  } else if(!type %in% c("dd","uu","du", "ud")){
    stop("type must be one of the following: dd, uu, ud, du, see R help file for details.")
  } else {
    ageneList <- switch(type, 
                        "dd" = RRHO_obj$genelist_dd,
                        "uu" = RRHO_obj$genelist_uu,
                        "du" = RRHO_obj$genelist_du,
                        "ud" = RRHO_obj$genelist_ud
    )
    atitle <- switch(type, 
                     "dd" = paste("Down",labels[1],"Down",labels[2]),
                     "uu" = paste("Up",labels[1],"Up",labels[2]),
                     "du" = paste("Down",labels[1],"Up",labels[2]),
                     "ud" = paste("Up",labels[1],"Down",labels[2])
    )
    
    
  }
  
  
  venn.plot <- draw.pairwise.venn(length(ageneList[[1]]),
                            length(ageneList[[2]]),
                            length(ageneList[[3]]), 
                            category=c(labels[1],labels[2]),
                            scaled=FALSE,
                            lwd=c(0,0),
                            fill=c("cornflowerblue", "darkorchid1"),
                            cex=1,
                            cat.cex=1.2,
                            cat.pos=c(0,0),
                            ext.text=FALSE,
                            ind=FALSE,
                            cat.dist=0.01)
  grid.draw(venn.plot);
  grid.text(atitle, y = 0.9)
  
 }