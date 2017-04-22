## Author: Jonathan Rosenblatt and Jason Stein 

## Suggest default step size
defaultStepSize <-function(list1, list2){
  n1<- dim(list1)[1]
  n2<- dim(list2)[1]
  result <- ceiling(min(sqrt(c(n1,n2))))	
  return(result)
}	
### Testing:
## list.length <- 100
## list.names <- paste('Gene',1:list.length, sep='')
## gene.list1<- data.frame(list.names, sample(100))
## gene.list2<- data.frame(list.names, sample(100))
## defaultStepSize(gene.list1, gene.list2)


