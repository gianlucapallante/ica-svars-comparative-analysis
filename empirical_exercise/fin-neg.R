## Find negative impacts

find_negative <- function(index.ica,index.dc,index.pml,init_exp){
  neg.ica       <- vector("list",length(index.ica))
  for (i in 1:length(index.ica)){
    neg.ica[[i]]       <- which(as.vector(init_exp$A.ID.ica[[index.ica[i]]]) <0)
  }
  neg.dc        <- vector("list",length(index.dc))
  for (i in 1:length(index.dc)){
    neg.dc[[i]]       <- which(as.vector(init_exp$A.ID.dc[[index.dc[i]]]) <0)
  }
  neg.pml       <- vector("list",length(index.pml))
  for (i in 1:length(index.pml)){
    neg.pml[[i]]       <- which(as.vector(init_exp$A.ID.pml[[index.pml[i]]]) <0)
  }
  
  results <- list(neg.ica = neg.ica,
                  neg.dc = neg.dc,
                  neg.pml = neg.pml)
  
  return(results)
}
