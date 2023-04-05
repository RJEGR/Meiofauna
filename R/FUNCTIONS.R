


m[,pos] <- paste0(into[pos], "__")

fill_ranks <- function(x, y) {
  
  # x: vector of ranks
  # y: backbone of ranks
  
  y <- c("p", "c", "o", "f", "g", "s")
  
  # m <- matrix(ncol = length(y))
  
  m <- NULL
  
  for(i in length(y):1) { 
    
    searchrank <- grepl(paste0("^", y[i]), x)
    
    if(sum(searchrank) > 0) {
      
      pos <- which(searchrank)
      
      rank <- x[pos]
      
      m[i] <- rank
      
    } else {
      
      
      # m[,i] <- paste0(y[i], "__")
      
      rank <- paste0(y[i], "__")
      
      m[i] <- rank
      
    }
    
  }
  
  m <- paste0(m, collapse = ";")
  
  return(m)
  
}


# x <-  c("p__Ascomycota","c__Eurotiomycetes","o__Chaetothyriales","f__Herpotrichiellaceae","g__Capronia", "s__mysp")
# x <- c("p__Ascomycota","c__Eurotiomycetes","f__Herpotrichiellaceae","g__Capronia")
x <- c("c__Eurotiomycetes")

fill_ranks(x)


paste0(vp, collapse = ";")
