kombo = function(vector)
{
  n = length(vector)
  if(n > 1){
    N = length(2:n)
    temp = vector("list", N)
    
    twice = matrix(NA, nrow = 2, ncol = n)
    for(i in 1:n){twice[,i] = rep(vector[i], 2)}
    
    temp[[1]] = twice
    for(j in 2:n){temp[[j]] = combn(vector, j)}
  }else{
    temp = vector
  }
  return(temp)
}

#Example:
#a = c(1,2,3,6)
#ka = kombo(a)