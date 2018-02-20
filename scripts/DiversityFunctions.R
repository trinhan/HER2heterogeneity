#functions (can write this in a script later)
## 1. rekonen index. a,b are vectors describing composition of a species of interest
renkonen = function(a, b) {
  stopifnot(length(a) == length(b))
  # Scale the vectors to be fractions [0,1]
  a = a/sum(a)
  b = b/sum(b)
  #Score.
  score = 0
  for ( i in 1:length(a) ) { score = score + min(a[i], b[i]) }
  return (score)
}

# coefficient of variation
CV <- function(x){
  (sd(x, na.rm=T)/mean(x, na.rm=T))*100
}

# KL index?
KL=function(g1, g2, p1, p2, gp1, gp2){
  
  nsamp=nrow(g1)
  kullback = matrix(NA, nrow=nsamp, ncol=3)
  colnames(kullback) = c("geno", "pheno", "g-p")
  rownames(kullback) = rownames(g1)
  
  for (i in 1:nsamp) {
    # browser()
    P = g1[i,] + 0.0001
    Q = g2[i,] + 0.0001
    ok = P > 0 & Q > 0
    kullback[i,1] = sum(P[ok] * log(P[ok]/Q[ok]))
    
    P = p1[i,] + 0.0001
    Q = p2[i,] + 0.0001
    ok = P > 0 & Q > 0
    kullback[i,2] = sum(P[ok] * log(P[ok]/Q[ok]))
    
    
    P = gp1[i,] + 0.0001
    Q = gp2[i,] + 0.0001
    ok = P > 0 & Q > 0
    kullback[i,3] = sum(P[ok] * log(P[ok]/Q[ok]))
  }
  kullback
}