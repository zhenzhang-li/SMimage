SMM = function(x)
{
  nf = length(x)
  k = ceiling(sqrt(nf))
  if(k%%2==0)
  {
    k = k + 1
  }
  matrixda = matrix(0, k, k)
  indexm = mindex(k)
  for(i in 1:nf)
  {
    matrixda[indexm[i, 1], indexm[i, 2]] = x[i]
  }
  return(matrixda)
}
















