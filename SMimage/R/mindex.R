mindex = function(k)
{
  if(k%%2==0)
  {
    k = k + 1
  }
  ixmatrix = matrix(0, 1, 2)
  ixmatrix[1, ] = c((k+1)/2, (k+1)/2)
  nc = (k+1)/2
  if(nc>=2)
  {
    for(i in 2:nc)
    {
      data = matrix(0, 8*(i-1), 2)
      m = dim(data)[1]/4
      for(j in 1:m)
      {
        data[j, ] = c((k+1)/2-(i-2)+(j-1), (k+1)/2+(i-1))
        data[m+j, ] = c((k+1)/2+(i-1), (k+1)/2+(i-2)-(j-1))
        data[2*m+j, ] = c((k+1)/2+(i-2)-(j-1), (k+1)/2-(i-1))
        data[3*m+j, ] = c((k+1)/2-(i-1), (k+1)/2-(i-2)+(j-1))
      }
      ixmatrix = rbind(ixmatrix, data)
    }
  } 
  return(ixmatrix)
}
















