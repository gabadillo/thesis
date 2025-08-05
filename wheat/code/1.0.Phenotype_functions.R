#Functions for Phenotype evaluation
# AUDPC
# Angle
# GDD50
# FHB_1

#AUDPC
AUDPC = function(x,y,add0 = T){
  if (add0){
    x = c(0,x)
    y = c(0,y)
  }
  output = 0
  n = length(x)
  for (i in c(1:(n-1))){
    output = output + (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return (output)
}

#Angle
GDD_angle.mid.point = function(x,y){
  n = length(x)
  j = n
  for (i in c(n:2)){
    if (y[i-1] >= y[i]){
      j = i-1
    }
    else{
      break
    }
  }
  midx = (x[j]+x[n])/2
  scale_y = y[j] - y[1]
  scale_x = midx - x[1]
  output = atan(scale_y/scale_x)*180/pi
  return (output)
}

#GDD50
GDD_interpol_a = function(x,y,ymean = 0.5,min.slope.tol = 1e-2){
  bound.flag = TRUE
  n = length(x)
  y.left = which(y<ymean)
  y.right = which(y>=ymean)
  
  if (length(y.left) > 0 & length(y.right) > 0){
    bound.flag = FALSE
    x.left = x[rev(y.left)[1]]
    x.right = x[y.right[1]]
    y.left = y[rev(y.left)[1]]
    y.right = y[y.right[1]]
  }
  else if (length(y.right) == 0){
    x.left = x[n-1]
    x.right = x[n]
    y.left = y[n-1]
    y.right = y[n]
  }
  else if (length(y.left) == 0){
    x.left = x[1]
    x.right = x[2]
    y.left = y[1]
    y.right = y[2]
  }
  m = (y.right-y.left)/(x.right-x.left)
  output = (ymean - y.left)/m + x.left
  if (m < min.slope.tol & bound.flag){
    return(max(x))
  }
  return(output)
}

max_var_obs = function(sc){
  all.ns = paste("n",0:6,sep="")
  arg.max = "n0"
  var.max = 0
  curr.ns = all.ns[sc+1]
  for (no in curr.ns){
    new.var = var(subset(obs.df,Obs==no)$Value)
    if (new.var > var.max){
      var.max = new.var
      arg.max = no
    }
  }
  return(arg.max)
}
