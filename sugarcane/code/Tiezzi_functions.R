Tiezzi = function(rs,ns){
  top = 0
  bot = 0
  for (i in 1:length(rs)){
    Vi = (1-rs[i])/(ns[i]-2)
    top = top + rs[i]/Vi
    bot = bot + 1/Vi
  }
  return(top/bot)
}


RMSE = function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  return(sqrt(mean((x - y)^2)))
}

