#Functions for non-linear logistic regression

sigmoid = function(x,a,b){
  return(1/(1+exp((a-x)/b)))
}

da_sigmoid = function(x,a,b){
  return(exp((a-x)/b)/(b*(exp((a-x)/b)+1)^2))
}

db_sigmoid = function(x,a,b){
  ((a-x)*exp((a+x)/b))/(b*(exp(a/b)+exp(x/b)))^2
}

primitive_sigmoid = function(x,a,b){
  output = b*(log(exp(a/b)+exp(x/b)))
  return(output)
}
