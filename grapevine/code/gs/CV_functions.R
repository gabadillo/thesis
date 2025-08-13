cv2.index = function(data,
                     K = 5,
                     Reps = 10,
                     GID.name = "GID",
                     Env.name = "Env"
){
  
  n.GID = length(unique(data[,GID.name]))
  n.Env = length(unique(data[,Env.name]))
  n.GxE = nrow(data)
  
  fold.size = rep(n.GxE%/%K,K)
  fold.size = fold.size + c(rep(1,n.GxE%%K),rep(0,K-n.GxE%%K))
  
  output = list()
  for (r in 1:Reps){
    full.index = sample(seq(n.GxE),n.GxE)
    folds = list()
    sorted.index = c()
    
    end = 0
    for (k in 1:K){
      start = end+1
      end = end + fold.size[k]
      folds[[paste("K",k,sep="")]] = sort(full.index[start:end])
      sorted.index = c(sorted.index, sort(full.index[start:end]))
    }
    
    output[[paste("R",r,sep="")]] = list(full.index = sorted.index,
                                         folds = folds)
  }
  return(output)
}


cv1.index = function(data,
                     K = 5,
                     Reps = 10,
                     GID.name = "GID",
                     Env.name = "Env"
){
  
  GIDs = unlist(unique(data[,GID.name]))
  n.GID = length(GIDs)
  
  fold.size = rep(n.GID%/%K,K)
  fold.size = fold.size + c(rep(1,n.GID%%K),rep(0,K-n.GID%%K))
  
  output = list()
  for (r in 1:Reps){
    full.index = sample(seq(n.GID),n.GID)
    folds = list()
    sort.index = c()
    
    end = 0
    for (k in 1:K){
      start = end+1
      end = end + fold.size[k]
      test.lines = GIDs[full.index[start:end]]
      folds[[paste("K",k,sep="")]] = sort(which(GIDs %in% test.lines))
      sort.index = c(sort.index,sort(which(GIDs %in% test.lines)))
    }
    
    output[[paste("R",r,sep="")]] = list(full.index = sort.index,
                                         folds = folds)
  }
  return(output)
}


split.set = function(cv.object,
                     cv.scheme,
                     r,
                     k,
                     e = NULL){
  
  if (cv.scheme %in% c('cv0','cv1','cv2')){
    test.set = cv.object[[r]]$folds[[k]]
    train.set = c()
    for (nk in 1:length(cv.object[[r]]$folds)){
      if (nk != k){
        train.set = c(train.set,cv.object[[r]]$folds[[nk]])
      }
    }
    train.set = sort(train.set)
  }
  
  else if (cv.scheme %in% c('cv00')){
    test.set = cv.object[[r]]$folds[[e]][[k]]
    train.set = c()
    for (ne in 1:length(cv.object[[r]]$folds)){
      if (ne != e){
        for (nk in 1:length(cv.object[[r]]$folds[[ne]])){
          if (nk != k){
            train.set = c(train.set, cv.object[[r]]$folds[[ne]][[nk]])
          }
        }
      }
    }
    train.set = sort(train.set)
  }
  else{
    stop("Invalid CV scheme. Allowed schemes: 'cv00', 'cv0', 'cv1', 'cv2'")
  }
  
  return(list(TRS.index = train.set,
              TS.index = test.set))
}
