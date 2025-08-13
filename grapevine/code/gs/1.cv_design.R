#design CV schemes

library(tidyverse)
source('CV_functions.R')

set.seed(0)
Pheno.clean = read_csv('../../data/clean/BLUEs.csv')
colnames(Pheno.clean)[1:2] = c('GID', 'group')
Pheno.clean$group[Pheno.clean$group == 'self'] = 'sp23'
Pheno.clean$group = factor(Pheno.clean$group)
Pheno.clean$Env = 'World'
# generate train-test splits
split_list = list()

####################################################
#leave-one-out scenarios without replicates

# scenario A -> Across groups -> 2 values
split_list[["A"]] = list(
  "JFI" = list(
    "train" = which(Pheno.clean$group %in% c("sp23","406","411","900","912")),
    "test"  = which(Pheno.clean$group %in% c("111","902","929","jardin"))
    ),
  "F23" = list(
    "train"  = which(Pheno.clean$group %in% c("111","902","929","jardin")),
    "test" = which(Pheno.clean$group %in% c("sp23","406","411","900","912"))
  )
)

# scenario B -> Across families within groups -> 9 values
split_list[["B"]] = list(
  # group F23
  "sp23" = list(
    "train" = which(Pheno.clean$group %in% c("406","411","900","912")),
    "test"  = which(Pheno.clean$group %in% c("sp23"))
  ),
  "406" = list(
    "train" = which(Pheno.clean$group %in% c("sp23","411","900","912")),
    "test"  = which(Pheno.clean$group %in% c("406"))
  ),
  "411" = list(
    "train" = which(Pheno.clean$group %in% c("sp23","406","900","912")),
    "test"  = which(Pheno.clean$group %in% c("411"))
  ),
  "900" = list(
    "train" = which(Pheno.clean$group %in% c("sp23","406","411","912")),
    "test"  = which(Pheno.clean$group %in% c("900"))
  ),
  "912" = list(
    "train" = which(Pheno.clean$group %in% c("sp23","406","411","900")),
    "test"  = which(Pheno.clean$group %in% c("912"))
  ),
  #### group FJI
  "111" = list(
    "train" = which(Pheno.clean$group %in% c("902","929","jardin")),
    "test"  = which(Pheno.clean$group %in% c("111"))
  ),
  "902" = list(
    "train" = which(Pheno.clean$group %in% c("111","929","jardin")),
    "test"  = which(Pheno.clean$group %in% c("902"))
  ),
  "929" = list(
    "train" = which(Pheno.clean$group %in% c("111","902","jardin")),
    "test"  = which(Pheno.clean$group %in% c("929"))
  ),
  "jardin" = list(
    "train" = which(Pheno.clean$group %in% c("111","902","929")),
    "test"  = which(Pheno.clean$group %in% c("jardin"))
  )
)

# scenario C -> Across families -> 9 values
split_list[["C"]] = list(
  # group F23
  "sp23" = list(
    "train" = which(!Pheno.clean$group %in% c("sp23")),
    "test"  = which(Pheno.clean$group %in% c("sp23"))
  ),
  "406" = list(
    "train" = which(!Pheno.clean$group %in% c("406")),
    "test"  = which(Pheno.clean$group %in% c("406"))
  ),
  "411" = list(
    "train" = which(!Pheno.clean$group %in% c("411")),
    "test"  = which(Pheno.clean$group %in% c("411"))
  ),
  "900" = list(
    "train" = which(!Pheno.clean$group %in% c("900")),
    "test"  = which(Pheno.clean$group %in% c("900"))
  ),
  "912" = list(
    "train" = which(!Pheno.clean$group %in% c("912")),
    "test"  = which(Pheno.clean$group %in% c("912"))
  ),
  #### group FJI
  "111" = list(
    "train" = which(!Pheno.clean$group %in% c("111")),
    "test"  = which(Pheno.clean$group %in% c("111"))
  ),
  "902" = list(
    "train" = which(!Pheno.clean$group %in% c("902")),
    "test"  = which(Pheno.clean$group %in% c("902"))
  ),
  "929" = list(
    "train" = which(!Pheno.clean$group %in% c("929")),
    "test"  = which(Pheno.clean$group %in% c("929"))
  ),
  "jardin" = list(
    "train" = which(!Pheno.clean$group %in% c("jardin")),
    "test"  = which(Pheno.clean$group %in% c("jardin"))
  )
)

########################################
# k-fold cross validation with reps
R = 30
set.seed(0)

# scenario D -> all data (1 x 30 values)
cv.object = cv2.index(Pheno.clean, R=30)
for (r in 1:R){
  split_list[["D"]][["all"]][[paste0("R",r)]] = cv.object[[r]]
}

# Scenarios E and F require subsets
Pheno.clean.subsets = list()

# scenario E -> groups (2 x 30 values)
Pheno.clean.subsets[["F23"]] = subset(Pheno.clean,group %in% c("sp23","406","411","900","912"))
cv.object = cv2.index(Pheno.clean.subsets[["F23"]], R=30)
for (r in 1:R){
  split_list[["E"]][["F23"]][[paste0("R",r)]] = cv.object[[r]]
}

Pheno.clean.subsets[["JFI"]] = subset(Pheno.clean,group %in% c("111","902","929","jardin"))
cv.object = cv2.index(Pheno.clean.subsets[["JFI"]], R=30)
for (r in 1:R){
  split_list[["E"]][["JFI"]][[paste0("R",r)]] = cv.object[[r]]
}

# scenario F -> families (9 x 30 values)
# doing with a loop
for (fam in levels(Pheno.clean$group)){
  print(fam)
  Pheno.clean.subsets[[fam]] = subset(Pheno.clean,group %in% fam)
  cv.object = cv2.index(Pheno.clean.subsets[[fam]], R=30)
  for (r in 1:R){
    split_list[["F"]][[fam]][[paste0("R",r)]] = cv.object[[r]]
  }
}

###### some tests -> it should return always TRUE
# 1. Scenario A: all individuals are tested but none overlapped
sum(sort(c(split_list$A$JFI$test,split_list$A$F23$test)) == 1:nrow(Pheno.clean)) == nrow(Pheno.clean)
sum(sort(c(split_list$A$JFI$train,split_list$A$F23$train)) == 1:nrow(Pheno.clean)) == nrow(Pheno.clean)
length(intersect(split_list$A$JFI$test,split_list$A$F23$test)) == 0
length(intersect(split_list$A$JFI$train,split_list$A$F23$train)) == 0

# 2. Scenario B: unique families in train are 4 for F23 and 3 for JFI and intersection with test is always null
length(unique(Pheno.clean$group[split_list$B$sp23$train])) == 4
length(unique(Pheno.clean$group[split_list$B$`900`$train])) == 4
length(unique(Pheno.clean$group[split_list$B$jardin$train])) == 3
length(unique(Pheno.clean$group[split_list$B$`929`$train])) == 3

length(intersect(unique(Pheno.clean$group[split_list$B$sp23$train]),Pheno.clean$group[split_list$B$sp23$test])) == 0
length(intersect(unique(Pheno.clean$group[split_list$B$jardin$train]),Pheno.clean$group[split_list$B$jardin$test])) == 0

# 3. Scenario C: same checks but now train should have 8 families in all cases
length(unique(Pheno.clean$group[split_list$C$sp23$train])) == 8
length(unique(Pheno.clean$group[split_list$C$`900`$train])) == 8
length(unique(Pheno.clean$group[split_list$C$jardin$train])) == 8
length(unique(Pheno.clean$group[split_list$C$`929`$train])) == 8

length(intersect(unique(Pheno.clean$group[split_list$C$sp23$train]),Pheno.clean$group[split_list$C$sp23$test])) == 0
length(intersect(unique(Pheno.clean$group[split_list$C$jardin$train]),Pheno.clean$group[split_list$C$jardin$test])) == 0

# 4. Scenario D: we should find almost all families (or all) present in both train and test (folds)
length(unique(Pheno.clean$group[split_list$D$all$R1$folds$K1])) == 9
length(unique(Pheno.clean$group[split_list$D$all$R5$folds$K3]))== 9
length(unique(Pheno.clean$group[split_list$D$all$R8$folds$K1]))== 9
length(unique(Pheno.clean$group[split_list$D$all$R12$folds$K3]))== 9
length(unique(Pheno.clean$group[split_list$D$all$R20$folds$K1]))== 9
length(unique(Pheno.clean$group[split_list$D$all$R29$folds$K3]))== 9

#### require using subset
# 5. Scenario E: we should find 5 or 4 families depending of the group
length(unique(Pheno.clean.subsets[["F23"]]$group[split_list$E$F23$R1$folds$K1])) == 5
length(unique(Pheno.clean.subsets[["F23"]]$group[split_list$E$F23$R3$folds$K5])) == 5
length(unique(Pheno.clean.subsets[["F23"]]$group[split_list$E$F23$R20$folds$K4])) == 5

length(unique(Pheno.clean.subsets[["JFI"]]$group[split_list$E$JFI$R1$folds$K1])) == 4
length(unique(Pheno.clean.subsets[["JFI"]]$group[split_list$E$JFI$R3$folds$K5])) == 4
length(unique(Pheno.clean.subsets[["JFI"]]$group[split_list$E$JFI$R20$folds$K4])) == 4

# 6. Scenario F: unique one family
for (fam in levels(Pheno.clean$group)){
  print(length(unique(Pheno.clean.subsets[[fam]]$group[split_list[["F"]][[fam]]$R1$folds$K3])) == 1)
}

write_csv(Pheno.clean, '../../data/gs/pheno.csv')
save(split_list, file = "../../data/gs/splits.RData")
