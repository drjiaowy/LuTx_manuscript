##Fig.1C & Fig.S3
JSD_total <- function(data){
  l <- length(data)
  matrix <- matrix(nrow = 1, ncol = l-1)  
  colnames(matrix) <- c(colnames(data)[-1])
  for (i in 1:(l-1)) {
    matrix[1,i] <- jensen_shannon(data[,i],data[,(i+1)])
  }
  return(matrix)
}

jensen_shannon <- function(p, q){
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  jsd = Hj - 0.5*(Hp+Hq)
  return(jsd)
}

shannon.entropy <- function(p){
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

##Fig.1D & Fig.S4A S4B
REP_total <- function(data, pre_donor, pre_reci){
  l <- length(data)
  Rep <- list()
  Rep[1] <- as.data.frame(pre_donor)
  Rep[2] <- as.data.frame(pre_reci)
  add <- c(Rep[[1]], Rep[[2]])
  for (i in 1:l) {
    Rep[i+2] <- as.data.frame(setdiff(rownames(data[data[,i]>0,]), add))
    add <- c(add, Rep[[i+2]])
  }
  matrix <- matrix(nrow = l, ncol = (l+2))  
  rownames(matrix) <- colnames(data)
  colnames(matrix) <- c(rep(c("pre donor", "pre reci", colnames(data)), 1))
  for (m in 1:(l+2)){
    for (n in 1:l) {
      matrix[n,m] <- sum(data[intersect(rownames(data[data[,n]>0,]), Rep[[m]]),n])/sum(data[data[,n]>0,n])*100
    }
  }
  return(matrix)
}

# e.g. Pt8
colnames(Pt8_new)
[1] "Pt8 R4U  ---  289102"              "Pt8 R8U  ---  191459"              "Pt8 R4L  ---  223318"              "Pt8 R8L  ---  33483"              
[5] "Pt8 D4U  ---  241749"              "Pt8 D8U  ---  249644"              "Pt8 D4L  ---  269183"              "Pt8 D8L  ---  219731"             
[9] "Pt8 POD0 R PBMC  ---  76722"       "Pt8 POD0 R Lu Percoll  ---  34072" "Pt8 POD0 D Lu Percoll  ---  26433" "Pt8 PBMC POD6  ---  88921"        
[13] "Pt8 PBMC POD15  ---  141289"       "Pt8 PBMC POD62  ---  512536"       "Pt8 PBMC POD167  ---  445799"      "Pt8 PBMC POD235  ---  141439"     
[17] "Pt8 PBMC POD356  ---  90838"       "Pt8 BAL POD15  ---  8165"          "Pt8 BAL POD56  ---  4853"          "Pt8 BAL POD91  ---  5145"         
[21] "Pt8 BAL POD154  ---  6937"         "Pt8 BAL POD259  ---  2023"         "Pt8 BAL POD294  ---  3286"         "Pt8 BAL POD364  ---  7116"        
[25] "Pt8 BW POD7  ---  3434"            "Pt8 Bx POD56  ---  461"            "Pt8 Bx POD91  ---  801"            "Pt8 Bx POD154  ---  138"          
[29] "Pt8 Bx POD294  ---  549"   

# e.g. Pt8 BAL
data <- Pt8_new[, 12:17]
pre_donor <- rownames(Pt8_new[Pt8_new[,5]>0|Pt8_new[,6]>0|Pt8_new[,7]>0|Pt8_new[,8]>0|Pt8_new[,11]>0, ])
pre_reci <- rownames(Pt8_new[Pt8_new[,1]>0|Pt8_new[,2]>0|Pt8_new[,3]>0|Pt8_new[,4]>0|Pt8_new[,9]>0|Pt8_new[,10]>0, ])
matrix <- REP_total(data, pre_donor, pre_reci)

##Fig.1E & Fig.S4C & Fig.S5
REP_remain_Fra_v2 <- function(pre_reci, data){
  data1 <- data[pre_reci,]  
  l <- length(data1)
  Rep <- list()
  Rep[[1]] <- rownames(data1[data1[,1]>0,])
  add <- Rep[[1]]
  for (i in 2:l) {
    Rep[i] <- as.data.frame(setdiff(rownames(data1[data1[,i]>0,]), add))
    add <- c(add, Rep[[i]])
  }
  
  matrix <- as.data.frame(matrix(nrow = l, ncol = l))  
  rownames(matrix) <- colnames(data) -> colnames(matrix) 
  
  for (m in 1:l){
    for (n in m:l) {
      matrix[m,n] <- length(which(data1[Rep[[m]], n]>0))/length(which(data1[Rep[[m]], m]>0))*100
    }
  }
  return(matrix)
}

# e.g. Pt8
CD4_reci <- rownames(Pt8_new[Pt8_new[,1]>0|Pt8_new[,3]>0, ])
CD8_reci <- rownames(Pt8_new[Pt8_new[,2]>0|Pt8_new[,4]>0, ])
BAL <- Pt8_new[,c(25, 18:24)]
PBMC <- Pt8_new[,c(12:17)]
Bx <- Pt8_new[,c(26:29)]
matrix <- REP_remain_Fra_v2(CD4_reci, BAL)

##Fig.2A & Fig.S6
Share_Lung_Fre <- function(data, Lym, PBMC, Lu, donor){
  matrix <- as.data.frame(matrix(nrow = 9, ncol = ncol(data)))
  rownames(matrix) <- c("Lym only", "PBMC only", "Lu only", "Lym+PBMC", "Lym+Lu", "PBMC+Lu", "Lym+PBMC+Lu", "unmappable", "donor")
  colnames(matrix) <- colnames(data)
  
  error <- intersect(donor, c(Lym, PBMC, Lu))
  Lym_only <- setdiff(Lym, c(PBMC, Lu, donor))
  PBMC_only <- setdiff(PBMC, c(Lym, Lu, donor))
  Lu_only <- setdiff(Lu, c(PBMC, Lym, donor))
  
  Lym_PBMC <- setdiff(intersect(Lym, PBMC), c(Lu, donor))
  Lym_Lu <- setdiff(intersect(Lym, Lu), c(PBMC, donor))
  PBMC_Lu <- setdiff(intersect(PBMC, Lu), c(Lym, donor))
  Lym_PBMC_Lu <- setdiff(intersect(Lym, intersect(PBMC, Lu)), donor)
  
  unmappable <- setdiff(rownames(data), c(Lym, PBMC, Lu, donor))
  donor <- setdiff(donor, error)
  
  for(m in 1:ncol(data)){
    matrix[1,m] <- sum(data[rownames(data) %in% Lym_only, m])
    matrix[2,m] <- sum(data[rownames(data) %in% PBMC_only, m])
    matrix[3,m] <- sum(data[rownames(data) %in% Lu_only, m])
    matrix[4,m] <- sum(data[rownames(data) %in% Lym_PBMC, m])   
    matrix[5,m] <- sum(data[rownames(data) %in% Lym_Lu, m])
    matrix[6,m] <- sum(data[rownames(data) %in% PBMC_Lu, m]) 
    matrix[7,m] <- sum(data[rownames(data) %in% Lym_PBMC_Lu, m])
    matrix[8,m] <- sum(data[rownames(data) %in% unmappable, m])
    matrix[9,m] <- sum(data[rownames(data) %in% donor, m])   
  }
  
  return(matrix)
}

# e.g. Pt8
data <- Pt8_new
Lym <- rownames(Pt8_new[Pt8_new[,1]>0|Pt8_new[,2]>0|Pt8_new[,3]>0|Pt8_new[,4]>0, ])
PBMC <- rownames(Pt8_new[Pt8_new[,9]>0, ])
Lu <- rownames(Pt8_new[Pt8_new[,10]>0, ])
donor <- rownames(Pt8_new[Pt8_new[,5]>0|Pt8_new[,6]>0|Pt8_new[,7]>0|Pt8_new[,8]>0|Pt8_new[,11]>0, ])
matrix <- Share_Lung_Fre(data, Lym, PBMC, Lu, donor)

##Fig.2B & Fig.2C
Share_Lung_detrat <- function(data, Lym, PBMC, Lu, donor){
  l <- ncol(data)
  matrix <- matrix(nrow = l, ncol = 38)  
  colnames(matrix) <- c( "Lym_only Pre-transplant clones detected", "Lym_only post-Tx clones detected", 
                         "PBMC_only Pre-transplant clones detected", "PBMC_only post-Tx clones detected",
                         "Lu_only Pre-transplant clones detected", "Lu_only post-Tx clones detected", 
                         "Lym_PBMC Pre-transplant clones detected", "Lym_PBMC post-Tx clones detected", 
                         "Lym_Lu Pre-transplant clones detected", "Lym_Lu post-Tx clones detected",
                         "PBMC_Lu Pre-transplant clones detected", "PBMC_Lu post-Tx clones detected", 
                         "Lym_PBMC_Lu Pre-transplant clones detected", "Lym_PBMC_Lu post-Tx clones detected", 
                         
                         "Lym_PBMC_Lu Lym_only Odds Ratio", "Lym_PBMC_Lu Lym_only P-value", 
                         "Lym_PBMC_Lu PBMC_only Odds Ratio", "Lym_PBMC_Lu PBMC_only P-value", 
                         "Lym_PBMC_Lu Lu_only Odds Ratio", "Lym_PBMC_Lu Lu_only P-value", 
                         
                         "Lym_PBMC Lym_only Odds Ratio", "Lym_PBMC Lym_only P-value",
                         "Lym_PBMC PBMC_only Odds Ratio", "Lym_PBMC PBMC_only P-value",
                         
                         "Lym_Lu Lym_only Odds Ratio", "Lym_Lu Lym_only P-value",
                         "Lym_Lu Lu_only Odds Ratio", "Lym_Lu Lu_only P-value",
                         
                         "PBMC_Lu PBMC_only Odds Ratio", "PBMC_Lu PBMC_only P-value",
                         "PBMC_Lu Lu_only Odds Ratio", "PBMC_Lu Lu_only P-value",
                         
                         "Lu_only Lym_only Odds Ratio", "Lu_only Lym_only P-value", 
                         "Lu_only PBMC_only Odds Ratio", "Lu_only PBMC_only P-value", 
                         "Lym_only PBMC_only Odds Ratio", "Lym_only PBMC_only P-value"
  )
  rownames(matrix) <- colnames(data)
  
  error <- intersect(donor, c(Lym, PBMC, Lu))
  
  Lym_only <- setdiff(Lym, c(PBMC, Lu, donor))
  PBMC_only <- setdiff(PBMC, c(Lym, Lu, donor))
  Lu_only <- setdiff(Lu, c(PBMC, Lym, donor))
  
  Lym_PBMC <- setdiff(intersect(Lym, PBMC), c(Lu, donor))
  Lym_Lu <- setdiff(intersect(Lym, Lu), c(PBMC, donor))
  PBMC_Lu <- setdiff(intersect(PBMC, Lu), c(Lym, donor))
  
  Lym_PBMC_Lu <- setdiff(intersect(Lym, intersect(PBMC, Lu)), donor)
  
  data <- data[!(rownames(data) %in% c(error, donor)),]
  
  matrix[,1] <- length(Lym_only)
  matrix[,3] <- length(PBMC_only)
  matrix[,5] <- length(Lu_only)
  matrix[,7] <- length(Lym_PBMC)
  matrix[,9] <- length(Lym_Lu)
  matrix[,11] <- length(PBMC_Lu)
  matrix[,13] <- length(Lym_PBMC_Lu)
  
  for (i in 1:l) {
    matrix[i,2] <- length(which(data[rownames(data) %in% Lym_only, i]>0))
    matrix[i,4] <- length(which(data[rownames(data) %in% PBMC_only, i]>0))
    matrix[i,6] <- length(which(data[rownames(data) %in% Lu_only, i]>0))
    matrix[i,8] <- length(which(data[rownames(data) %in% Lym_PBMC, i]>0))
    matrix[i,10] <- length(which(data[rownames(data) %in% Lym_Lu, i]>0))
    matrix[i,12] <- length(which(data[rownames(data) %in% PBMC_Lu, i]>0))
    matrix[i,14] <- length(which(data[rownames(data) %in% Lym_PBMC_Lu, i]>0))
    
    a <- fisher.test(matrix(c(matrix[i,14], matrix[i,13], matrix[i,2], matrix[i,1]),nrow=2))
    matrix[i,15] <- a[[3]]
    matrix[i,16] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,14], matrix[i,13], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,17] <- b[[3]]
    matrix[i,18] <- b[[1]]
    c <- fisher.test(matrix(c(matrix[i,14], matrix[i,13], matrix[i,6], matrix[i,5]),nrow=2))
    matrix[i,19] <- c[[3]]
    matrix[i,20] <- c[[1]]
    
    a <- fisher.test(matrix(c(matrix[i,8], matrix[i,7], matrix[i,2], matrix[i,1]),nrow=2))
    matrix[i,21] <- a[[3]]
    matrix[i,22] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,8], matrix[i,7], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,23] <- b[[3]]
    matrix[i,24] <- b[[1]]
    
    a <- fisher.test(matrix(c(matrix[i,10], matrix[i,9], matrix[i,2], matrix[i,1]),nrow=2))
    matrix[i,25] <- a[[3]]
    matrix[i,26] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,10], matrix[i,9], matrix[i,6], matrix[i,5]),nrow=2))
    matrix[i,27] <- b[[3]]
    matrix[i,28] <- b[[1]]
    
    a <- fisher.test(matrix(c(matrix[i,12], matrix[i,11], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,29] <- a[[3]]
    matrix[i,30] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,12], matrix[i,11], matrix[i,6], matrix[i,5]),nrow=2))
    matrix[i,31] <- b[[3]]
    matrix[i,32] <- b[[1]]
    
    a <- fisher.test(matrix(c(matrix[i,6], matrix[i,5], matrix[i,2], matrix[i,1]),nrow=2))
    matrix[i,33] <- a[[3]]
    matrix[i,34] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,6], matrix[i,5], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,35] <- b[[3]]
    matrix[i,36] <- b[[1]]
    c <- fisher.test(matrix(c(matrix[i,2], matrix[i,1], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,37] <- c[[3]]
    matrix[i,38] <- c[[1]]
  }
  return(matrix)
}

# e.g. Pt8
data <- Pt8_new
Lym <- rownames(Pt8_new[Pt8_new[,1]>0|Pt8_new[,2]>0|Pt8_new[,3]>0|Pt8_new[,4]>0, ])
PBMC <- rownames(Pt8_new[Pt8_new[,9]>0, ])
Lu <- rownames(Pt8_new[Pt8_new[,10]>0, ])
donor <- rownames(Pt8_new[Pt8_new[,5]>0|Pt8_new[,6]>0|Pt8_new[,7]>0|Pt8_new[,8]>0|Pt8_new[,11]>0, ])
matrix <- Share_Lung_detrat(data, Lym, PBMC, Lu, donor)

##Fig.2D & Fig.2E & Fig.3
error <- intersect(donor, c(Lym, PBMC, Lu))
Lym_only <- setdiff(Lym, c(PBMC, Lu, donor))
PBMC_only <- setdiff(PBMC, c(Lym, Lu, donor))
Lu_only <- setdiff(Lu, c(PBMC, Lym, donor))
Lym_PBMC <- setdiff(intersect(Lym, PBMC), c(Lu, donor))
Lym_Lu <- setdiff(intersect(Lym, Lu), c(PBMC, donor))
PBMC_Lu <- setdiff(intersect(PBMC, Lu), c(Lym, donor))
Lym_PBMC_Lu <- setdiff(intersect(Lym, intersect(PBMC, Lu)), donor)
unmappable <- setdiff(rownames(data), c(Lym, PBMC, Lu, donor))

rcd4 <-  Pt8_new[,c(1,3)]
rcd4 <-  normalize(rcd4)
CD4HVG <-  rownames(rcd4[rcd4[,2]>0.00002 & rcd4[,2] > rcd4[,1]*2,])
CD4NonHVG <- setdiff(rownames(Pt8_new[Pt8_new[,1]>0|Pt8_new[,3]>0,]), CD4HVG)
rcd8 <-  Pt8_new[,c(2,4)]
rcd8 <-  normalize(rcd8)
CD8HVG <-  rownames(rcd8[rcd8[,2]>0.00002 & rcd8[,2] > rcd8[,1]*2,])
CD8NonHVG <- setdiff(rownames(Pt8_new[Pt8_new[,2]>0|Pt8_new[,4]>0,]), CD8HVG)
dcd4 <-  Pt8_new[,c(5,7)]
dcd4 <-  normalize(dcd4)
CD4GVH <-  rownames(dcd4[dcd4[,2]>0.00002 & dcd4[,2] > dcd4[,1]*2,])
CD4NonGVH <- setdiff(rownames(Pt8_new[Pt8_new[,5]>0|Pt8_new[,7]>0,]), CD4GVH)
dcd8 <-  Pt8_new[,c(6,8)]
dcd8 <-  normalize(dcd8)
CD8GVH <-  rownames(dcd8[dcd8[,2]>0.00002 & dcd8[,2] > dcd8[,1]*2,])
CD8NonGVH <- setdiff(rownames(Pt8_new[Pt8_new[,6]>0|Pt8_new[,8]>0,]), CD8GVH)

Lym_only_data <- data[rownames(data) %in% Lym_only, ]

# e.g. Pt8 BAL Lym_only
BAL <- Pt8_Lym_only[c(8, 9, 12, 14), c(8, 9, 12, 14)]*100
colnames(BAL) <- c("POD15", "POD56", "POD259", "POD364") -> rownames(BAL)
BAL$Name <- colnames(BAL)
melt <- melt(BAL, id.vars = c("Name"), variable.name="YName",value.name="overlap")
mycolors <- c("#99CC00","#9900CC","#996600","#0033FF")
melt$color <- c(mycolors)
melt <- melt[!(melt$Name == melt$YName),]
melt$angle <- c(3.14/12*1, 3.14/12*3, 3.14/12*5, 3.14/12*7, 
                3.14/12*9, 3.14/12*11, 3.14/12*13, 3.14/12*15,
                3.14/12*17, 3.14/12*19, 3.14/12*21, 3.14/12*23)

magnitude <- c(melt$overlap)
angle <- c(melt$angle)
directionlabels <- c(1:4)
colors <- c(melt$color)
par(cex.axis=0.7)
par(cex.lab=1)
par(mar=c(0,0,0,0))

radial.plot(c(magnitude), 
            c(angle), 
            lwd=4, 
            show.grid.label=T, #put the concentric circle labels going down
            show.radial.grid=F,
            line.col=colors,
            labels=directionlabels,
            radial.lim = c(0,100), #range of grid circle
            #main="circular diagram!",
            label.prop = 0,
            rad.col = c("white")
)

##Fig.4 & Fig.S7 & Fig.5 & Fig.S8
# e.g. Pt8 path
matrix <- as.data.frame(matrix(nrow=9, ncol=24))
rownames(matrix) <- c("Lym_only","PBMC_only","Lu_only","Lym_PBMC","Lym_Lu","PBMC_Lu","Lym_PBMC_Lu","unmappable","donor")
colnames(matrix) <- c("CD4HVG", "CD4NonHVG", "CD4GVH", "CD4NonGVH", "CD8HVG", "CD8NonHVG", "CD8GVH", "CD8NonGVH", 
                      "CMV", "EA", "EBV", "EC", "HCV", "HIV", "Influenza", "KP", "MTB", "SA", "Salmonella", "Shi",                       
                      "YFV", "COV", "total path", "total path but COV")

list <- list(Lym_only, PBMC_only, Lu_only, Lym_PBMC, Lym_Lu, PBMC_Lu, Lym_PBMC_Lu, unmappable, donor)

for(n in 1:9){
  matrix[n,1] <- length(intersect(list[[n]], CD4HVG))/length(list[[n]])*100
  matrix[n,2] <- length(intersect(list[[n]], CD4NonHVG))/length(list[[n]])*100
  matrix[n,3] <- length(intersect(list[[n]], CD4GVH))/length(list[[n]])*100
  matrix[n,4] <- length(intersect(list[[n]], CD4NonGVH))/length(list[[n]])*100
  matrix[n,5] <- length(intersect(list[[n]], CD8HVG))/length(list[[n]])*100
  matrix[n,6] <- length(intersect(list[[n]], CD8NonHVG))/length(list[[n]])*100
  matrix[n,7] <- length(intersect(list[[n]], CD8GVH))/length(list[[n]])*100
  matrix[n,8] <- length(intersect(list[[n]], CD8NonGVH))/length(list[[n]])*100
  
  matrix[n,9] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$CMV)]))/length(list[[n]])*100
  matrix[n,10] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$EA)]))/length(list[[n]])*100
  matrix[n,11] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$EBV)]))/length(list[[n]])*100
  matrix[n,12] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$EC)]))/length(list[[n]])*100
  matrix[n,13] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$HCV)]))/length(list[[n]])*100
  matrix[n,14] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$HIV)]))/length(list[[n]])*100
  matrix[n,15] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$Influenza)]))/length(list[[n]])*100
  matrix[n,16] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$KP)]))/length(list[[n]])*100
  matrix[n,17] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$MTB)]))/length(list[[n]])*100
  matrix[n,18] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$SA)]))/length(list[[n]])*100
  matrix[n,19] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$Salmonella)]))/length(list[[n]])*100
  matrix[n,20] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$Shi)]))/length(list[[n]])*100
  matrix[n,21] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$YFV)]))/length(list[[n]])*100
  matrix[n,22] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$COV)]))/length(list[[n]])*100
  matrix[n,23] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                   !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                   !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                   !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                   !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                   !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                   !is.na(Pt8_path$YFV)| !is.na(Pt8_path$COV)] ))/length(list[[n]])*100
  matrix[n,24] <- length(intersect(list[[n]], rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                   !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                   !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                   !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                   !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                   !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                   !is.na(Pt8_path$YFV)  ] ))/length(list[[n]])*100
}

matrix <- as.data.frame(matrix(nrow=ncol(Pt8_new), ncol=24))
rownames(matrix) <- colnames(Pt8_new)
colnames(matrix) <- c("CD4HVG", "CD4NonHVG", "CD4GVH", "CD4NonGVH", "CD8HVG", "CD8NonHVG", "CD8GVH", "CD8NonGVH", 
                      "CMV", "EA", "EBV", "EC", "HCV", "HIV", "Influenza", "KP", "MTB", "SA", "Salmonella", "Shi",                       
                      "YFV", "COV", "total path", "total path but COV")

for(n in 1:ncol(Pt8_new)){
  matrix[n,1] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD4HVG, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,2] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD4NonHVG, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,3] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD4GVH, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,4] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD4NonGVH, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,5] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD8HVG, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,6] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD8NonHVG, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,7] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD8GVH, n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,8] <- length(which(Pt8_new[rownames(Pt8_new) %in% CD8NonGVH, n]>0))/length(which(Pt8_new[, n]>0))*100
  
  matrix[n,9] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,10] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EA)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,11] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EBV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,12] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EC)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,13] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$HCV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,14] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$HIV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,15] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Influenza)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,16] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$KP)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,17] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$MTB)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,18] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$SA)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,19] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Salmonella)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,20] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Shi)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,21] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$YFV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,22] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$COV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,23] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                                   !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                                   !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                                   !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                                   !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                                   !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                                   !is.na(Pt8_path$YFV)| !is.na(Pt8_path$COV)], n]>0))/length(which(Pt8_new[, n]>0))*100
  matrix[n,24] <- length(which(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                                   !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                                   !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                                   !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                                   !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                                   !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                                   !is.na(Pt8_path$YFV)], n]>0))/length(which(Pt8_new[, n]>0))*100
}

for(n in 1:ncol(Pt8_new)){
  matrix[n,1] <- sum(Pt8_new[rownames(Pt8_new) %in% CD4HVG, n])/sum(Pt8_new[, n])*100
  matrix[n,2] <- sum(Pt8_new[rownames(Pt8_new) %in% CD4NonHVG, n])/sum(Pt8_new[, n])*100
  matrix[n,3] <- sum(Pt8_new[rownames(Pt8_new) %in% CD4GVH, n])/sum(Pt8_new[, n])*100
  matrix[n,4] <- sum(Pt8_new[rownames(Pt8_new) %in% CD4NonGVH, n])/sum(Pt8_new[, n])*100
  matrix[n,5] <- sum(Pt8_new[rownames(Pt8_new) %in% CD8HVG, n])/sum(Pt8_new[, n])*100
  matrix[n,6] <- sum(Pt8_new[rownames(Pt8_new) %in% CD8NonHVG, n])/sum(Pt8_new[, n])*100
  matrix[n,7] <- sum(Pt8_new[rownames(Pt8_new) %in% CD8GVH, n])/sum(Pt8_new[, n])*100
  matrix[n,8] <- sum(Pt8_new[rownames(Pt8_new) %in% CD8NonGVH, n])/sum(Pt8_new[, n])*100
  
  matrix[n,9] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)], n])/sum(Pt8_new[, n])*100
  matrix[n,10] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EA)], n])/sum(Pt8_new[, n])*100
  matrix[n,11] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EBV)], n])/sum(Pt8_new[, n])*100
  matrix[n,12] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$EC)], n])/sum(Pt8_new[, n])*100
  matrix[n,13] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$HCV)], n])/sum(Pt8_new[, n])*100
  matrix[n,14] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$HIV)], n])/sum(Pt8_new[, n])*100
  matrix[n,15] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Influenza)], n])/sum(Pt8_new[, n])*100
  matrix[n,16] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$KP)], n])/sum(Pt8_new[, n])*100
  matrix[n,17] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$MTB)], n])/sum(Pt8_new[, n])*100
  matrix[n,18] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$SA)], n])/sum(Pt8_new[, n])*100
  matrix[n,19] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Salmonella)], n])/sum(Pt8_new[, n])*100
  matrix[n,20] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$Shi)], n])/sum(Pt8_new[, n])*100
  matrix[n,21] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$YFV)], n])/sum(Pt8_new[, n])*100
  matrix[n,22] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$COV)], n])/sum(Pt8_new[, n])*100
  matrix[n,23] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                          !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                          !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                          !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                          !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                          !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                          !is.na(Pt8_path$YFV)| !is.na(Pt8_path$COV)], n])/sum(Pt8_new[, n])*100
  matrix[n,24] <- sum(Pt8_new[rownames(Pt8_new) %in% rownames(Pt8_path)[!is.na(Pt8_path$CMV)| !is.na(Pt8_path$EA)| 
                                                                          !is.na(Pt8_path$EBV)| !is.na(Pt8_path$EC)| 
                                                                          !is.na(Pt8_path$HCV)| !is.na(Pt8_path$HIV)| 
                                                                          !is.na(Pt8_path$Influenza)| !is.na(Pt8_path$KP)| 
                                                                          !is.na(Pt8_path$MTB)| !is.na(Pt8_path$SA)| 
                                                                          !is.na(Pt8_path$Salmonella)| !is.na(Pt8_path$Shi)|
                                                                          !is.na(Pt8_path$YFV)], n])/sum(Pt8_new[, n])*100
}



