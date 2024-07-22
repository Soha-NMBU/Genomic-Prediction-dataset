#install pls package (if not already installed)
#install.packages("pls")
#load pls package
rm(list = ls())
library(pls)
library(SKM)
library(BGLR)
#make this example reproducible
set.seed(1)
Summary_All=data.frame()
Summary_All=data.frame()

data <- load("Norway_data_2022.RData", verbose = TRUE)
ls()
Data_names=list("EC1","EC2","EC3","EC4","EC5","EC6")
Data_names[1]
for(w in 1:6){
#w=1
name <-Data_names[w]

Pheno=Pheno
head(Pheno)
Gids=sort(Markers1$Line)
Gids
Markers_Raw=Markers1[Gids,]
Markers=Markers_Raw[,-1]
Pos_var=which(apply(Markers,2,var)>0)
#length(Pos_var)
Markers_Scaled=scale(Markers[,Pos_var])

Geno=Markers_Scaled%*%t(Markers_Scaled)/ncol(Markers_Scaled)

Pheno = Pheno[order(Pheno$Env, Pheno$Line),] #ordenar pheno por environment
rownames(Pheno) =1:nrow(Pheno)
#head(Pheno)
dim(Pheno)
dim(Geno)

# Data preparation
Line <- model.matrix(~0 + as.factor(Line), data = Pheno)
dim(Line)
X_Env <- model.matrix(~0 + Env, data = Pheno)
dim(X_Env)
colnames(X_Env)
#####Environmenal covariates
EC<-EC_ls[[w]]
head(EC[,1:6])
EC_Ord= EC[order(EC$Env, EC$Line),]
EC=as.matrix(EC_Ord[,-c(1:2)])
head(EC[,1:6])
dim(EC)
EC_Scaled=scale(EC)
X_Env=X_Env
dim(X_Env)
K.E=X_Env%*%t(X_Env)/ncol(X_Env)
Geno1=Geno
diag(Geno1)=diag(Geno1)+0.001
L_G <- t(chol(Geno1))
X_Line<- Line%*%L_G
Geno11=as.matrix(Geno1)
K.L=Line%*%Geno11%*%t(Line)
K.LE=K.L*K.E
X_LE=t(chol(K.LE))
X=cbind(X_Env,X_Line)
head(Pheno)
Trait_names=colnames(Pheno)[4:ncol(Pheno)]
Trait_names
All_summary=data.frame()
  
for (t in 1:2){
 #t=4
 Trait_t=Trait_names[t]
 y<- Pheno[, 3+t]
y[4]=NA
 pos_Na=which(is.na(y)==TRUE)
 pos_Na
if (length(pos_Na)>0) {
    y[ pos_Na]=median(y)
} else {
  y=y
}
  Envs=unique(Pheno$Env)
  Predictions <- data.frame()
  for (e in 1:length(Envs)){
    cat("*** Fold:", e, " ***\n")
    #  e=1
    tst_e=Envs[e]
    pos_test_e=which(Pheno$Env==tst_e)
    ETA=list(A=list(K=K.E , model="RKHS"),B=list(K=K.L, model="RKHS"), D=list(X=EC_Scaled, model="BRR"))
    y2=y
    y2[pos_test_e]=NA
    model <-BGLR(
      ETA=ETA,
      y=y2,
      nIter =10000,
      burnIn= 2500,
      verbose = F)
    y_tst=y[pos_test_e]
    
    ###summary(model_Final, what = "validation")
    Y_predicted <-model$yHat[pos_test_e] 
    Predictions_e=data.frame(Line=Pheno$Line[pos_test_e], Fold=e, Env=Pheno$Env[pos_test_e], Observed=y_tst,Predicted=c(Y_predicted))
    Predictions=rbind(Predictions,Predictions_e)
  }
  Predictions
  Pred_Env_Summary=gs_summaries(Predictions)$env
  Pred_Env_Summary
  All_summary=rbind(All_summary,data.frame(Trait_t=Trait_t, Pred_Env_Summary))
}

Saving_name=paste(name, "Pred_BGLR_CVO_Traits_E+G.csv", sep="_")
write.csv(All_summary,file=Saving_name[1])
}

