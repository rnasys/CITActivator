###
#   SPLSDA cross validation accuracy:
#   sparisy removal
###
### 2 class : Activator, Others
# load packages
library("sda")
library("tm")
library("crossval")
library("measures")
library("mixOmics")
library("qwraps2")
library("pROC")
library("caret")

# load data training & test data
data <- read.table("110RBP_amino_seq_feature_matrix_total_length.txt",sep = "\t",header = T,row.names = 1)
dim(data)
TE <- read.table("110RBP.txt",sep = "\t",header = T,row.names = 1)
df <- merge(data,TE,by="row.names")
rownames(df)<- df$Row.names
data_TE <- df$TE_AVE
df <- df[,-1]
data <- df[,c(1:ncol(data))]
data_Class <- ifelse(data_TE>1.3,"Activator","Neutral & Inhibitor")
summary(factor(data_Class))
X <- data
Y <- factor(data_Class)

list.keepX <- c(seq(50, 300, 25))
accuracy <- c()
act_accuracy <- c()
oth_accuracy <- c()
pred.splsda = function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=FALSE)
{ 
  myTdm1 <- as.DocumentTermMatrix(slam::as.simple_triplet_matrix(Xtrain), weighting = weightTf)
  dim(as.matrix(removeSparseTerms(myTdm1,0.95)))
  Xtrain <- as.matrix(removeSparseTerms(myTdm1,0.95))
  Ytrain <- Ytrain
  #dim(X)
  tune.splsda <- tune.splsda(Xtrain, Ytrain, ncomp = 5, validation = 'Mfold', folds = 10, 
                             progressBar = TRUE, dist = 'centroids.dist',
                             test.keepX = list.keepX, nrepeat = 3, cpus = 2)
  error <- tune.splsda$error.rate  # error rate per component for the keepX grid
  ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
  ncomp
  select.keepX <- tune.splsda$choice.keepX[1:ncomp] 
  freq.splsda <- mixOmics::splsda(Xtrain, Ytrain, ncomp = 2,keepX = select.keepX)  # set ncomp to 10 for performance assessment later
  selVars <- selectVar(freq.splsda)$name
  
  freq.splsda <- mixOmics::splsda(Xtrain[,selVars], Ytrain, ncomp = 10) 
  plotIndiv(freq.splsda , comp = 1:2,
            group = Ytrain, ind.names = T, ellipse.level = 0.9,cex = 3,
            ellipse = T,legend = TRUE, title = 'PLSDA',col = c("#F08080","grey"))
  plotVar(freq.splsda, var.names=T,cutoff = 0.3,cex = 3.5) 
  
  test.predict <- predict(freq.splsda, Xtest[,selVars], dist = "centroids.dist")
  ynew = test.predict$class$centroids.dist[,2]
  ynew_prob <- as.numeric(test.predict$predict[,1,10])
  x <- ifelse(as.character(ynew) == "Activator",1,0)
  y <- ifelse(as.character(Ytest) == "Activator",1,0)
  print(x)
  print(y)
  
  print(cft <- table(x, y))
  if((nrow(as.matrix(cft))[1]==2) & (ncol(as.matrix(cft))==2)){
    tp <- cft[2, 2]
    tn <- cft[1, 1]
    fp <- cft[2, 1]
    fn <- cft[1, 2]
    cm <- confusionMatrix(cft, positive = "1")
    #print(cm)
    # performance
    Accuracy <- (tp + tn)/(tp + tn + fp + fn)
    Sensitivity <- tp/(tp + fn)
    Specificity <- tn/(tn + fp)
    Precision <- tp/(tp+fp) 
    others_Precision <- tn/(tn+fn) 
    Recall <- as.numeric(cm$byClass[6])
    F1 <- as.numeric(cm$byClass[7])
    Balanced_Accuracy <- as.numeric(cm$byClass[11])
    roc_obj <- roc(y,ynew_prob)
    AUC <- auc(roc_obj)
    print(Accuracy)
    print(Precision)
    print(others_Precision)
    print(F1)
    print(Balanced_Accuracy)
    print(AUC)
    return(c(Accuracy,Precision,others_Precision,F1,Balanced_Accuracy,AUC))
  }else{
    return(rep(NA,6))
  }
}
cv.lda.plsda = crossval(pred.splsda, X, Y, K=7, B=3)
cv.lda.plsda$stat
colMeans(na.omit(cv.lda.plsda$stat.cv))  
write.table(cv.lda.plsda$stat.cv,"twoclass_cv_1.3_0.95_7.results",quote = F)


df <- read.table("twoclass_cv_1.3_0.95_7.results",sep=" ",row.names = 1)
df <- na.omit(df)
df_ggplot <- data.frame(acc = as.numeric(c(df[,1],df[,2],df[,3],df[,4])),
                        class = c(rep("ACC: Classification Accuracy",nrow(df)),
                                  rep("PPV: Positive Predictive Value",nrow(df)),
                                  rep("NPV: Negative Predictive Value",nrow(df)),
                                  rep("AUC: Area Under Curve",nrow(df))))
library("ggplot2")
library("ggthemes")
library("Rmisc")
df_ggplot_sum <- summarySE(df_ggplot, measurevar="acc", groupvars=c("class"))
write.table(df_ggplot_sum,"comp10_twoclass_cv_1.3.summary",sep = "\t",quote = F)

pdf("Evaluation_7fold_1.3.pdf",height = 4,width = 6)
ggplot(df_ggplot_sum,aes(x = class,y = acc,fill = class))+#geom_boxplot()+
  geom_bar(position = position_dodge(),stat="identity",width = 0.45,color = "grey20",alpha= 0.6)+
  geom_errorbar(aes(ymin=acc, ymax=acc+sd),
                width=.15,               # Width of the error bars
                position=position_dodge(.9),
                color = "grey20") + 
  ylab("")+
  ggtitle("Model Evaluation") + 
  theme_few() +
  scale_x_discrete(limits = c("ACC: Classification Accuracy",
                              "PPV: Positive Predictive Value",
                              "NPV: Negative Predictive Value",
                              "AUC: Area Under Curve"),
                   labels = c("ACC","PPV","NPV","AUC"))+ 
  scale_fill_manual(values = c("#E2AB7F","#C05850","grey60","#527590"),
                    limits = c("ACC: Classification Accuracy",
                               "PPV: Positive Predictive Value",
                               "NPV: Negative Predictive Value",
                               "AUC: Area Under Curve"))+
  theme(axis.text.x = element_text(size=10,angle=0,colour = "black",face="bold",hjust = 0.5,vjust = 1),
        axis.text.y = element_text(size=11,angle=0,colour = "black",hjust = 1,vjust = 1),
        plot.title = element_text(lineheight=50, face="bold",size = 13,hjust = 0.5),
        axis.title.y = element_text(size=12,face="bold",color="black"),
        axis.title.x = element_text(size=10,color="white"),
        legend.text = element_text(size=11,color="black"),
        #legend.position= c(0.25,0.9),
        legend.title = element_text(size=10,color="white"))
dev.off()