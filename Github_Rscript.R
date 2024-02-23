# input, myabun: abundance table, n samples x m genomens
# input, mydata: metadata table, n samples x s variables

### Beta Diveristy Analysis ###
mybray=vegdist(myabun,method = "bray")
mypcoa=pcoa(mybray)
a=diag(t(mypcoa$vectors)%*%mypcoa$vectors)
explain=a/sum(a)*100
explain=round(explain,2)

### PERMANOVA test ###
set.seed(315) 
adonis2(as.matrix(mybray)~Group,mydata) #single factor
set.seed(315)
adonis2(as.matrix(mybray)~Age,mydata) # single factor
set.seed(315)
adonis2(as.matrix(mybray)~Group+Age,mydata,by="margin") # margin

### RDA analysis ###
myabun_hellinger=decostand(myabun,method = "hellinger")
decorana(myabun_hellinger)
tmpmod=rda(myabun_hellinger~Group,mydata)
set.seed(315)
anova(tmpmod)
RsquareAdj(tmpmod)
tmpgoodness=goodness(tmpmod)

### Random Forest ###
selectgroup=c("EDG","LDG")
ind=which(mygroup$Group%in%selectgroup)
mydata=mygroup[ind,]
mydata$Group=factor(mydata$Group)
myabun_54genome=myabun[ind,genome54]
tmptraindata=mydata
tmptrainabun=myabun_54genome
train.control=trainControl(method="LOOCV", savePredictions=TRUE,classProbs = TRUE)
set.seed(315)
registerDoMC(2) # parallel 
model=train(x=tmptrainabun_rel,y=tmptraindata$Group,method = "rf",trControl = train.control)

### Guild Index calculation ###
name1=mynamicmods$ID_raw[which(mynamicmods$group=="Guild 1")]
name2=mynamicmods$ID_raw[which(mynamicmods$group=="Guild 2")]
b1=apply(myabun[,name1],1,function(x){vegan::diversity(x,index = "simpson")}) # simpson index
b2=apply(myabun[,name2],1,function(x){vegan::diversity(x,index = "simpson")}) # simpson index
c1=apply(myabun[,name1],1,sum)
c2=apply(myabun[,name2],1,sum)
mygenome54_Index=data.frame(Guild1index=c1*100*b1,Guild2index=c2*100*b2)

### AUROC ###
selectgroup=c("EDG","LDG")
ind=which(mygroup$Group%in%selectgroup)
g=roc(predictor=mygenome54_Index$Guild1index[ind],response =mydata$Group[ind],auc = TRUE)
ci95=plot.roc(x=mydata$Group[ind],predictor = mygenome54_Index$Guild1index[ind],ci=TRUE,print.auc=TRUE)
g=roc(predictor=mygenome54_Index$Guild2index[ind],response =mydata$Group[ind],auc = TRUE)
ci95=plot.roc(x=mydata$Group[ind],predictor = mygenome54_Index$Guild2index[ind],ci=TRUE,print.auc=TRUE)
### AUPRC ###
library(PRROC)
mypredict=predict(mymodel$finalModel,myWang2020allabun,type = "prob")
myp="ESRD"
myn="CON"
yintercept=round(length(which(mypredict$group==myp))/nrow(mypredict),2)
sscurves=pr.curve(scores.class0 =mypredict$Case[mypredict$group==myp],
                  scores.class1 =mypredict$Case[mypredict$group==myn],curve=T)
AUPRC_draw=data.frame(sscurves$curve)