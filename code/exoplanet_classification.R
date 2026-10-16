library(survival)
library(plyr)
library(dplyr)
library(funModeling)
library(VIM)
library(mctest)
library(mice)
library(MASS)
library(lmtest)
library(gam)
library(car) 
library(sandwich)
library(robustbase)
library(robust)
require(corrgram)
library(caret)
library(sampling)
library(rpart)
library(rpart.plot)
library(missForest)
library(fastAdaboost)
library(partykit)
library(pls)
library(gbm)
library(nnet)
library(pROC)
library(ROCR)

DATA_PATH <- file.path("data", "esopianeti.csv")
b <- read.csv(DATA_PATH, sep=",",
              dec = ".",
              stringsAsFactors=TRUE, na.strings=c("NA","NaN", ""))

##PREPROCESSING
#Distribuzioni di alcune covariate
table(b$sy_snum)
table(b$sy_pnum)
table(b$disc_telescope)
table(b$disc_facility)

#Distribuzioni dalla target
table(b$ttv_flag)
table(b$ttv_flag)/nrow(b)

# Si tolgono gli identificativi e formati rastr e decstr perch� sono in formato sessagesimale. Gli ultimi due sono gli intervalli di 
#confidenza della misurazione della Gaia Magnitude
b = b[, -c(1:2, 27, 29, 35:36)]

# analisi preliminare: zero variance e factor
b$ttv_flag =as.factor(b$ttv_flag)
isfactor <- sapply(b, function(x) is.factor(x))
isfactor

# selezione delle variabili categoriali
factordata <- b[, isfactor]
str(factordata)
sapply(factordata, function(x) length(levels(x)))
b = b[, -c(3)]

# selezione delle variabili continue
numeric <- sapply(b, function(x) is.numeric(x))
numeric <-b[, numeric]
str(numeric)

# visualizzazione delle variabili quantitative per osservarne la varianza
sapply(numeric, function(x) var(x,na.rm=TRUE))

## MISSING IMPUTATION
# visualizzazione del numero di missing per ogni variabile
sapply(b, function(x)(sum(is.na(x))))
status=df_status(b, print_results = F)
status 
sapply(b, function(x)(sum(is.na(x))))
status<-df_status(b, print_results = F)
status %>% arrange(-q_na)

# si tolgono le colonne con na > 20% e il target
b = b[, -c(13:15, 17:19)]

# si ottiene un vettore con solo le variabili categoriali
isfactor <- sapply(b, function(x) is.factor(x))
factordata <- b[, isfactor]
str(factordata)

# si ottiene un vettore con solo le variabili continue
numeric <- sapply(b, function(x) is.numeric(x))
numeric <-b[, numeric]
str(numeric)

# ottengo una lista di covariate
Formula <- paste(colnames(factordata), collapse="','")
Formula
Formula1 <- paste(colnames(numeric), collapse="','")
Formula1
list=c('disc_telescope','pl_orbper','pl_dens','sy_dist', 'sy_snum','sy_pnum','pl_orbper','pl_orbsmax','pl_rade','pl_radj','pl_bmasse','pl_bmassj','pl_dens','pl_orbeccen','ttv_flag','st_teff','st_mass','st_lum','st_logg','ra','dec','sy_dist','sy_vmag','sy_kmag','sy_gaiamag')
variabili = b[,list]
head(variabili)

# rappresentazione grafica dei dati mancanti
mancanti <- aggr(variabili, col=c('navyblue', 'yellow'), numeber = FALSE, sortvars=TRUE, labels = names(variabili), cex.axis=1, gap=1)

# Si toglie la colonna disc_telescope, perch� � sottogruppo della variabile disc_facility
set.seed(1)  
grup.tree <- rpart(ttv_flag ~ factor(disc_telescope), cp=0, data = b, method = "class")
grup.tree$cptable
plot(grup.tree)

library(partykit)
z <- as.party(grup.tree)
b$node=predict(z, b, type = "node")
head(b)

aggr = b %>% group_by(node) %>% summarise(
  n = dplyr::n(),
  prob = mean(as.numeric(as.character(ttv_flag)))
)
aggr
b$disc_telescope = NULL

#missing imputation con random forest
set.seed(196) 
b2 = b[,-c(12)]
ds_imputato <- missForest(b2, maxiter = 10, ntree = 30) 
imputedRF=data.frame(ds_imputato$ximp) 
data_imputedRF <- cbind(imputedRF, b$ttv_flag)
names(data_imputedRF)[names(data_imputedRF) == 'b$ttv_flag'] <- 'ttv_flag'
sapply(imputedRF, function(x)(sum(is.na(x)))) # NA counts
status<-df_status(imputedRF, print_results = F) #perc NA
status %>% arrange(-q_na)

##COLLINEARITA'
# si costruisce un dataset con solo covariate continue
numeric2 <- sapply(data_imputedRF, function(x) is.numeric(x))
numeric2 <-data_imputedRF[, numeric2]
Formula2 <- paste(colnames(numeric2), collapse="'+'")
Formula2

# si calcolano i coefficienti di correlazione lineare
rho_num = cor(numeric2)
rho_num

# matrice di correlazione
corrgram(numeric2, lower.panel = panel.cor, cex=1, cex.labels = 1)
cor(data_imputedRF$pl_rade, data_imputedRF$pl_radj)
cor(data_imputedRF$pl_bmasse, data_imputedRF$pl_bmassj)
cor(data_imputedRF$sy_vmag, data_imputedRF$sy_kmag)
cor(data_imputedRF$sy_kmag, data_imputedRF$sy_gaiamag)
cor(data_imputedRF$sy_vmag, data_imputedRF$sy_gaiamag)

# vengono tolte le variabili collineari e non di interesse
data_imputedRF$pl_radj <- NULL
data_imputedRF$pl_bmassj <- NULL
data_imputedRF$sy_kmag <- NULL
data_imputedRF$sy_gaiamag <- NULL

## SEPARAZIONE DEL DATASET: TRAINING, VALIDATION, SCORE
# undersampling con le vere prior (12,6)
df <- data_imputedRF
set.seed(1000)
id <- strata(df, stratanames = "ttv_flag", size = c(2400, 350), method = "srswr")
df_under <- df |> slice(id$ID_unit)
table(df_under$ttv_flag)

# suddivisione del dataset
fractionTraining   <- 0.60
fractionValidation <- 0.30
fractionScore      <- 0.10
sampleSizeTraining   <- floor(fractionTraining   * nrow(df_under))
sampleSizeValidation <- floor(fractionValidation * nrow(df_under))
sampleSizeScore      <- floor(fractionScore      * nrow(df_under))
indicesTraining    <- sort(sample(seq_len(nrow(df_under)), size=sampleSizeTraining))
indicesNotTraining <- setdiff(seq_len(nrow(df_under)), indicesTraining)
indicesValidation  <- sort(sample(indicesNotTraining, size=sampleSizeValidation))
indicesScore       <- setdiff(indicesNotTraining, indicesValidation)
Training   <- df_under[indicesTraining, ]
Validation <- df_under[indicesValidation, ]
Score      <- df_under[indicesScore, ]

# controllo dei dataset
table(data_imputedRF$ttv_flag)/nrow(data_imputedRF)
table(Training$ttv_flag)/nrow(Training)
table(Validation$ttv_flag)/nrow(Validation)
table(Score$ttv_flag)/nrow(Score)

##STEP 1: TUNING DEI MODELLI
# ricodifica del target
Training$u=ifelse(Training$ttv_flag==1,"yes","no")
Training$ttv_flag=NULL
Validation$u=ifelse(Validation$ttv_flag==1,"yes","no")
Validation$ttv_flag=NULL
Score$u=ifelse(Score$ttv_flag==1,"yes","no")
Score$ttv_flag=NULL

# albero decisionale
set.seed(1)
metric <- "Sens"
Ctrl <- trainControl(method = "cv" , number=10, summaryFunction = twoClassSummary, classProbs = TRUE) # this gives ROC as measure, without this is accuracy
rpartTune2 <- train(u ~ ., data = Training, method = "rpart", 
                    tuneLength = 15, trControl = Ctrl, minsplit = 5, metric=metric)
rpartTune2

# tuning con la metrica ROC
metric <- "ROC"
set.seed(1)
cvCtrl <- trainControl(method = "cv", number=10, search="grid", classProbs = TRUE,
                       summaryFunction = twoClassSummary)
tree <- train(u ~ ., data = Training, method = "rpart",
              tuneLength = 10, metric=metric,
              trControl = cvCtrl)
tree

# visualizzazione comportamento della ROC rispetto alla complessit� dell'albero
plot(tree)

# matrice di confusione dell'albero decisionale
getTrainPerf(tree)
confusionMatrix(tree)

# importanza e selezione delle variabili nell'albero
plot(varImp(object=tree),main="train tuned - Variable Importance")
vi=as.data.frame(tree$finalModel$variable.importance)
viname=row.names(vi)
head(viname)

# dataset di training, validation e score con le covariate selezionate dall'albero
Train_modsel=Training[,c("dec","ra","node","pl_orbper", "pl_orbsmax","pl_rade","sy_pnum","sy_vmag", "pl_bmasse")]
Train_modsel=cbind(Training$u, Train_modsel)
head(Train_modsel)
names(Train_modsel)[1] <- "u"
head(Train_modsel)
Valid_modsel=Validation[,c("dec","ra","node","pl_orbper", "pl_orbsmax","pl_rade","sy_pnum","sy_vmag", "pl_bmasse")]
Valid_modsel=cbind(Validation$u, Valid_modsel)
head(Valid_modsel)
names(Valid_modsel)[1] <- "u"
head(Valid_modsel)
Score_modsel=Score[,c("dec","ra","node","pl_orbper", "pl_orbsmax","pl_rade","sy_pnum","sy_vmag", "pl_bmasse")]
Score_modsel=cbind(Score$u, Score_modsel)
head(Score_modsel)
names(Score_modsel)[1] <- "u"
head(Score_modsel)
Data_modsel=rbind(Train_modsel,Valid_modsel, Score_modsel)

# random forest
seed <- 7
set.seed(seed)
metric <- "ROC"
control <- trainControl(method="cv", number=10, search="grid", summaryFunction = twoClassSummary, classProbs = TRUE)
tunegrid <- expand.grid(.mtry=c(1:5))
rf <- train(u~., data=Training, method="rf", metric=metric, tuneGrid=tunegrid, ntree=250, trControl=control)
rf

# grafico del ROC in base al numero di covariate prese in modo casuale
plot(rf)

# oggetto nel random forest
ls(rf)
# performance dato mtry
rf$results
# modello finale e le sue performance
getTrainPerf(rf)
rf$resample
# media crossvalidata con 10 fold
mean(rf$resample$ROC)
mean(rf$resample$Spec)
z=data.frame(rf$resample)
head(z)
# matrice cross validata
confusionMatrix(rf)

# importanza e selezione delle variabili in una random forest 
Vimportance <- varImp(rf)
head(Vimportance)
Vimportance
VImP=as.data.frame(Vimportance$importance)
head(VImP)
dim(VImP)
V=subset(VImP, Overall>15)
V2=t(V)
row.names(V) 

# dataset di training, validation e score con le covariate selezionate dal random forest
u=Training[19] 
u_v=Validation[19] 
Xselected=Training[,c("sy_pnum", "pl_orbper","pl_orbsmax","pl_rade","pl_bmasse","pl_bmassprov", "pl_dens", "pl_orbeccen", "st_teff", "st_mass", "st_lum", "st_logg", "ra", "dec", "sy_dist", "sy_vmag", "node")]
Xselected_v=Validation[,c("sy_pnum", "pl_orbper","pl_orbsmax","pl_rade","pl_bmasse", "pl_bmassprov", "pl_dens", "pl_orbeccen", "st_teff", "st_mass", "st_lum", "st_logg", "ra","dec","sy_dist","sy_vmag", "node" )]
train_rf=cbind(Xselected,u) #newDAta
head(train_rf)
valid_rf=cbind(Xselected_v,u_v)
colnames(train_rf)
colnames(valid_rf)
row.names(V)

# modello logistico
set.seed(1)
metric <- "ROC"
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary)
glm=train(u~.,
          data=train_rf,method = "glm",
          trControl = ctrl, tuneLength=5, trace=TRUE, na.action = na.pass, metric=metric)
glm
confusionMatrix(glm)

#Lasso
set.seed(1)
metric <- "ROC"
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary)
grid = expand.grid(.alpha=1,.lambda=seq(0, 1, by = 0.01))
lasso=train(u~.,
            data=Training,method = "glmnet",
            trControl = ctrl, tuneLength=5, na.action = na.pass,metric=metric,
            tuneGrid=grid)
lasso
plot(lasso)
confusionMatrix(lasso)

#Naive Bayes con il dataset ricavato con l'albero
set.seed(1)
metric <- "ROC"
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary)
naivebayes_tree=train(u~.,
                      data=Train_modsel,method = "naive_bayes", metric=metric,
                      trControl = ctrl, tuneLength=5, na.action = na.pass) 
naivebayes_tree
confusionMatrix(naivebayes_tree)

#Naive Bayes con il dataset ricavato dal Random Forest
set.seed(1)
metric <- "ROC"
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary)
naivebayes_rf=train(u~.,
                    data=train_rf,method = "naive_bayes", metric=metric,
                    trControl = ctrl, tuneLength=5, na.action = na.pass) 
naivebayes_rf
confusionMatrix(naivebayes_rf)

#pls
set.seed(1234)
metric <- "ROC"
Control=trainControl(method= "cv",number=10, classProbs=TRUE,
                     summaryFunction=twoClassSummary)
pls=train(u~. , data=train_rf, method = "pls", metric=metric,
          trControl = Control, tuneLength=50)
pls
plot(pls)
confusionMatrix(pls)

# Bagging tree
set.seed(1)
metric <- "ROC"
Control=trainControl(method= "cv",number=10, classProbs=TRUE,
                     summaryFunction=twoClassSummary)
bagg <- train( u~. , data=Training, method="treebag", trControl=Control,
               importance=TRUE, metric=metric)

bagg
confusionMatrix(bagg)

#Gradient Boosting
set.seed(1234)
metric <- "ROC"
Control=trainControl(method= "cv",number=10, classProbs=TRUE,
                     summaryFunction=twoClassSummary)

gbm <- train(u~. , data=Training, method = "gbm", metric = metric,
    tuneLength = 5 , trControl = Control)
  
gbm
confusionMatrix(gbm)

#AdaBoost
set.seed(1)
metric <- "ROC"
Control=trainControl(method= "cv",number=10, classProbs=TRUE,
                     summaryFunction=twoClassSummary)
ada <- train( u~. , data=Training, method="adaboost", trControl=Control,
              importance=TRUE, metric=metric)

ada
confusionMatrix(ada)

#Reti neurali con il dataset ricavato con l'albero
#preProcess = "range" per normalizzare le covariate tra 0 e 1
set.seed(7)
metric <- "ROC"
ctrl = trainControl(method="cv", number=10, search = "grid", classProbs=TRUE, 
                    summaryFunction = twoClassSummary)
nnet_tree<- train(u~. , data=Train_modsel,
                  method = "nnet",
                  preProcess = "range", 
                  metric=metric, trControl=ctrl,
                  trace = TRUE, 
                  maxit = 300)
print(nnet_tree) 
nnet_tree
plot(nnet_tree)
getTrainPerf(nnet_tree) 

confusionMatrix(nnet_tree)

#Reti neurali con il dataset ottenuto con il random forest
#preProcess = "range" per normalizzare le covariate tra 0 e 1
library(caret)
library(nnet)
set.seed(7)
metric <- "ROC"
ctrl = trainControl(method="cv", number=10, search = "grid", classProbs=TRUE, 
                    summaryFunction = twoClassSummary)
nnet_rf<- train(u~. , data=train_rf,
                method = "nnet",
                preProcess = "range", 
                metric=metric, trControl=ctrl,
                trace = TRUE,
                maxit = 300)
print(nnet_rf) 
nnet_rf
plot(nnet_rf)
getTrainPerf(nnet_rf) 

confusionMatrix(nnet_rf)

##STEP 2: confronto modelli
results <- resamples(list(glm=glm, lasso=lasso, naivebayes_tree=naivebayes_tree, naivebayes_rf=naivebayes_rf,rf=rf, pls=pls, tree=tree, gbm=gbm, bagg=bagg, nnet_tree=nnet_tree, nnet_rf=nnet_rf, ada=ada))
summary(results)
bwplot(results)

#Valori previsti sul validation
head(Validation)
valid_rf$glm=predict(glm,valid_rf, "prob")[,2]
Validation$lasso=predict(lasso,Validation, "prob")[,2]
Valid_modsel$naivebayes_tree=predict(naivebayes_tree,Valid_modsel, "prob")[,2]
valid_rf$naivebayes_rf=predict(naivebayes_rf,valid_rf, "prob")[,2]
valid_rf$pls=predict(pls,valid_rf, "prob")[,2]
Validation$rf=predict(rf,Validation, "prob")[,2]
Validation$tree=predict(tree,Validation, "prob")[,2]
Validation$gbm=predict(gbm,Validation, "prob")[,2]
Validation$bagg=predict(bagg,Validation, "prob")[,2]
Valid_modsel$nnet_tree=predict(nnet_tree,Valid_modsel, "prob")[,2]
valid_rf$nnet_rf=predict(nnet_rf,valid_rf, "prob")[,2]
Validation$ada=predict(ada,Validation, "prob")[,2]
p_prova = predict(ada      , Validation, "prob")
head(p_prova)

#Curve ROC
roc.glm=roc(u ~ glm, data = valid_rf)
roc.lasso=roc(u ~ lasso, data = Validation)
roc.naivebayes_tree=roc(u ~ naivebayes_tree, data = Valid_modsel)
roc.naivebayes_rf=roc(u ~ naivebayes_rf, data = valid_rf)
roc.pls=roc(u ~ pls, data = valid_rf)
roc.rf=roc(u ~ rf, data = Validation)
roc.tree=roc(u ~ tree, data = Validation)
roc.bagg=roc(u ~ bagg, data = Validation)
roc.gbm=roc(u ~ gbm, data = Validation)
roc.nnet_tree=roc(u ~ nnet_tree, data = Valid_modsel)
roc.nnet_rf=roc(u ~ nnet_rf, data = valid_rf)
roc.ada=roc(u ~ ada, data = Validation)
roc.glm
roc.lasso
roc.naivebayes_tree
roc.naivebayes_rf
roc.pls
roc.tree
roc.rf
roc.bagg
roc.gbm
roc.nnet_tree
roc.nnet_rf
roc.ada

#Grafici curve ROC
plot(roc.glm)
plot(roc.lasso,add=T,col="blue")
plot(roc.naivebayes_tree,add=T,col="yellow")
plot(roc.pls,add=T,col="green")
plot(roc.tree,add=T,col="pink")
plot(roc.rf,add=T,col="violet")
plot(roc.naivebayes_rf,add=T,col="red")
plot(roc.bagg,add=T,col="orange")
plot(roc.gbm,add=T,col="black")
plot(roc.nnet_tree,add=T,col="purple")
plot(roc.nnet_rf,add=T,col="deepskyblue")
plot(roc.ada,add=T,col="brown")
legend( "bottomright", c("glm", "lasso", "naivebayes_tree", "pls", "tree", "rf", "naivebayes_rf", "bagg", "gbm", "nnet_tree", "nnet_rf", "ada"), lty=1, bty='n', cex=.65,
        col=c("grey", "blue", "yellow", "green", "pink", "violet", "red", "orange", "black", "purple", "deepskyblue", "brown") )

#curve liftchart
gain_lift(data = Validation, score = 'ada', target = 'u')
gain_lift(data = Validation, score = 'rf', target = 'u')
gain_lift(data = Validation, score = 'gbm', target = 'u')

#STEP 3: modello vincente
gbm
confusionMatrix(gbm)
pred_y<-ifelse(Validation$gbm>0.5, "yes","no")
pred_y=as.factor(pred_y)
head(pred_y)
actual=Validation$u
actual=as.factor(actual)
y=Validation$u
y=ifelse(y=="yes",1,0)
table(y)
predProbT=Validation$gbm

predRoc <- prediction(predProbT,y)
class(predRoc)

roc.perf = performance(predRoc, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1) 

#studio della soglia
sens.perf = performance(predRoc, measure = "sens")
plot(sens.perf)

sens=as.data.frame(sens.perf@y.values)
colnames(sens)="spec"
head(sens)

performance(predRoc,measure="spec")

Validation$pR = predict(gbm      , Validation, "prob")[,2]
head(Validation$pR)

gain_lift(data = Validation, score = 'pR', target = 'u')

plot(sens.perf)

##STEP 4: score del modello vincente
Score$prob = predict(gbm, Score, "prob")
head(Score$prob)
probY=Score$prob[,2]
Score$pred_y=ifelse(probY>0.2, "Y","N")
head(Score)
table(Score$pred_y)
table(Score$pred_y)/nrow(Score)
