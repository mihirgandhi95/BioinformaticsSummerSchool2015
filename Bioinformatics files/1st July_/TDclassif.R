rm(list=ls())

library(spikeslab)
library(glmnet)
library(ROCR)

#################################################################
# Data set
#################################################################
# Golub data
data(leukemia)
n = nrow(leukemia); p = ncol(leukemia)-1
X = as.matrix(leukemia[, -1]); Y = as.vector(leukemia[, 1])
Gene = names(leukemia)[-1]

#################################################################
# Train / test sets
#################################################################
set.seed(1)
Train = sort(sample(1:n, round(2*n/3))); ntrain = length(Train)
Test = sort(setdiff(1:n, Train)); ntest = length(Test)
Xtrain = X[Train, ]; Ytrain = Y[Train]
Xtest = X[Test, ]; Ytest = Y[Test]

# Dimensions 
cat(n, p, ntrain, ntest, '\n')

#################################################################
# Elastic-net logistic regression
#################################################################
# alpha = 1 -> Lasso; alpha = 0 -> Ridge
alpha = 1

# Solution path
EN = glmnet(Xtrain, Ytrain, family="binomial", alpha=1)
plot(EN)
# Cross-validation
ENcv = cv.glmnet(Xtrain, Ytrain, family="binomial", alpha=1)
plot(ENcv)
EN = glmnet(Xtrain, Ytrain, family="binomial", alpha=1, lambda=ENcv$lambda.1se)
EN = glmnet(Xtrain, Ytrain, family="binomial", alpha=1, lambda=ENcv$lambda.min)
# Selected genes
Gene[which(EN$beta!=0)]
# Prediction
plot(predict(EN, Xtrain), Ytrain)
points(predict(EN, Xtest), Ytest, col=2)
# ROC curve
plot(performance(prediction(predict(EN, Xtrain), Ytrain), "tpr", "fpr"))
plot(performance(prediction(predict(EN, Xtest), Ytest), "tpr", "fpr"))

#################################################################
# Stability selection
#################################################################
B = 1e3
Sel.min = rep(0, p); Sel.1se = Sel.min
for (b in 1:B){
   if(b %% round(sqrt(B))==0){cat(b, '')}
   # Training sample
   Train = sort(sample(1:n, round(2*n/3))); ntrain = length(Train)
   Xtrain = X[Train, ]; Ytrain = Y[Train]
   ## Cross validation
   #ENcv = cv.glmnet(Xtrain, Ytrain, family="binomial", alpha=1)
   # Selected genes
   EN = glmnet(Xtrain, Ytrain, family="binomial", alpha=1, lambda=ENcv$lambda.1se)
   Sel.1se[which(EN$beta!=0)] = Sel.1se[which(EN$beta!=0)]+1
   EN = glmnet(Xtrain, Ytrain, family="binomial", alpha=1, lambda=ENcv$lambda.min)
   Sel.min[which(EN$beta!=0)] = Sel.min[which(EN$beta!=0)]+1
}
plot(sort(Sel.1se))
Gene[which(Sel.1se > B/2)]
Sel.1se[which(Sel.1se > B/2)]
plot(sort(Sel.min))
Gene[which(Sel.min > B/2)]
Sel.min[which(Sel.min > B/2)]

#################################################################
# Gene filtering
#################################################################
Pval = rep(0, p)
invisible(sapply(1:n, function(j){
#    Pval[j] <<- t.test(X[, j], Y)$p.value
   Pval[j] <<- wilcox.test(X[, j], Y)$p.value
}))
hist(Pval, breaks=sqrt(p))
PvalAdj = p.adjust(Pval, method="bonferroni")
length(Gene[which(PvalAdj < .05)])
