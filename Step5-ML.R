#LASSO回归，输入只能是矩阵。行为样本，列为基因，分组只能是（0,1），多组可以转化MultinomialExample
library(glmnet)
# data("BinomialExample")
# x <- BinomialExample$x
# y <- BinomialExample$y
x <- input[,-1]
levels(group) <- c(0, 1)
y <- as.integer(as.character(Group))
fit <- glmnet(x, y, family = "binomial", nlambda = 100, alpha = 1)
plot(fit, xvar = "lambda", label = TRUE)
cvfit <- cv.glmnet(data.matrix(x), as.numeric(y), nfolds = 10)
plot(cvfit)
coef(cvfit, s = "lambda.min")
cvfit$lambda.min
cvfit$lambda.1se

pdf("lasso.pdf", height = 8, width = 10)
plot(cvfit)
dev.off()
coef <- coef(cvfit, s = "lambda.min")
factors <- as.matrix(coef)
factors[factors == 0] <- NA
factors <- na.omit(factors)
factors
write.csv(factors, "lasso.csv")

#计算诊断评分
predict_model<-glm(Group~venn_gene,data=input)
predict_score<-predict(predict_model)
predict_score<-data.frame(predict_score)
