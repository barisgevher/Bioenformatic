library(GEOquery)
library(genefilter)
library(caret)
library(ggplot2)
library(rpart)

okunan1 <- getGEO("GSE6008", GSEMatrix = TRUE)
eset1 <- okunan1[[1]]
kopya <- exprs(eset1)
colnames(pData(eset1))
dim(eset1)


BiocManager::install("hgu133plus2.db",force = TRUE)
annotation(eset1) <- "hgu133plus2.db"

durum <- factor(pData(eset1)$characteristics_ch1)

filtrelenmis <- varFilter(eset1, var.cutoff = 0.9)
sonveri <- data.frame(t(exprs(filtrelenmis)))
sonveri$durum <- durum
dim(sonveri)

set.seed(123)
train_indices <- createDataPartition(durum, p = 0.7, list = FALSE)
train_data <- sonveri[train_indices, ]
test_data <- sonveri[-train_indices, ]


model_rf <- train(durum ~ ., data = train_data, method = "rf")
predictions_rf <- predict(model_rf, newdata = test_data)
conf_matrix_rf <- confusionMatrix(predictions_rf, test_data$durum)
accuracy_rf <- conf_matrix_rf$overall["Accuracy"]
print(paste("Random Forest Doğruluk:", accuracy_rf))
feature_importance_rf <- varImp(model_rf)
print(feature_importance_rf)
print(conf_matrix_rf)


tree_model <- rpart(durum ~ ., data = sonveri, method = "class")
predictions_tree <- predict(tree_model, newdata = test_data, type = "class")
conf_matrix_tree <- confusionMatrix(predictions_tree, test_data$durum)
accuracy_tree <- conf_matrix_tree$overall["Accuracy"]
print(paste("Karar Ağacı Doğruluk:", accuracy_tree))
importance_tree <- tree_model$variable.importance
top_features <- importance_tree[order(-importance_tree, decreasing = TRUE)][1:20]
print(top_features)
print(conf_matrix_tree)


svm_model <- train(durum ~ ., data = train_data, method = "svmLinear")
predictions_svm <- predict(svm_model, newdata = test_data)
conf_matrix_svm <- confusionMatrix(predictions_svm, test_data$durum)
accuracy_svm <- conf_matrix_svm$overall["Accuracy"]
print(paste("SVM Doğruluk:", accuracy_svm))
feature_importance_svm <- varImp(svm_model)
print(feature_importance_svm)
print(conf_matrix_svm)
