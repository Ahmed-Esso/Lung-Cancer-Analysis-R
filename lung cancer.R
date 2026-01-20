# 1 - Imports and data loading
library(tidyverse); library(caret); library(corrplot); library(rpart); library(rpart.plot)
library(randomForest); library(xgboost); library(e1071); library(nnet); library(class)

file_path <- "cancer patient data sets.csv"

df <- read.csv(file_path, stringsAsFactors = FALSE)
str(df); summary(df); head(df)

id_cols <- c("index", "Patient.Id", "Patient Id")
id_cols_present <- intersect(id_cols, names(df))
df <- df[, !(names(df) %in% id_cols_present)]

if ("Level" %in% names(df)) df$Level <- as.factor(df$Level)
colSums(is.na(df))

# 2 - Exploratory data analysis and visuals

# New color palette
level_palette <- c(Low = "#003049", Medium = "#f77f00", High = "#d62828")

# Global theme
theme_set(theme_minimal(base_size = 14) +
            theme(
              plot.background   = element_rect(fill = "white", color = NA),
              panel.background  = element_rect(fill = "white", color = NA),
              panel.grid.major  = element_line(color = "#e0e0e0", linewidth = 0.3),
              panel.grid.minor  = element_blank(),
              axis.text         = element_text(color = "#333333"),
              axis.title        = element_text(color = "#263238", face = "bold"),
              plot.title        = element_text(color = "#1b1f32", face = "bold", hjust = 0.5),
              legend.title      = element_text(face = "bold", color = "#263238"),
              legend.text       = element_text(color = "#333333"),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key        = element_rect(fill = "white", color = NA)
            )
)

# 2.1 Level distribution - donut chart
if ("Level" %in% names(df)) {
  level_counts <- df %>% count(Level) %>%
    mutate(prop = n / sum(n), ypos = cumsum(prop) - 0.5 * prop)
  
  ggplot(level_counts, aes(x = 2, y = prop, fill = Level)) +
    geom_col(color = "white") + coord_polar(theta = "y") + xlim(0.5, 2.5) +
    scale_fill_manual(values = level_palette) +
    geom_text(aes(y = ypos, label = paste0(Level, " (", n, ")")), color = "white", size = 4) +
    theme_void() + labs(title = "Lung Cancer Level Distribution")
}

# 2.2 Age distribution by level (side-by-side histograms)
if (all(c("Age", "Level") %in% names(df))) {
  ggplot(df, aes(x = Age, fill = Level)) +
    geom_histogram(binwidth = 5, color = "black", alpha = 0.85) +
    scale_fill_manual(values = level_palette, guide = "none") +
    labs(title = "Age Distribution by Cancer Level", x = "Age", y = "Count") +
    facet_wrap(~ Level, nrow = 1)
}

# 2.3 Smoking vs level
if (all(c("Smoking", "Level") %in% names(df))) {
  chi_smoking <- chisq.test(table(df$Smoking, df$Level)); p_smoking <- chi_smoking$p.value
  smoking_plot_data <- df %>%
    group_by(Smoking, Level) %>% summarise(n = n(), .groups = "drop") %>%
    group_by(Smoking) %>% mutate(prop = n / sum(n), label = paste0(round(prop * 100), "%"))
  
  ggplot(smoking_plot_data, aes(x = as.factor(Smoking), y = prop, fill = Level)) +
    geom_col(position = "fill", color = "black") +
    geom_text(color = "white", aes(label = label), position = position_fill(vjust = 0.5), size = 3) +
    scale_fill_manual(values = level_palette) +
    labs(title = paste0("Smoking Score vs Cancer Level (Chi-square p = ", signif(p_smoking, 3), ")"),
         x = "Smoking score", y = "Proportion")
}

# 2.4 Correlation heatmap
num_vars <- df %>% select(where(is.numeric))
if (ncol(num_vars) > 1) {
  corr_mat <- cor(num_vars)
  calm_pal <- colorRampPalette(c("#003049", "#f77f00", "#d62828"))
  corrplot(corr_mat, method = "color", type = "full", diag = TRUE,
           tl.cex = 0.7, tl.col = "#263238", tl.srt = 45,
           col = calm_pal(200), addCoef.col = "black",
           number.digits = 2, number.cex = 0.6, bg = "white")
}

# 2.5 Age boxplot by level
if (all(c("Age", "Level") %in% names(df))) {
  anova_age <- aov(Age ~ Level, data = df); p_age <- summary(anova_age)[[1]][["Pr(>F)"]][1]
  ggplot(df, aes(x = Level, y = Age, fill = Level)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.6, color = "black") +
    scale_fill_manual(values = level_palette) +
    labs(title = paste0("Age Distribution by Cancer Level (ANOVA p = ", signif(p_age, 3), ")"),
         x = "Cancer level", y = "Age")
}

# 2.6 Key symptoms vs level
symptom_vars <- c("Chest.Pain", "Coughing.of.Blood", "Shortness.of.Breath", "Wheezing")
if (all(symptom_vars %in% names(df))) {
  symptoms_long <- df %>% select(Level, all_of(symptom_vars)) %>%
    pivot_longer(cols = all_of(symptom_vars), names_to = "Symptom", values_to = "Score")
  
  ggplot(symptoms_long, aes(x = as.factor(Score), fill = Level)) +
    geom_bar(position = "fill", color = "black") +
    scale_fill_manual(values = level_palette) +
    facet_wrap(~ Symptom, ncol = 2) +
    labs(title = "Symptom Severity Scores by Cancer Level", x = "Severity score", y = "Proportion")
}

# 2.7 Gender vs level
if (all(c("Gender", "Level") %in% names(df))) {
  df$GenderLabel <- factor(df$Gender, levels = c(1, 2), labels = c("Male", "Female"))
  gender_plot_data <- df %>% group_by(GenderLabel, Level) %>% summarise(n = n(), .groups = "drop")
  ggplot(gender_plot_data, aes(x = GenderLabel, y = n, fill = Level)) +
    geom_col(position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = n), color = "black",
              position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = level_palette) +
    labs(title = "Gender Distribution by Cancer Level", x = "Gender", y = "Count")
}

# 2.8 Genetic risk vs level
if (all(c("Genetic.Risk", "Level") %in% names(df))) {
  ggplot(df, aes(x = as.factor(Genetic.Risk), fill = Level)) +
    geom_bar(position = "fill", color = "black") +
    scale_fill_manual(values = level_palette) +
    labs(title = "Genetic Risk Score vs Cancer Level", x = "Genetic risk score", y = "Proportion")
}

# 2.9 Statistical tests (numerical output)
if (all(c("Age", "Level") %in% names(df))) {
  anova_age <- aov(Age ~ Level, data = df); print(summary(anova_age))
}
if (all(c("Smoking", "Level") %in% names(df))) {
  chi_smoking_num <- chisq.test(table(df$Smoking, df$Level)); print(chi_smoking_num)
}
if (all(c("Genetic.Risk", "Level") %in% names(df))) {
  chi_genetic <- chisq.test(table(df$Genetic.Risk, df$Level)); print(chi_genetic)
}
if (all(c("Shortness.of.Breath", "Level") %in% names(df))) {
  kw_breath <- kruskal.test(Shortness.of.Breath ~ Level, data = df); print(kw_breath)
}

# 2.10 Shortness of breath vs level
if (all(c("Shortness.of.Breath", "Level") %in% names(df))) {
  kw_breath <- kruskal.test(Shortness.of.Breath ~ Level, data = df); p_kw <- kw_breath$p.value
  ggplot(df, aes(x = Level, y = Shortness.of.Breath, fill = Level)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.6, color = "black") +
    scale_fill_manual(values = level_palette) +
    labs(title = paste0("Shortness of Breath by Cancer Level (Kruskal p = ", signif(p_kw, 3), ")"),
         x = "Cancer level", y = "Shortness of breath score")
}

# 3 - Data preprocessing
if ("Level" %in% names(df)) df$Level <- factor(df$Level) else stop("Target column 'Level' not found in the data")

set.seed(123)
train_index <- createDataPartition(df$Level, p = 0.8, list = FALSE)
train_data <- df[train_index, ]; test_data <- df[-train_index, ]

train_data <- train_data %>% select(-any_of("GenderLabel"))
test_data  <- test_data  %>% select(-any_of("GenderLabel"))

x_train <- train_data %>% select(-Level); y_train <- train_data$Level
x_test  <- test_data %>% select(-Level);  y_test  <- test_data$Level

preproc <- preProcess(x_train, method = c("center", "scale"))
x_train_scaled <- predict(preproc, x_train)
x_test_scaled  <- predict(preproc, x_test)

label_levels <- levels(df$Level)
y_train_xgb <- as.numeric(y_train) - 1
y_test_xgb  <- as.numeric(y_test) - 1

dtrain <- xgb.DMatrix(data = as.matrix(x_train), label = y_train_xgb)
dtest  <- xgb.DMatrix(data = as.matrix(x_test), label = y_test_xgb)

# 4 - Models
eval_model <- function(true, pred, model_name) {
  true <- factor(true); pred <- factor(pred, levels = levels(true))
  cm <- confusionMatrix(pred, true)
  list(model = model_name, accuracy = cm$overall["Accuracy"], cm = cm)
}

results_list <- list()

# 4.1 Logistic regression
logit_model <- multinom(Level ~ ., data = train_data, trace = FALSE)
logit_pred  <- predict(logit_model, newdata = test_data)
res_logit   <- eval_model(y_test, logit_pred, "Logistic Regression")
results_list[["Logistic Regression"]] <- res_logit

# 4.2 Decision tree
par(bg = "white")
tree_model <- rpart(Level ~ ., data = train_data, method = "class")
rpart.plot(tree_model, main = "Decision tree for lung cancer level",
           box.palette = list("#d62828", "#003049", "#f77f00"),
           shadow.col = "#d0d0d0", col = "white", nn = TRUE)

tree_pred <- predict(tree_model, newdata = test_data, type = "class")
res_tree  <- eval_model(y_test, tree_pred, "Decision Tree")
results_list[["Decision Tree"]] <- res_tree

# 4.3 K Nearest Neighbors
k_value <- 5
knn_pred <- knn(train = as.matrix(x_train_scaled), test = as.matrix(x_test_scaled),
                cl = y_train, k = k_value)
res_knn <- eval_model(y_test, knn_pred, paste0("KNN (k=", k_value, ")"))
results_list[["KNN"]] <- res_knn

# 4.4 Random Forest
rf_model <- randomForest(x = x_train, y = y_train, ntree = 300,
                         mtry = floor(sqrt(ncol(x_train))), importance = TRUE)
rf_pred <- predict(rf_model, newdata = x_test)
res_rf  <- eval_model(y_test, rf_pred, "Random Forest")
results_list[["Random Forest"]] <- res_rf

par(bg = "white")
varImpPlot(rf_model, main = "Random forest variable importance", col = "black")

# 4.5 XGBoost
xgb_params <- list(objective = "multi:softmax", num_class = length(label_levels),
                   eval_metric = "merror", max_depth = 4, eta = 0.1,
                   subsample = 0.8, colsample_bytree = 0.8)

xgb_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 150,
                       watchlist = list(train = dtrain, test = dtest), verbose = 0)

xgb_pred_num <- predict(xgb_model, dtest)
xgb_pred     <- factor(label_levels[xgb_pred_num + 1], levels = label_levels)
res_xgb      <- eval_model(y_test, xgb_pred, "XGBoost")
results_list[["XGBoost"]] <- res_xgb

# 4.6 Naive Bayes
nb_model <- naiveBayes(Level ~ ., data = train_data)
nb_pred  <- predict(nb_model, newdata = test_data)
res_nb   <- eval_model(y_test, nb_pred, "Naive Bayes")
results_list[["Naive Bayes"]] <- res_nb

# 4.7 Neural Network
nn_train <- cbind(x_train_scaled, Level = y_train)
nn_model <- nnet(Level ~ ., data = nn_train, size = 5, maxit = 300, trace = FALSE)
nn_pred  <- predict(nn_model, newdata = x_test_scaled, type = "class")
res_nn   <- eval_model(y_test, nn_pred, "Neural Network")
results_list[["Neural Network"]] <- res_nn

# 5 - Trying all features with cross validated Random Forest
control <- trainControl(method = "cv", number = 5)
rf_cv <- train(Level ~ ., data = train_data, method = "rf", trControl = control)
rf_cv_pred <- predict(rf_cv, newdata = test_data)
res_rf_cv  <- eval_model(y_test, rf_cv_pred, "Random Forest CV")
results_list[["Random Forest CV"]] <- res_rf_cv

# 6 - Evaluation summary
model_names <- names(results_list)
accuracies  <- sapply(results_list, function(x) as.numeric(x$accuracy))

eval_table <- data.frame(Model = model_names, Accuracy = round(accuracies, 4))
print(eval_table[order(-eval_table$Accuracy), ])

# 6.1 Visual comparison of model accuracies
plot_eval <- eval_table %>%
  arrange(desc(Accuracy)) %>%
  ggplot(aes(x = reorder(Model, Accuracy), y = Accuracy, fill = Accuracy)) +
  geom_col(color = "black") +
  geom_text(aes(label = Accuracy), vjust = -0.3, size = 4) +
  scale_fill_gradient(low = "#003049", high = "#d62828") +
  labs(title = "Model Accuracy Comparison", x = "Model", y = "Accuracy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

print(plot_eval)
