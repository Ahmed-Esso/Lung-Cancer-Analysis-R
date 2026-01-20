# 🫁 Lung Cancer Level Prediction Using R

A comprehensive machine learning project for predicting lung cancer severity levels (Low, Medium, High) using various classification algorithms in R.

## 📋 Project Overview

This project analyzes patient data to predict lung cancer severity levels using multiple machine learning models. It includes extensive exploratory data analysis (EDA), statistical testing, and model comparison to identify the best-performing algorithm.

## 📊 Dataset

The dataset contains patient information with various health indicators including:
- **Demographics**: Age, Gender
- **Risk Factors**: Smoking, Genetic Risk, Air Pollution exposure
- **Symptoms**: Chest Pain, Coughing of Blood, Shortness of Breath, Wheezing
- **Target Variable**: Cancer Level (Low, Medium, High)

## 🛠️ Technologies & Libraries

```r
library(tidyverse)      # Data manipulation & visualization
library(caret)          # Machine learning framework
library(corrplot)       # Correlation visualization
library(rpart)          # Decision trees
library(rpart.plot)     # Decision tree visualization
library(randomForest)   # Random Forest algorithm
library(xgboost)        # XGBoost algorithm
library(e1071)          # Naive Bayes & SVM
library(nnet)           # Neural Networks
library(class)          # KNN algorithm
```

## 📈 Analysis Pipeline

### 1. Data Preprocessing
- Loading and cleaning patient data
- Handling missing values
- Feature encoding and scaling
- Train-test split (80/20)

### 2. Exploratory Data Analysis
- Cancer level distribution (Donut chart)
- Age distribution by cancer level
- Smoking vs cancer level relationship
- Correlation heatmap of numeric features
- Symptom severity analysis
- Gender and genetic risk analysis

### 3. Statistical Testing
- **ANOVA**: Age differences across cancer levels
- **Chi-square tests**: Smoking and genetic risk associations
- **Kruskal-Wallis test**: Shortness of breath severity

### 4. Machine Learning Models
| Model | Description |
|-------|-------------|
| Logistic Regression | Multinomial classification |
| Decision Tree | CART algorithm with visualization |
| K-Nearest Neighbors | k=5 with scaled features |
| Random Forest | 300 trees with variable importance |
| XGBoost | Gradient boosting with tuned hyperparameters |
| Naive Bayes | Probabilistic classifier |
| Neural Network | Single hidden layer (5 neurons) |
| Random Forest CV | 5-fold cross-validated RF |

## 📁 Project Structure

```
Lung-Cancer-Analysis-R/
├── README.md
├── lung-cancer.R              # Main analysis script
├── cancer-patient-data-sets.csv   # Dataset
└── Lung-Cancer-Report.pptx    # Presentation report
```

## 🚀 Getting Started

### Prerequisites
Make sure you have R (≥ 4.0) and RStudio installed.

### Installation

1. Clone this repository:
```bash
git clone https://github.com/Ahmed-Esso/Lung-Cancer-Analysis-R.git
cd Lung-Cancer-Analysis-R
```

2. Install required packages:
```r
install.packages(c("tidyverse", "caret", "corrplot", "rpart", 
                   "rpart.plot", "randomForest", "xgboost", 
                   "e1071", "nnet", "class"))
```

3. Run the analysis:
```r
source("lung-cancer.R")
```

## 📊 Key Visualizations

- 🍩 Cancer level distribution donut chart
- 📊 Age distribution histograms by level
- 🔥 Correlation heatmap
- 📦 Boxplots for continuous variables
- 🌲 Decision tree visualization
- 📈 Model accuracy comparison bar chart
- 🎯 Variable importance plots

## 🎯 Results

The project compares all models and displays:
- Accuracy metrics for each model
- Confusion matrices
- Variable importance rankings
- Visual comparison of model performance

## 👤 Author

**Ahmed Essam**
- GitHub: [@Ahmed-Esso](https://github.com/Ahmed-Esso)
- Website: [ahmed-essam.framer.website](https://ahmed-essam.framer.website/)
- Twitter: [@Ahmed__Esso](https://twitter.com/Ahmed__Esso)

## 📄 License

This project is open source and available under the [MIT License](LICENSE).

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the issues page.

---
⭐ If you found this project helpful, please give it a star!
