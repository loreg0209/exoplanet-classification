# Exoplanet Transit Timing Variation Classification

Machine learning project for predicting **Transit Timing Variations (TTV)** in exoplanets using supervised classification models in **R**.

The objective is to identify whether an exoplanet exhibits variations in its transit timing (`ttv_flag`) based on its physical, orbital and stellar characteristics. Detecting these variations can provide indirect evidence of additional bodies in a planetary system that may not be directly observable.

---

## Objective

The target variable is:

```text
ttv_flag
```

where:

- `1` indicates the presence of Transit Timing Variations
- `0` indicates no observed variation

The goal is to identify the most effective classifier for recognizing exoplanets with TTV while limiting false negatives.

The original dataset contains **5,514 exoplanets** and a wide set of orbital, stellar and discovery-related attributes.

---

## Preprocessing

The preprocessing pipeline includes:

- removal of identifier variables
- removal of redundant representations
- analysis of zero-variance features
- analysis of missing values
- exclusion of variables with more than 20% missingness
- optimal grouping of high-cardinality categorical variables
- Random Forest imputation
- collinearity analysis
- removal of highly correlated predictors
- class balancing through undersampling
- stratified Train / Validation / Score partitioning

The final split is:

- **60% Training**
- **30% Validation**
- **10% Score**

---

## Missing Data

Variables with a high percentage of missing observations are excluded from the analysis.

For the remaining predictors, missing values are estimated using **Random Forest imputation** through `missForest`.

This approach allows the analysis to preserve observations while handling both numeric and categorical predictors.

---

## Collinearity

Strong correlations were identified among several pairs of variables, including:

- planetary radius expressed in Earth and Jupiter units
- planetary mass expressed in Earth and Jupiter units
- different stellar magnitude measures

Redundant variables were removed before model training to reduce unnecessary collinearity.

---

## Class Imbalance

The target variable is strongly imbalanced, with the majority of observations belonging to the non-TTV class.

To improve the ability of the classifiers to recognize positive cases, the majority class is reduced using **stratified undersampling** before creating the final Training, Validation and Score datasets.

---

## Competitive Models

The project compares a broad set of supervised learning algorithms:

1. Decision Tree
2. Random Forest
3. Logistic Regression
4. Lasso Logistic Regression
5. Naive Bayes with Tree-based feature selection
6. Naive Bayes with Random-Forest-based feature selection
7. Partial Least Squares (PLS)
8. Bagging
9. Gradient Boosting
10. AdaBoost
11. Neural Network with Tree-based feature selection
12. Neural Network with Random-Forest-based feature selection

Models are tuned and compared primarily using **ROC-based evaluation**.

---

## Model Validation

Model performance is evaluated on a dedicated validation set.

ROC curves are used to compare the discriminative ability of all candidate models.

The strongest classifiers are further compared using **Lift curves**, focusing on the first two population deciles.

The results show that:

- Gradient Boosting captures **90.6%** of the target population within the first two deciles
- AdaBoost captures **89.6%**
- Random Forest captures **88.5%**

Based on this comparison, **Gradient Boosting** is selected as the final model.

---

## Threshold Selection

After selecting the final model, the classification threshold is chosen by focusing on **Sensitivity**.

This is particularly important because a false negative corresponds to an exoplanet with TTV being incorrectly classified as a non-TTV case.

The selected probability threshold is:

```text
0.20
```

This operating point prioritizes the correct identification of positive TTV cases.

---

## Final Model

The final classifier is:

```text
Gradient Boosting
```

with a decision threshold of:

```text
0.20
```

The model is then applied to the independent Score dataset to obtain the final classifications.

---

## Key Takeaways

The project demonstrates a complete supervised machine learning workflow including:

- data preprocessing
- missing-value management
- feature screening
- class-imbalance handling
- model comparison
- cross-validation
- ROC analysis
- Lift analysis
- threshold optimization
- final scoring

A key result is that model selection based only on overall accuracy would not be sufficient for this task. Evaluating ranking performance and sensitivity provides a more useful framework for identifying the minority TTV class.

---

## Repository Structure

```text
exoplanet-ttv-classification/
├── code/
│   └── exoplanet_ttv_classification.R
├── data/
│   └── esopianeti.csv
├── docs/
│   ├── Exoplanet_TTV_Machine_Learning_Report.pdf
└── README.md
```

