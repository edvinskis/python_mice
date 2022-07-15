
### Can a `Python` package do what `mice` can?

Most real world datasets contain at least some missing values and complicate the analysis of data. The de-facto standard for imputation method of incomplete data in `R` is `mice`, that solves the missing data problem iteratively on a variable-by-variable basis. It is well established that mice can yield valid inferences under many missing data conditions. However, there is not such of an ubiquitous choice for handling missing data in `Python`.

This repository contains code for a model-based simulation study that is used to evaluate different `Python` imputation methods under different missingness mechanisms and proportions to whether they can produce valid inferences. The `Python` imputation methods `KNNImputer`, `IterativeImputer`, `miceforest` and `MIDASpy` are considered in this study.