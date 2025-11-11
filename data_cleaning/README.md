
# Data Cleaning

* `methylation_pkl2feather.ipynb` changes data format from pickle to feather so that they can be read into R.
* `geneexp_cleaning.R` apply normalization and log transformation to gene counts, then split into train and test set.
* `methylation_data_split.R` split each methylation study's data into training and test set. Only studies with at least 80 samples were retained.
