# AgePred
Experiments that compare different machine learning methods' age prediction capability.

* `data_cleaning` folder contains data cleaning codes that split each study/tissue into training and test set. For gene expressions, I used DeSeq2 for normalization first.
* `existing_methods` folder contains training codes for competitors such as elastic net, random forest and lightGBM.
* `new_methods` folder contains training codes for deconvolution.
