---
Package Name: "MWAS R Package v0.9.3"
Developers: "Hu Huang,Emmanuel Montassier, Pajau Vangay, Gabe Al Ghalith, Dan Knights"
Date: "03-01-2015"
---

****


**MWAS (microbiome-wide association study)** package is a R-based toolbox for microbiome study, developped by the members of the Knights Lab at the University of Minnesota, Twin Cities. It provides three main functional modules: learning a predictve model, predicting an unknown microbiome data, and visualization of different results. The latest update is version 0.9.3 (03-2015). 

MWAS is developed in R, however, it also provides a Unix command-line interface as a simplified application for those who are not familliar with R.  

****

#### Quick Index
  
#### 1. [Installing MWAS](http://rpubs.com/hwangtiger/install_mwas)  

Click the above link (section title) for detailed information.  

* Use the following command to set `MWAS_DIR` in the Terminal (or an equivalent command window; `/MWAS_directory` should be your actual directory):   
`echo "export MWAS_DIR=$HOME/MWAS_directory" >> ~/.bash_profile`  

* You might need to install dependencies seperately, if it cannot install or load the required packages. Most of the dependencies would be installed when running the corresponding function commands, except one pacakge `optparse`. Follow the steps below to install this package:  
   + Open R Console in Terminal (or use RStudio)  
   + Install the pacakge: `install.packages("optparse")`  
   + You should be able to use the MWAS functions now.   
   
*(Detailed testing information is available [here](http://rpubs.com/hwangtiger/install_mwas).)*

***

#### 2. [MWAS "*learn*" Module](http://rpubs.com/hwangtiger/mwas_learn)  

* **Command-line version (in Terminal)**     
`Rscript MWAS_DIR/bin/mwas_analysis.R -w learn -M SVM -C linear -i MWAS_DIR/data/taxa/GG_100nt_even10k-adults_L7.biom -m MWAS_DIR/data/gg-map-adults.txt -o example/svm_output -c COUNTRY -f -v FDR -s 0.05`

`-w`: learn mode  
`-M`: classifier type  
`-C`: kernel type for SVM  
`-i`: input file  
`-m`: mapfile  
`-c`: category name  
`-o`: output directory  
`-f`: proceed feature selection
`-v`: feature selection method: `fdr` or `rf`  
`-s`: threshold for feature selection (determines the number of features)

* **R version (in R Console)**

If you are familiar with R, you could manipulate your data in a more flexible way. Here is the same example as shown in the command-line version.  

`opts <- list()`  
`opts$mode <- "learn"`  
`opts$method <- "SVM"`  
`opts$input_fp <- "data/taxa/GG_100nt_even10k-adults_L7.biom"`  
`opts$map_fp <- "data/gg-map-adults.txt"`  
`opts$category <- "COUNTRY"`  
`opts$outdir <- "example/svm_learn"`  
`opts$nfolds <- 5`  
`opts$method_param <- "linear"`  
`opts$ftMethod <- "FDR"`  
`opts$is_feat <- TRUE`  
`opts$feat_param <- 0.05`  

`train_params <- import.train.params(opts)`  
`best_model <- train.mwas(train_params)`  


***

#### 3. [MWAS "*predict*" Module](http://rpubs.com/hwangtiger/mwas_predict)

***

#### 4. [MWAS Visualization Module](http://rpubs.com/hwangtiger/mwas_visualization)  

***

#### 5. [MWAS feature "*statistics*" Module](http://rpubs.com/hwangtiger/MWAS_feat_stats)

***

#### 6. Example 1: Learning a predictive model   




***
#### 7. Example 2: Prediction from an unknown dataset  

***

#### 8. Example 3: Taxon Statistical Analysis and Visualization  

***

#### 9. Common Errors and Solutions

***

### References  
Breiman, L. (2001). **Random forests**. Machine learning. 45(1), 5-32.  
Leo Breiman and Adele Cutler. (2003) **Random Forest - Classification Description**. Retrieved on November 1, 2014 from [http://www.math.usu.edu/~adele/forests/cc_home.htm](http://www.math.usu.edu/~adele/forests/cc_home.htm)  
Hastie, T., Tibshirani, R., Friedman, J., Hastie, T., Friedman, J., & Tibshirani, R. (2009). **The Elements of Statistical Learning** (Vol. 2, No. 1). New York: springer.  
Adbi, H., & Williams, L. J. (2010). **Jackknife**. In: Neil Salkind (Ed.), _Encyclopedia of Research Design_. Thousand Oaks, CA: Sage.   
Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J. C., & Müller, M. (2011). **pROC: an open-source package for R and S+ to analyze and compare ROC curves**. BMC bioinformatics, 12(1), 77.  
Chang, C. C., & Lin, C. J. (2011). **LIBSVM: a library for support vector machines**. ACM Transactions on Intelligent Systems and Technology (TIST), 2(3), 27.  
Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). **kernlab-an S4 package for kernel methods in R**. Journal of Statistical Software, 11(9)  
Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). **Regularization paths for generalized linear models via coordinate descent**. Journal of Statistical Software, 33(1), 1-22.  
Ben-Hur, A., & Weston, J. (2010). **A user’s guide to support vector machines**. In: O. Carugo, F. Eisenhaber (eds.), _Data mining techniques for the life sciences_ (pp. 223-239). Humana Press.  
Cawley, G. C., & Talbot, N. L. (2010). **On over-fitting in model selection and subsequent selection bias in performance evaluation**. The Journal of Machine Learning Research, 11, 2079-2107.  
Storey, J. D. (2002). **A direct approach to false discovery rates**. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(3), 479-498.  
Noble, W. S. (2009). **How does multiple testing correction work?** Nature biotechnology, 27(12), 1135-1137.  
Storey, J.D. (2010). **False discovery rate**. Retrieved on Feb. 1, 2015, from [http://www.genomine.org/papers/Storey_FDR_2010.pdf](http://www.genomine.org/papers/Storey_FDR_2010.pdf)  
Hu Huang, Emmanuel Montassier, Pajau Vangay, Gabe Al Ghalith, Dan Knights. "**Robust statistical models for microbiome phenotype prediction with the MWAS package**" (in preparation)

