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

#### 2. MWAS "*learn*" Module  

* **Command-line version (in Terminal)**     
`Rscript bin/mwas_analysis.R -w learn -M SVM -C linear -i data/taxa/GG_100nt_even10k-adults_L7.biom -m data/gg-map-adults.txt -o example/svm_output -c COUNTRY -f -v FDR -s 0.05`

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

***

#### 3. MWAS "*predict*" Module

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

Storey, J. D. (2002). **A direct approach to false discovery rates**. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(3), 479-498.  
Noble, W. S. (2009). **How does multiple testing correction work?** Nature biotechnology, 27(12), 1135-1137.  
Storey, J.D. (2010). **False discovery rate**. Retrieved on Feb. 1, 2015, from [http://www.genomine.org/papers/Storey_FDR_2010.pdf](http://www.genomine.org/papers/Storey_FDR_2010.pdf)  
Hu Huang, Emmanuel Montassier, Pajau Vangay, Gabe Al Ghalith, Dan Knights. "**Robust statistical models for microbiome phenotype prediction with the MWAS package**" (in preparation)

