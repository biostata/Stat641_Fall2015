## To choose a CRAN Mirror
#chooseCRANmirror(graphics=FALSE)
## Choose 57 for Berkeley, 60 for Iowa, 68 for Revolution R

#install.packages("rneos", repos="http://r-forge.r-project.org/projects/rneos/")

fn_getPackages <- function(){
  Packages <- c(
    'alr3', 
    'animint',            ## For animated ggplot2 graphics
    'biganalytics', 
    'biglm', 
    'bigmemory', 
    'boot',  
    'car', 
    'chron', 
    'cluster', 
    'clusterSim',  
    'clValid', 
    'DAAG',
    #'devtools',           ## Only works on R3.2.2 or higher
    'doSNOW', 
    'dtw', 
    'epicalc', 
    'faraway', 
    'fAssets', 
    'fBasics', 
    'fda', 
    'fdakma',
    'fdasrvf',
    'fields', 
    'fOptions',
    'foreach',
    'forecast',            ## for time series modeling  
    'foreign', 
    'fpc',                 ## For cluster comparisons ## Wont install on lmcg
#     'fPortfolio', 
#     'fRegression', 
#     'fSeries', 
#     'fUtilities', 
    'gam', 
    'gdata', 
    'glmnet', 
#     'gmaps', 
    'gmodels', 
    'ggplot2',  
    'graphics', 
    'grid',                ## Needed by animint
    'gridExtra', 
    'gsubfn', 
    'gtools', 
    'hexbin',             ## Needed by animint
    'HiddenMarkov', 
    'HMM', 
    'Hmisc', 
    'hydroTSM', 
    'inline', 
    'iterators', 
#   	'jpeg',			## To read/write jpg images
    'kernlab', 
    'knitr', 
    'leaps', 
    'limma', 
    'lme4', 
    'lmtest',  
    'longitudinalData', 
    'lpSolve', 
    'lubridate', 
    'maps', 
    'maptools', 
    'mAr', 
    'matrixStats', 
    'mclust', 
    'MEMSS', 
    'mgcv', 
    'modeest', 
    'nlme', 
    'nortest', 
    'numDeriv',           ## To calculate gradient, hessian, jacobian of functions
    'nws',
    'orthogonalsplinebasis',
    'pdfCluster',   	    ## Cluster pdfs from mixtures, non-parametric
    'PerformanceAnalytics', 
    'phangorn', 
    'plotrix', 
    'pls', 
    'plyr',  
#     'png',                ## To read and write png images
    'profr',  
    'proto',              ## Needed by animint
    'psych', 
    'quantmod', 
    'quantreg', 
    'R2HTML', 
    'rbenchmark', 
    'Rcpp',
    'readr',  ## Read flat/tabular text files from disk
    'reshape',  
    'reshape2',  
    'RMySQL',   
    'robustbase', 
    'robustX',             ## For multivariate median and other experimental stats
#     'RJSONIO',             ## Required by animint
#     'Rsymphony', 
    'sandwich', 
    'scales',              ## Needed by animint
    'SenSrivastava',
    'servr',               ## To plot html & js objects produced by animint
    'sets',                ## Set operations      
    'slam',   
    'snow', 
    'sos',
    'speff2trial',
    'st', 
    'stats', 
    'survey', 
    'svMisc', 
    'tseries', 
    'TTR', 
    'urca',
    'UsingR',    
    'vrtest', 
    'wle', 
#     'wordcloud', 
#     'WriteXLS',
    'xgboost',  ## Extreme Gradient Boosting
#     'XML', 
    'xtable', 
    'xts', 
    'zipcode', 
    'zoo'
  )
  return(Packages)  
}

fn_Install_Packages_CRAN <- function(Packages = fn_getPackages()){
  ## Choose USA (IA) as the CRAN mirror
  Mirrors <- getCRANmirrors(all = FALSE, local.only = TRUE)
  chooseCRANmirror(graphics = F, ind = which(Mirrors$Name == 'USA (IA)'))

  for(Package in Packages){
    if(require(package=Package, character.only=T) == F){
      try(install.packages(Package, dependencies = TRUE))
    } else{
      print(paste(Package, 'already exists'))
    }
  }
}

#fn_Install_Packages_CRAN(Packages = fn_getPackages())
