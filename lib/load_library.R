"load.library" <- function (myPackages)  {   
  
  # start loop to determine if each package is installed
  for(package in myPackages){
    
    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))
    
    # if package is not installed locally, download, then load
    else {
      install.packages(package)
      do.call("library", list(package))
    }
  } 
}