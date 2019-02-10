#!/usr/bin/env Rscript


packages <- c("ape", "ggplot2")

# Function to check whether package is installed
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

for(x in packages){
  if(!is.element(x, installed.packages()[,1]))
    {install.packages(x, repos="http://cran.fhcrc.org")
  } else {print(paste(x, " library already installed"))}
}

