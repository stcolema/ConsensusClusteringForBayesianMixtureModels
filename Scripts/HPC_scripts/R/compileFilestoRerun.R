#!/bin/Rscript

# Compile the various files to rerun into one data.frame

library(stringr)
library(dplyr)

files=list.files("./", pattern=".csv$")
files=str_sort(files, numeric=T)

n_files=length(files)

my_data=list()
list_ind=1
for(i in 1:n_files){
  
  x = read.csv(files[i], header = T, row.names=1)
  if(!is.null(x) & nrow(x) != 0){
    x$Simulation = i
    x$Rerun = as.character(x$Rerun)
    my_data[[list_ind]] = x
    list_ind=list_ind+1
  }
}
if(list_ind==1){
  cat("\nNo files to rerun.")
  stop("No problems.")
}

file_df = my_data %>% 
  unlist() %>% 
  matrix(ncol=2, byrow = T) %>% 
  as.data.frame()

colnames(file_df)= c("Rerun", "Simulation")
print(file_df)
