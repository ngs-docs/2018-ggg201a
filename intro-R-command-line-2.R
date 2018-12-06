#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

data(iris)
print(iris)

write.csv(x = iris, file = args[1], quote = F, row.names = F)