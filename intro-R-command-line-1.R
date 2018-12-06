data(iris)
print(iris)

write.csv(x = iris, file = "data2/iris.csv", quote = F, row.names = F)