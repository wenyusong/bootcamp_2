#Data Handling
#Problem B
apply(my_data, 1, mean)

#Problem C
apply(my_data,1, min)
apply(my_data, 1, max)

#Problem D
gene_sd <- apply(my_data, 1, sd)
head(gene_sd)
max(gene_sd)
order(gene_sd)
sort(gene_sd)
which.max(apply(my_data,1,sd))


#Problem E
max(my_data)
which.max(apply(my_data,1,max))
which.max(apply(my_data,2,max))

#Ploting A
my_data["EGFR",]
hist(as.numeric(my_data["EGFR",]),main="Histogram of EGFR expression", xlab="EGFR expression", ylab="Frequency(# of samples)", breaks=50)

#Ploting B
plot(as.numeric(my_data["EGFR",]), as.numeric(my_data["IDH1",]),main="EGFR versus IDH1", xlab="EGFR expression", ylab="IDH1 expression")

#C
znf_row <- grep("ZNF",rownames(my_data))
grepl("ZNF", rownames(my_data))

znf_row
znf_gene <- my_data[znf_row,]
mean_znf_gene <- apply(znf_gene, 2, mean)
boxplot(mean_znf_gene)
mean_all_gene <- apply(my_data, 2, mean)
boxplot(mean_all_gene)
boxplot(mean_znf_gene,mean_all_gene)

#Function
my_vector <- c(10,5,2,6,8,4,1,9,3,7)
sort_fun <- function(x) {
  for(i in 1:(length(x)-1)){
    imin = i
      for(j in (i+1):length(x)){
        if(x[j]<x[imin]) {
          imin = j
          }
       }
    temp <- x[i]
    x[i] <- x[imin]
    x[imin] <- temp
    }
 return(x)
}
sort_fun(my_vector)

position=numeric()
index_fun <- function(x) {
  for(i in 1:(length(x)-1)){
    imin = i
    for(j in (i+1):length(x)){
      if(x[j]<x[imin]) {
        imin = j
      }
      position <- append(position,imin)
    }
  }
 return(position) 
}

index_fun(my_vector)

library(curl)
my_data <- read.table(curl("https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/unifiedScaledFiltered.txt"),sep="\t",header=1)

vec1 <- c(1:5)
vec2 <- seq(2, 10, 2)

cor(vec1, vec2)
w <- cor(t(my_data))
dim(w)
w2 <- as.data.frame(as.table(w))
head(w2)
dim(w2)


myFunciton <- function(var1, var2){
  for(i in 1:length(var1)){
    
  }
  output <- var1*var2
  return(output)
}
#will's solution
w3 <- w2[w2$Var1!=w2$Var2, ]
w4 <- w3[order(-w3$Freq),]
head(w4, 10)

#D. By default, the R function cor() uses the Pearson method for calculating correlation 
#coefficients. If we change to using a rank-based method like Spearman’s rho, 
#do we find a different pair of genes that are the most highly correlated?
rank_base_cor <- cor(t(my_data), method = "spearman")
rank_base_matrix <- as.data.frame(as.table(rank_base_cor))
rank_base_odered <- rank_base_matrix[rank_base_matrix$Var1!=rank_base_matrix$Var2,]
rank_base_odered <- rank_base_odered[order(-rank_base_odered$Freq),]
head(rank_base_odered, 10)

#yes, the most correlated pair of gene is different

#E. Plot a matrix of the top 20 most correlated genes 
#(hint: consider using the library corrplot)
#install corrplot package
install.packages("corrplot")
library(corrplot)

top20_gene <- rank_base_odered[1:20,]
top20_cor_matirx <- rank_base_cor[top20_gene$Var1,]
corrplot(top20_cor_matirx, method="circle")


#A. Cluster the top 20 most correlated genes using complete hierarchical 
#clustering and visualize the output tree

d <- dist(as.matrix(my_data[top20_gene$Var1,]))
hc <- hclust(d)
plot(hc)

e <- as.matrix(my_data[top20_gene$Var1,])
kc <- kmeans(e, 3)
library(cluster)

clusplot(e, kc$cluster)

#C. If we cluster these genes with complete hierarchical clustering and divide the tree into 4 groups, 
#how many genes are in the largest group? The smallest?
#8 genes
