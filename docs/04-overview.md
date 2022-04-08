# Introduction

## The difficulties of high-dimensional datasets

Biological data, especially molecular data, is inherently high-dimensional. For example, expression matrices can have data for tens of thousands of genes, and each gene represents a different dimension of the data. The sample sizes necessary to fit statistical models to such a high-dimensional dataset are astronomically high and impractical for most modern studies. High-dimensional data also greatly increases the computational resources and time necessary for any bioinformatic analysis.


## What is dimensionality reduction?

Dimensionality reduction methods allow scientists to get around these issues. They are techniques that reduce the number of input variables in a dataset. These methods are possible for biological data because so much of the data variables are correlated. If we go back to our example using expression data, while each gene represents a different dimension, genes involved in similar metabolic pathways will be upregulated or downregulated together and thus those dimensions are not independent. Dimensional reduction methods can identify these correlations and combine them in such a way to reduce the number of data dimensions to a workable amount while minimizing the information loss. As an added benefit, combining correlated dimensions can result in less noisy data, since the data are averaged across multiple dimensions to provide a better representation of patterns. Dimensional reduction also makes it easier for scientists to visualize data in ways that are comprehensible to the human brain. 

