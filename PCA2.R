#data I gerenaretd using https://www.librec.net/datagen.html, at approx 45 dgrees : Will help verify components
df <- read.csv(text = 
"x,y,lab
3, 3.05, 1
3.6, 3.05, 1
4.35, 3.75, 1
5.7, 5.75, 1
6.2, 4.85, 1
7.9, 5.45, 1
7.85, 6.35, 1
7.45, 7.05, 1
8.9, 7.6, 1
10.4, 8.9, 1
9.85, 7.15, 1
9.9, 7.85, 1
9.35, 8.55, 1
8, 7.8, 1
6.55, 6.55, 1
4.9, 4.95, 1
4.15, 4.75, 1
3.5, 3.9, 1
5.15, 4.3, 1
6.55, 5.6, 1
8.75, 6.75, 1
8.75, 8.25, 1
6.85, 7.35, 1
7.15, 6.55, 1
7.85, 7.05, 1
8.5, 7.3, 1
10.45, 8.3, 1
7.25, 5.55, 1
6.85, 4.85, 1
6, 4.35, 1
5.6, 5.15, 1", header = TRUE)

# Plot to have a sense -------------------------------------------------
head(df)
data.matrix <- as.matrix(df[,c("x","y")])
rownames(data.matrix) <- paste("Sample-", 1:dim(data.matrix)[1], sep="")
#data.matrix <- scale(data.matrix)
plot(data.matrix[,"x"], data.matrix[,"y"])


# PCA using prcomp -----------------------------------------------------
pca         <- prcomp(data.matrix, scale= FALSE, center = FALSE)

# PCA outputs (sdev, x, rotation)---------------------------------------
# sdev: standard deviation along component axis
pca$sdev
(pca.var <- pca$sdev^2)
(pca.var.per <- round(pca.var/sum(pca.var) * 100, 1))
barplot(pca.var.per, xlab = "principal component", ylab="% variation")
# x: pricinpal components : (ie for each sample, we have transformed cordinates along PC1, PC2, ... PCN)
pca$x
plot(pca$x[,1], pca$x[,2])
# rotation: the columns are PC directions. columns are unit vectors along PC directions. THey also will be orthogonal
pca$rotation
apply( pca$rotation, 2, function(e) sum(e^2))#All must sum to 1
#some call these PC vector components as loading scores, because you can take component of original 
#data with loading scores(as a vector dot product and get transfomed vector in terms of PC directions)
loading_scores <- pca$rotation[,1]; loading_scores
#Plot PCA transformed PCA data
pca.data <- data.frame(Sample=rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])
ggplot(data = pca.data, aes(x = X, y = Y, label=Sample)) + 
  geom_point() + 
  #geom_text() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_bw() + 
  ggtitle("My PCA graph")

#if we talke original data series and take dot product of each row on PC1 pf rotation, we will get 
projections.unscaled = data.matrix %*% pca$rotation
head(projections.unscaled)
head(pca$x) #This matches with above if scale=FALSE, center = FALSE in prcomp. or data.matrix not scaled(cetered + scaled)



#SCALE = TRUE ----------------------------------------------------------
##If we use ceter and scaling, then data matrix is scaled and sentered before processing.
pca         <- prcomp(data.matrix, scale= TRUE)
projections = data.matrix %*% pca$rotation
#These wont be same, the two below: because data matrix is scaled and sentered before processing.
head(projections)
head(pca$x)
projections.sc = scale(data.matrix) %*% pca$rotation
#Now these will be equal
head(projections.sc)
head(pca$x)




#########################################################################
# USING SVD
# svd output: 
# svd output is 3 matrix : U D V  s.t  data matrix = A = U D V' 
# U supposed to be eigenvector of AA', 
# D the diagonal mat ofeigenvalues
# V supposed to be eigenvector of A'A 
#########################################################################
#d	:a vector containing the singular values of x, of length min(n, p), sorted decreasingly.
#u	:a matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu).
#v	:a matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension c(p, nv).

res.svd <- svd(data.matrix)
#res.svd$u %*% diag(res.svd$d) %*% t(res.svd$v)  should be same as data.matrix
#Let us verify
sum(apply(res.svd$u %*% diag(res.svd$d) %*% t(res.svd$v) - data.matrix, 1, function(e) sum(e^2)))

#if A be data matrix (m x n), U supposed to be eigenvector of AA'
#eigen(data.matrix %*% t(data.matrix))$vectors to be same as res.svd$u
#So, the SVD result and eigen vector result must be parallel (dot product should 1 or -1 here)
aat <- data.matrix %*% t(data.matrix)
sum(eigen(aat)$vectors[,1] * res.svd$u[,1])
sum(eigen(aat)$vectors[,2] * res.svd$u[,2])
#similarly for v
ata <- t(data.matrix) %*% data.matrix
sum(eigen(ata)$vectors[,1] * res.svd$v[,1])
sum(eigen(ata)$vectors[,2] * res.svd$v[,2])

#Eigen values from AA' or A'A will be squared of svd
# This result will be almost 0
eigen(ata)$value  - (res.svd$d)^2  #Note that while taking variance we need to do res.svd$d ^2 before

#A = UDV' => AV = UD   => it means individually A * v[,i] = d[i] * u[,i] 
#ie. V is like rotation matrix : THe below must all be 0
summary(data.matrix %*% res.svd$v[,1]  - res.svd$d[1] * res.svd$u[,1])

#Find components:  data * rotation, ie AV is new projection.
head(data.matrix %*% res.svd$v)
head(projections.unscaled)
# Or AV -> UD is new projection. AV is actual full projection
# if we want top reduce dimension, then we can drop the dimestions
# from d matrix to get projection
head(res.svd$u %*% diag(res.svd$d))

#If we want projections to mxn data to mxk, where k < n ie dimension reduction,
# we can take first 
svd.vars <- eigen(ata)$value
svd.cum.vars <- cumsum(svd.vars)/sum(svd.vars) #%variation
svd.cum.vars
#[1] 0.9966629 1.0000000
#We see above that 99% variation is explained by first component, => k = 1
#So, we can take first k (k=1) components of UD as PC components
k = 1 #in this case
pccomps <- res.svd$u %*% diag(res.svd$d)[,1:k, drop=FALSE]
head(pccomps) #ie first column of head(res.svd$u %*% diag(res.svd$d))


#Suppose we are given components, say pccomps above, how do we get back the original data (with some)
# We will definitely need rotation matrix
# AV = UD, so this means UD is GIVEN, we recover A from equation AV = GIVEN, ie A = GIVEN * V' 
n = dim(data.matrix)[2]
recovered = cbind(pccomps, rep(0, ncol=n-k)) %*% t(res.svd$v)
colnames(recovered) <- c("x","y")
head(recovered)
head(data.matrix)

#Plot to have a look 
tmp <- rbind(data.frame(recovered, label= "recovered"), data.frame(data.matrix,label= "original" ))
ggplot(tmp, aes(x=x, y=y, color=label)) + geom_point()
