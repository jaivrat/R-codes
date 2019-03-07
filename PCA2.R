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
projections = data.matrix %*% pca$rotation
head(projections)
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


