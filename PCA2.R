#data I gerenaretd using https://www.librec.net/datagen.html, at approx 45 dgrees
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

head(df)
data.matrix <- as.matrix(df[,c("x","y")])
#data.matrix <- scale(data.matrix)
rownames(data.matrix) <- paste("Sample-", 1:dim(data.matrix)[1], sep="")
pca         <- prcomp(data.matrix, scale= FALSE, center = FALSE)
pca$sdev
pca$x
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var) * 100, 1)
pca.var.per
barplot(pca.var.per, xlab = "principal component", ylab="% variation")

pca.data <- data.frame(Sample=rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])
pca.data


ggplot(data = pca.data, aes(x = X, y = Y, label=Sample)) + 
  geom_point() + 
  #geom_text() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_bw() + 
  ggtitle("My PCA graph")

loading_scores <- pca$rotation[,1]; loading_scores
#All must sum to 1
apply( pca$rotation, 2, function(e) sum(e^2))

#if we talke original data series and take dot product of each row on PC1 pf rotation, we will get 
projections = data.matrix %*% pca$rotation
#means <- apply(projections,2, mean)
#head(apply(projections, 2, function(e) e - means))
head(projections)
head(pca$x)

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


