d = data.frame(x=c(2.5,0.5,2.2,1.9,3.1,2.3,2.0,1.0,1.5,1.1),
               y=c(2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9))

d

# mean-adjusted values 
d$x_adj = d$x - mean(d$x)
d$y_adj = d$y - mean(d$y)


# calculate covariance matrix and eigenvectors/values
(cm = cov(d[,1:2]))
(e = eigen(cm))

# principal component vector slopes
s1 = e$vectors[1,1] / e$vectors[2,1] # PC1
s2 = e$vectors[1,2] / e$vectors[2,2] # PC2

plot(d$x_adj, d$y_adj, asp=T, pch=16, xlab='x', ylab='y')
abline(a=0, b=s1, col='red')
abline(a=0, b=s2)

# PCA data = rowFeatureVector (transposed eigenvectors) * RowDataAdjust (mean adjusted, also transposed)
feat_vec = t(e$vectors)
row_data_adj = t(d[,3:4])
final_data = data.frame(t(feat_vec %*% row_data_adj)) # ?matmult for details
names(final_data) = c('x','y')


plot(final_data, asp=T, xlab='PCA 1', ylab='PCA 2', pch=16)

#First principal component v(dot product)e1, second component v(dot product)e1
as.matrix(d[,3:4]) %*% e$vectors[,1, drop = FALSE] #these are componnet on line e1
compOnE1 <- e$vectors[,1, drop = FALSE] %*% t(as.matrix(d[,3:4]) %*% e$vectors[,1, drop = FALSE])
compOnE1 <- data.frame(t(compOnE1))

compOnE2 <- e$vectors[,2, drop = FALSE] %*% t(as.matrix(d[,3:4]) %*% e$vectors[,2, drop = FALSE])
compOnE2 <- data.frame(t(compOnE2))

plot(d$x_adj, d$y_adj, asp=T, pch=1, xlab='x', ylab='y', col = 'black')
abline(a=0, b=s1, col='red')
abline(a=0, b=s2)
points(compOnE1, asp=T, pch=16, col='red')
points(compOnE2, asp=T, pch=16, col='green')
