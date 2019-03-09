d = data.frame(x=c(2.5,0.5,2.2,1.9,3.1,2.3,2.0,1.0,1.5,1.1),
               y=c(2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9))

d

# mean-adjusted values 
d$x_adj = d$x - mean(d$x)
d$y_adj = d$y - mean(d$y)
d

# calculate covariance matrix and eigenvectors/values
(cm = cov(d[,1:2]))
(e = eigen(cm))

#----- Checks
#Egen vectors are reported as columns of matrix e$vectors
# First -> e$vectors[,1], second e$vectors[,2]
cm %*% e$vectors[,1]
e$values[1] * e$vectors[,1]

cm %*% e$vectors[,2]
e$values[2] * e$vectors[,2]
#------


# principal component vector slopes
s1 = e$vectors[2,1] / e$vectors[1,1] # PC1
s2 = e$vectors[2,2] / e$vectors[1,2] # PC2

plot(d$x_adj, d$y_adj, asp=T, pch=16, xlab='x', ylab='y')
abline(a=0, b=s1, col='red')
abline(a=0, b=s2)

# PCA data = rowFeatureVector (transposed eigenvectors) * RowDataAdjust (mean adjusted, also transposed)
final_data   = data.frame(as.matrix(d[,3:4]) %*% e$vectors)
names(final_data) = c('x','y')


plot(final_data, asp=T, xlab='PCA 1', ylab='PCA 2', pch=16)

#First principal component v(dot product)e1, second component v(dot product)e1
as.matrix(d[,3:4]) %*% e$vectors[,1, drop = FALSE] #these are componnet on line e1 (ACTUAL COMPONENT)
#but we want to draw projections (red and green points)
#so we multiply with vector e$vectors[,1, drop = FALSE] before to get the projections
compOnE1 <- e$vectors[,1, drop = FALSE] %*% t(as.matrix(d[,3:4]) %*% e$vectors[,1, drop = FALSE])
compOnE1 <- data.frame(t(compOnE1))
#Similarly projection on E2
compOnE2 <- e$vectors[,2, drop = FALSE] %*% t(as.matrix(d[,3:4]) %*% e$vectors[,2, drop = FALSE])
compOnE2 <- data.frame(t(compOnE2))

plot(d$x_adj, d$y_adj, asp=T, pch=1, xlab='x', ylab='y', col = 'black')
abline(a=0, b=s1, col='red')
abline(a=0, b=s2)
points(compOnE1, asp=T, pch=16, col='red')
points(compOnE2, asp=T, pch=16, col='green')

#================================= SVD =========================================+#
svd.res <- svd(d[,1:2])
e$values
e$vectors

# A = UDV'
# U : eigen vectors of AA'  and V: eigen vectors of A'A
# D = diag(sqrt eigenvalues of A'A)
aat = as.matrix(d[,1:2]) %*% t(as.matrix(d[,1:2]))
aat.eig <- eigen(aat)
aat.eig$values[abs(aat.eig$values) < 1e-10] <- 0
expected.D = sqrt(aat.eig$values)
expected.U = ata.eig$vectors
#seems it is slightl;y different from e$vectors
#Slopes must match
expected.U[2,1]/expected.U[1,1] - s1
expected.U[2,2]/expected.U[1,2]


#A'A
ata.eig   =  eigen( t(as.matrix(d[,1:2])) %*% as.matrix(d[,1:2]))
(expected.D.new = sqrt(aat.eig$values))
expected.V.new = ata.eig$vectors
expected.V.new

eigen(t(as.matrix(d[,1:2])) %*% as.matrix(d[,1:2]))


