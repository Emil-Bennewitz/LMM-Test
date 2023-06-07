
# Create a test block matrix to test swap.block, swap.column on
M1 <- matrix(1:16, nrow = 4, ncol = 4)
M2 <- matrix(1:9 , nrow=3)
M3 <- matrix(1:4, nrow=2)
test.matrix<- bdiag(M1,M2,M3)
test.matrix.swap<-as.matrix(swap.block(test.matrix,5,7))
as.matrix(swap.column(test.matrix,5 ,7))