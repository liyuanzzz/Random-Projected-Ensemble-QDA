QDAbag <- function(xtrain,ytrain,xtest,ytest,B,p){
  n1.train <- nrow(xtrain); n2.train <- nrow(ytrain)
  n1.test <- nrow(xtest); n2.test <- nrow(ytest)
  
  #ztest <- rbind(xtest,ytest)
  labelz <- c(rep(1,n1.test),rep(2,n2.test))
  n=n1.test+n2.test
  IDX_QDAbb <- matrix(0,nrow=n,ncol=B)
  IDX_QDA <- rep(0,n)
  
  
  c1.train <- p/n1.train; c2.train <- p/n2.train
  pori <- ncol(xtrain)
  
  for(bb in 1:B){
    
    #random projection
    Q=qr.Q(qr(matrix(rnorm(pori*p),pori)))
    R=qr.Q(qr(matrix(rnorm(pori*p),pori)))
    L=diag(diag(R)/abs(diag(R)))
    R1=Q%*%L
    xt <- xtrain%*%R1; yt <- ytrain%*%R1
    ztest <- rbind(xtest%*%R1,ytest%*%R1)
    
    #pindex <- sample(1:pori,p)
    #xt <- xtrain[,pindex]; yt <- ytrain[,pindex]
    
    
    hatmux <- colMeans(xt)
    hatmuy <- colMeans(yt)
    
    hatSigmaX <- cov(xt)
    hatSigmaY <- cov(yt)
    error01 <- (hatmux%*%inv(hatSigmaX)%*%hatmux)*(1-c1.train)-(hatmuy%*%inv(hatSigmaY)%*%hatmuy)*(1-c2.train)
    error04 <- log(det(hatSigmaX))-log(det(hatSigmaY))
    
    
    for (i in c(1:n)){
      z <- ztest[i,]
      error02 <- (z%*%inv(hatSigmaX)%*%(z))*(1-c1.train)-(z%*%inv(hatSigmaY)%*%(z))*(1-c2.train)
      #z <- ztest[i,]%*%R1
      #error02 <- (z%*%inv(hatSigmaX)%*%t(z))*(1-c1.train)-(z%*%inv(hatSigmaY)%*%t(z))*(1-c2.train)
      error03 <- (z%*%inv(hatSigmaX)%*%hatmux)*(1-c1.train)-(z%*%inv(hatSigmaY)%*%hatmuy)*(1-c2.train)
      #regori <- t(z-mu1)
      IDX_QDAbb[i,bb]<- -(error04+error02-2*error03+error01)/p
      #IDX_oriQDA[i]<- (regori)%*%(inv(Sigma2)-inv(Sigma1))%*%t(regori)-2*(mu2-mu1)%*%inv(Sigma2)%*%(z-mu1/2-mu2/2)-log(det(Sigma1))+log(det(Sigma2))
    }
  }
  IDX_QDA <- rowMeans(IDX_QDAbb)
  IDX_QDA <- (IDX_QDA<=1e-06)+1
  
  
  eL <- round(sum(abs(IDX_QDA-labelz))/n,4)
  #eO <- round(sum(abs(IDX_oriQDA-labelz))/n,4)
  return(eL)
  
}

LDAbag <- function(xtrain,ytrain,xtest,ytest,B,p){
  n1.train <- nrow(xtrain); n2.train <- nrow(ytrain)
  n1.test <- nrow(xtest); n2.test <- nrow(ytest)
  
  #ztest <- rbind(xtest,ytest)
  labelz <- c(rep(1,n1.test),rep(2,n2.test))
  n=n1.test+n2.test
  IDX_LDAbb <- matrix(0,nrow=n,ncol=B)
  IDX_LDA <- rep(0,n)
  
  
  c1.train <- p/n1.train; c2.train <- p/n2.train
  pori <- ncol(xtrain)
  
  for(bb in 1:B){
    
    #random projection
    Q=qr.Q(qr(matrix(rnorm(pori*p),pori)))
    R=qr.Q(qr(matrix(rnorm(pori*p),pori)))
    L=diag(diag(R)/abs(diag(R)))
    R1=Q%*%L
    xt <- xtrain%*%R1; yt <- ytrain%*%R1
    ztest <- rbind(xtest%*%R1,ytest%*%R1)
    
    #pindex <- sample(1:pori,p)
    #xt <- xtrain[,pindex]; yt <- ytrain[,pindex]
    
    
    hatmux <- colMeans(xt)
    hatmuy <- colMeans(yt)
    
    #hatSigmaX <- cov(xt)
    #hatSigmaY <- cov(yt)
    hatSigmaXY <- cov(rbind(xt,yt))
    error01 <- (hatmux%*%inv(hatSigmaXY)%*%hatmux)*(1-c1.train)-(hatmuy%*%inv(hatSigmaXY)%*%hatmuy)*(1-c2.train)
    error04 <- 0
    
    
    for (i in c(1:n)){
      z <- ztest[i,]
      error02 <- (z%*%inv(hatSigmaXY)%*%(z))*(1-c1.train)-(z%*%inv(hatSigmaXY)%*%(z))*(1-c2.train)
      #z <- ztest[i,]%*%R1
      #error02 <- (z%*%inv(hatSigmaX)%*%t(z))*(1-c1.train)-(z%*%inv(hatSigmaY)%*%t(z))*(1-c2.train)
      error03 <- (z%*%inv(hatSigmaXY)%*%hatmux)*(1-c1.train)-(z%*%inv(hatSigmaXY)%*%hatmuy)*(1-c2.train)
      #regori <- t(z-mu1)
      IDX_LDAbb[i,bb]<- -(error04+error02-2*error03+error01)/p
      #IDX_oriQDA[i]<- (regori)%*%(inv(Sigma2)-inv(Sigma1))%*%t(regori)-2*(mu2-mu1)%*%inv(Sigma2)%*%(z-mu1/2-mu2/2)-log(det(Sigma1))+log(det(Sigma2))
    }
  }
  
  IDX_LDA <- rowMeans(IDX_LDAbb)
  IDX_LDA <- (IDX_LDA<=1e-06)+1
  
  eL <- round(sum(abs(IDX_LDA-labelz))/n,4)
  #eO <- round(sum(abs(IDX_oriQDA-labelz))/n,4)
  return(eL)
  
}

