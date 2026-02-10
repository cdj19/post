####### post #######

setGeneric("post", 
           function(model,x1name=NULL,x1vals=NULL,x2name=NULL,x2vals=NULL,holds=NULL,
                    n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2){
             standardGeneric("post")
           }
)

setClassUnion("arrayORNULL", c("array","NULL"))
setClassUnion("listORcharacter", c("list","character"))

setClass("post", 
         slots = c(est = "array", 
                   did = "arrayORNULL", 
                   sims = "array", 
                   model = "character", 
                   link = "listORcharacter", 
                   quantiles = "numeric", 
                   call = "call")
)


post.glm <- function(model,x1name=NULL,x1vals=NULL,x2name=NULL,x2vals=NULL,holds=NULL,
                     n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2){
  
  call <- match.call()
  
  sims <- postSim(model, n.sims=n.sims)
  
  model.link <- family(model)$link
  if (model.link=="identity"){link <- identity} 
  else if (model.link=="probit"){link <- pnorm}
  else if (model.link=="logit"){link <- plogis}
  else if (model.link=="log"){link <- exp}
  else if (model.link=="cloglog"){link <- function(x){1-exp(-exp(x))}}
  else {stop("Link function is not supported")}
  
  n.obs <- nrow(model.matrix(model))
  if (is.null(weights)){wi <- rep(1, n.obs)} else if (length(weights) != n.obs){
    stop("weights must have the same length as the estimation sample")
  } else{wi <- weights}
  k <- length(model.matrix(model)[1,])
  n.q <- length(quantiles)
  
  if (is.null(x1name)){
    X <- array(NA, c(n.obs,k))
    newdata <- data.frame(model$model)
    if (!is.null(holds)){
      for (i in 1:length(holds)){
        newdata[ ,names(holds)[i]] <- holds[[i]]
      }
    }
    X <- aperm(model.matrix(lm(formula(model), data=newdata)))
    l1 <- array(NA, c(nrow(sims@coef),1))
    l1[,1] <- apply(link(sims@coef %*% X), 1, function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(1,n.q+1))
    l2[1,1] <- mean(l1)
    l2[1,2:(n.q+1)] <- quantile(l1, probs=quantiles)
    colnames(l2) <- c("mean",quantiles)
    
    ans <- new("post", 
               est=round(l2, digits=digits),
               did=NULL,
               sims=l1, 
               model=class(model), 
               link=model.link, 
               quantiles=quantiles, 
               call=call)
    return(ans)
  }
  
  else if (is.null(x2name)){
    
    n.x1 <- length(x1vals)
    X <- array(NA, c(n.obs,k,n.x1))
    
    for (i in 1:(n.x1)){
      newdata <- data.frame(model$model)
      if (!is.null(holds)){
        for (j in 1:length(holds)){
          newdata[ ,names(holds)[j]] <- holds[[j]]
        }
      }
      newdata[ ,x1name] <- x1vals[i]
      X[ , ,i] <- model.matrix(lm(formula(model), data=newdata))
    }
    
    X <- aperm(X, c(2,1,3))
    l1 <- apply(apply(X, c(2,3), function(x) link(sims@coef %*% x)), c(1,3), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1))
    l2[1:n.x1,1] <- apply(l1, 2, mean)
    l2[1:n.x1,2:(n.q+1)] <- aperm(apply(l1, 2, function(x) quantile(x, probs=quantiles)))
    
    l2[nrow(l2),1] <- mean(l1[ ,n.x1] - l1[ ,1])
    l2[nrow(l2),2:(n.q+1)] <- quantile(l1[ ,n.x1] - l1[ ,1], probs=quantiles)
    rownames(l2) <- c(paste(c(rep(paste(x1name,"="),n.x1),
                              paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),
                            c(x1vals,"")))  
    colnames(l2) <- c("mean",quantiles)
    
    ans <- new("post", 
               est=round(l2, digits=digits), 
               did=NULL, 
               sims=l1, 
               model=class(model), 
               link=model.link, 
               quantiles=quantiles, 
               call=call)
    return(ans)
  } 
  
  else{
    
    n.x1 <- length(x1vals)
    n.x2 <- length(x2vals)
    X <- array(NA, c(n.obs,k,n.x1,n.x2))
    
    for (j in 1:n.x2){
      for (i in 1:n.x1){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (k in 1:length(holds)){
            newdata[ ,names(holds)[k]] <- holds[[k]]
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        newdata[ ,x2name] <- x2vals[j]
        X[ , ,i,j] <- model.matrix(lm(formula(model), data=newdata))
      }
    }
    
    X <- aperm(X, c(2,1,3,4))
    l1 <- apply(apply(X, c(2,3,4), function(x) link(sims@coef %*% x)), c(1,3,4), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
    l2[1:n.x1,1,1:n.x2] <- apply(l1,c(2,3),mean)
    l2[1:n.x1,2:(n.q+1),1:n.x2] <- aperm(apply(l1, c(2,3), function(x) quantile(x, probs=quantiles)), c(2,1,3))
    l2[nrow(l2),1,1:n.x2] <- apply(l1[ ,n.x1,1:n.x2] - l1[ ,1,1:n.x2], 2, mean)
    l2[nrow(l2),2:(n.q+1),1:n.x2] <- apply(l1[ ,n.x1,1:n.x2] - l1[ ,1,1:n.x2], 2, function(x) quantile(x, probs=quantiles))
    dimnames(l2) <- list(paste(c(rep(paste(x1name,"="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                         c("mean",quantiles),
                         paste(c(rep(paste(x2name,"="),n.x2)),
                               c(x2vals)))
    
    if (is.null(did)){did <- c(x2vals[1],x2vals[n.x2])} else{did <- did}
    l3 <- array(NA, c(1,n.q+1))
    l3[1,1] <- mean((l1[ ,n.x1,match(did[2],x2vals)] - l1[ ,1,match(did[2],x2vals)]) -  (l1[ ,n.x1,match(did[1],x2vals)] - l1[ ,1,match(did[1],x2vals)]))
    l3[1,2:(n.q+1)] <- quantile((l1[ ,n.x1,match(did[2],x2vals)] - l1[ ,1,match(did[2],x2vals)]) -  (l1[ ,n.x1,match(did[1],x2vals)] - l1[ ,1,match(did[1],x2vals)]), probs=quantiles)
    dimnames(l3) <- list("did",c("mean",quantiles)) 
    
    ans <- new("post", 
               est=round(l2, digits=digits), 
               did=round(l3, digits=digits), 
               sims=l1, 
               model=class(model), 
               link=model.link, 
               quantiles=quantiles, 
               call=call)
    return(ans)
  }
}


post.polr <- function(model,x1name=NULL,x1vals=NULL,x2name=NULL,x2vals=NULL,holds=NULL,
         n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2){
  
  call <- match.call()
  
  sims <- suppressMessages(postSim(model, n.sims=n.sims))
  
  n.obs <- length(model$model[,1])
  if (is.null(weights)){wi <- rep(1, n.obs)} else if (length(weights) != n.obs){
    stop("weights must have the same length as the estimation sample")
  } else{wi <- weights}
  
  if (model$method=="probit"){link <- pnorm}
  else if (model$method=="logistic"){link <- plogis}
  else if (model$method=="cloglog"){link <- function(x){1-exp(-exp(x))}}
  else {stop("Link function is not supported")}
  
  k <- length(model.matrix(polr(getCall(model)$formula, model$model))[1,])
  n.q <- length(quantiles)
  n.y <- length(levels(model$model[,1]))
  n.z <- length(model$zeta)
  tau <- array(NA, c(n.sims,n.z+2))
  tau[,1] <- -Inf
  tau[,2:(ncol(tau)-1)] <- sims@zeta[,1:n.z]
  tau[,ncol(tau)] <- Inf
  beta <- sims@coef
  
  if (is.null(cut)){
    
    if (is.null(x1name)){
      
      X_temp <- array(NA, c(n.obs,k))
      X <- array(NA, c(n.obs,k-1))
      
      newdata <- data.frame(model$model)
      if (!is.null(holds)){
        for (j in 1:length(holds)){
          newdata[ ,names(holds)[j]] <- holds[[j]]
        }
      }
      X_temp[ , ] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
      X[ , ] <- X_temp[,-1]
      X <- aperm(X)
      
      l1 <- array(NA, c(n.sims, n.obs, n.y))
      for (z in 1:n.y){
        l1[,,z] <- link(tau[,z+1] - beta %*% X) - link(tau[,z] - beta %*% X)
      }
      
      l2 <- apply(l1, c(1,3), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.y,n.q+1))
      for (i in 1:n.y){
        l3[i,1] <- mean(l2[,i])
        l3[i,2:(n.q+1)] <- quantile(l2[,i], probs=quantiles)
      }
      rownames(l3) <- paste(c(rep("Y =",n.y)), c(1:n.y))
      colnames(l3) <- c("mean",quantiles)
      
      ans <- new("post", 
                 est=round(l3, digits=digits), 
                 did=NULL, 
                 sims=l2, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    }
    
    else if (is.null(x2name)){
      
      n.x1 <- length(x1vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1))
      X <- array(NA, c(n.obs,k-1,n.x1))
      
      for (i in 1:(n.x1)){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (j in 1:length(holds)){
            newdata[ ,names(holds)[j]] <- holds[[j]]
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        X_temp[ , ,i] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.y))
      X <- aperm(X, c(2,1,3))
      for (z in 1:n.y){
        l1[,,,z] <- apply(X, c(2,3), function(x) (link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
      }
      
      l2 <- apply(l1, c(1,3,4), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.x1+1, n.q+1, n.y))
      for (j in 1:n.y){
        for (i in 1:n.x1){
          l3[i,1,j] <- mean(l2[,i,j])
          l3[i,2:(n.q+1),j] <- quantile(l2[,i,j], probs=quantiles)
        }
        l3[nrow(l3),1,j] <- mean(l2[ ,n.x1,j] - l2[ ,1,j])
        l3[nrow(l3),2:(n.q+1),j] <- quantile(l2[ ,n.x1,j] - l2[ ,1,j], probs=quantiles)
      }
      dimnames(l3) <- list(paste(c(rep(paste(x1name,"="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep("Y =",length(levels(model$model[,1])))),
                                 c(1:length(levels(model$model[,1]))))) 
      
      ans <- new("post", 
                 est=round(l3, digits=digits), 
                 did=NULL, 
                 sims=l2, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    }
    
    else{
      
      n.x1 <- length(x1vals)
      n.x2 <- length(x2vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1,n.x2))
      X <- array(NA, c(n.obs,k-1,n.x1,n.x2))
      
      for (j in 1:n.x2){
        for (i in 1:(n.x1)){
          newdata <- data.frame(model$model)
          if (!is.null(holds)){
            for (k in 1:length(holds)){
              newdata[ ,names(holds)[k]] <- holds[[k]]
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
      
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.x2, n.y))
      for (z in 1:n.y){
        l1[,,,,z] <- apply(X, c(2,3,4), function(x) (link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
      }
      
      l2 <- apply(l1, c(1,3,4,5), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.x1+1, n.q+1, n.x2, n.y))
      for (k in 1:n.y){
        for (j in 1:n.x2){
          for (i in 1:n.x1){
            l3[i,1,j,k] <- mean(l2[,i,j,k])
            l3[i,2:(n.q+1),j,k] <- quantile(l2[,i,j,k], probs=quantiles)
          }
          l3[n.x1+1,1,j,k] <- mean(l2[,n.x1,j,k] - l2[,1,j,k])
          l3[n.x1+1,2:(n.q+1),j,k] <- quantile(l2[,n.x1,j,k] - l2[,1,j,k], probs=quantiles)
        }
      }
      dimnames(l3) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep(paste(x2name,"="),n.x2)),x2vals),
                           paste(c(rep("Y =",n.y)), c(1:n.y))) 
      
      if (is.null(did)){did <- c(x2vals[1],x2vals[n.x2])} else{did <- did}
      l4 <- array(NA, c(n.y,n.q+1))
      for (i in 1:n.y){
        l4[i,1] <- mean((l2[ ,n.x1,match(did[2],x2vals),i] - l2[ ,1,match(did[2],x2vals),i]) -  (l2[ ,n.x1,match(did[1],x2vals),i] - l2[ ,1,match(did[1],x2vals),i]))
        l4[i,2:(n.q+1)] <- quantile((l2[ ,n.x1,match(did[2],x2vals),i] - l2[ ,1,match(did[2],x2vals),i]) -  (l2[ ,n.x1,match(did[1],x2vals),i] - l2[ ,1,match(did[1],x2vals),i]), probs=quantiles)
      }
      yvals <- 1:n.y
      dimnames(l4) <- list(paste(c(rep(paste("Y","="),n.y)),yvals),c("mean",quantiles)) 
      
      ans <- new("post", 
                 est=round(l3, digits=digits), 
                 did=round(l4, digits=digits), 
                 sims=l2, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    }
  }
  
  else{
    
    if (is.null(x1name)){
      
      X_temp <- array(NA, c(n.obs,k))
      X <- array(NA, c(n.obs,k-1))
      
      newdata <- data.frame(model$model)
      if (!is.null(holds)){
        for (j in 1:length(holds)){
          newdata[ ,names(holds)[j]] <- holds[[j]]
        }
      }
      X_temp[ , ] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
      X[ , ] <- X_temp[,-1]
      X <- aperm(X)
      
      l1 <- apply(link(-tau[,cut+1] + beta %*% X), 1, function(x) weighted.mean(x, wi))
      l2 <- array(NA, c(1,n.q+1))
      l2[1,1] <- mean(l1)
      l2[1,2:(n.q+1)] <- quantile(l1, probs=quantiles)
      colnames(l2) <- c("mean",quantiles)
      
      ans <- new("post", 
                 est=round(l2, digits=digits), 
                 did=NULL, 
                 sims=l1, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    }
    
    
    else if (is.null(x2name)){
      
      n.x1 <- length(x1vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1))
      X <- array(NA, c(n.obs,k-1,n.x1))
      
      for (i in 1:(n.x1)){
        newdata <- data.frame(model$model)
        if (!is.null(holds)){
          for (j in 1:length(holds)){
            newdata[ ,names(holds)[j]] <- holds[[j]]
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        X_temp[ , ,i] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      X <- aperm(X, c(2,1,3))
      l1 <- apply(apply(X, c(2,3), function(x) link(-tau[,cut+1] + beta %*% x)), 
                  c(1,3), function(x) weighted.mean(x, wi))
      l2 <- array(NA, c(n.x1+1,n.q+1))
      for (i in 1:n.x1){
        l2[i,1] <- mean(l1[,i])
        l2[i,2:(n.q+1)] <- quantile(l1[,i], probs=quantiles)
      }
      l2[nrow(l2),1] <- mean(l1[ ,ncol(l1)] - l1[ ,1])
      l2[nrow(l2),2:(n.q+1)] <- quantile(l1[ ,ncol(l1)] - l1[ ,1], probs=quantiles)
      rownames(l2) <- c(paste(c(rep(paste(x1name,"="),n.x1),
                                paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),
                              c(x1vals,"")))  
      colnames(l2) <- c("mean",quantiles)
      
      ans <- new("post", 
                 est=round(l2, digits=digits), 
                 did=NULL, 
                 sims=l1, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    } 
    
    else{
      
      n.x1 <- length(x1vals)
      n.x2 <- length(x2vals)
      
      X_temp <- array(NA, c(n.obs,k,n.x1,n.x2))
      X <- array(NA, c(n.obs,k-1,n.x1,n.x2))
      
      for (j in 1:n.x2){
        for (i in 1:n.x1){
          newdata <- data.frame(model$model)
          if (!is.null(holds)){
            for (k in 1:length(holds)){
              newdata[ ,names(holds)[k]] <- holds[[k]]
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- suppressWarnings(model.matrix(polr(getCall(model)$formula, data=newdata)))
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
      l1 <- apply(apply(X, c(2,3,4), function(x) link(-tau[,cut+1] + beta %*% x)), c(1,3,4), function(x) weighted.mean(x, wi))
      l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
      for (j in 1:n.x2){
        for (i in 1:n.x1){
          l2[i,1,j] <- mean(l1[,i,j])
          l2[i,2:(n.q+1),j] <- quantile(l1[,i,j], probs=quantiles)
        }
        l2[nrow(l2),1,j] <- mean(l1[ ,n.x1,j] - l1[ ,1,j])
        l2[nrow(l2),2:(n.q+1),j] <- quantile(l1[ ,n.x1,j] - l1[ ,1,j], probs=quantiles)
      }
      dimnames(l2) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep(paste(x2name," ="),n.x2)),
                                 c(x2vals)))
      
      if (is.null(did)){did <- c(x2vals[1],x2vals[n.x2])} else{did <- did}
      l3 <- array(NA, c(1,n.q+1))
      l3[1,1] <- mean((l1[ ,n.x1,match(did[2],x2vals)] - l1[ ,1,match(did[2],x2vals)]) -  (l1[ ,n.x1,match(did[1],x2vals)] - l1[ ,1,match(did[1],x2vals)]))
      l3[1,2:(n.q+1)] <- quantile((l1[ ,n.x1,match(did[2],x2vals)] - l1[ ,1,match(did[2],x2vals)]) -  (l1[ ,n.x1,match(did[1],x2vals)] - l1[ ,1,match(did[1],x2vals)]), probs=quantiles)
      dimnames(l3) <- list("did",c("mean",quantiles)) 
      
      ans <- new("post", 
                 est=round(l2, digits=digits), 
                 did=round(l3, digits=digits), 
                 sims=l1, 
                 model=class(model), 
                 link=model$method, 
                 quantiles=quantiles, 
                 call=call)
      return(ans)
    }
  }
}

setMethod("post", signature(model = "lm"), post.glm)
setMethod("post", signature(model = "glm"), post.glm)
setMethod("post", signature(model = "svyglm"), post.glm)  
setMethod("post", signature(model = "polr"), post.polr)  

    



