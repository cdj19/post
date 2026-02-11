####### post #######

setGeneric("post",
           function(model,x1name=NULL,x1vals=NULL,x2name=NULL,x2vals=NULL,holds=NULL,
                    n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2,
                    dist=c("normal","t")){
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
                     n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2,
                     dist=c("normal","t")){

  call <- match.call()
  dist <- match.arg(dist)

  sims <- postSim(model, n.sims=n.sims, dist=dist)
  
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

  if (!is.null(x1name) && !x1name %in% names(model$model)) {
    stop(sprintf("x1name='%s' is not a variable in the model", x1name))
  }
  if (!is.null(x2name) && !x2name %in% names(model$model)) {
    stop(sprintf("x2name='%s' is not a variable in the model", x2name))
  }

  if (!is.null(x1name) && grepl("(", x1name, fixed = TRUE)) {
    stop(sprintf("x1name='%s' appears to use an inline transformation. Define the variable in your data before fitting the model (e.g., data$x <- factor(x)).", x1name))
  }
  if (!is.null(x2name) && grepl("(", x2name, fixed = TRUE)) {
    stop(sprintf("x2name='%s' appears to use an inline transformation. Define the variable in your data before fitting the model (e.g., data$x <- factor(x)).", x2name))
  }

  if (!is.null(x1name) && is.factor(model$model[[x1name]])) {
    x1vals <- as.character(x1vals)
    if (!all(x1vals %in% levels(model$model[[x1name]]))) {
      stop(sprintf("x1vals must be valid levels of factor '%s': %s",
                   x1name, paste(levels(model$model[[x1name]]), collapse = ", ")))
    }
  }
  if (!is.null(x2name) && is.factor(model$model[[x2name]])) {
    x2vals <- as.character(x2vals)
    if (!all(x2vals %in% levels(model$model[[x2name]]))) {
      stop(sprintf("x2vals must be valid levels of factor '%s': %s",
                   x2name, paste(levels(model$model[[x2name]]), collapse = ", ")))
    }
  }

  if (!is.null(holds)) {
    bad <- setdiff(names(holds), names(model$model))
    if (length(bad) > 0) {
      stop(sprintf("holds variable(s) not found in model: %s",
                   paste(bad, collapse = ", ")))
    }
    for (i in seq_along(holds)) {
      nm <- names(holds)[i]
      if (is.factor(model$model[[nm]])) {
        val <- as.character(holds[[i]])
        if (!val %in% levels(model$model[[nm]])) {
          stop(sprintf("holds value '%s' is not a valid level of factor '%s': %s",
                       val, nm, paste(levels(model$model[[nm]]), collapse = ", ")))
        }
      }
    }
  }

  if (is.null(x1name)){
    X <- array(NA, c(n.obs,k))
    newdata <- data.frame(model$model)
    if (!is.null(holds)){
      for (i in 1:length(holds)){
        newdata[ ,names(holds)[i]] <- holds[[i]]
      }
    }
    X <- aperm(model.matrix(terms(model), data=newdata, xlev=model$xlevels))
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
      X[ , ,i] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
    }
    
    X <- aperm(X, c(2,1,3))
    l1 <- apply(apply(X, c(2,3), function(x) drop(link(sims@coef %*% x))), c(1,3), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1))
    l2[1:n.x1,1] <- apply(l1, 2, mean)
    q_arr <- apply(l1, 2, function(x) quantile(x, probs=quantiles))
    if (n.q == 1) q_arr <- array(q_arr, dim = c(1, length(q_arr)))
    l2[1:n.x1,2:(n.q+1)] <- aperm(q_arr)
    
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
          for (h in 1:length(holds)){
            newdata[ ,names(holds)[h]] <- holds[[h]]
          }
        }
        newdata[ ,x1name] <- x1vals[i]
        newdata[ ,x2name] <- x2vals[j]
        X[ , ,i,j] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
      }
    }
    
    X <- aperm(X, c(2,1,3,4))
    l1 <- apply(apply(X, c(2,3,4), function(x) drop(link(sims@coef %*% x))), c(1,3,4), function(x) weighted.mean(x, wi))
    l2 <- array(NA, c(n.x1+1,n.q+1,n.x2))
    l2[1:n.x1,1,1:n.x2] <- apply(l1,c(2,3),mean)
    q_arr <- apply(l1, c(2,3), function(x) quantile(x, probs=quantiles))
    if (n.q == 1) dim(q_arr) <- c(1, dim(q_arr))
    l2[1:n.x1,2:(n.q+1),1:n.x2] <- aperm(q_arr, c(2,1,3))
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
         n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2,
         dist=c("normal","t")){

  call <- match.call()
  dist <- match.arg(dist)

  sims <- suppressMessages(postSim(model, n.sims=n.sims, dist=dist))
  
  n.obs <- length(model$model[,1])
  if (is.null(weights)){wi <- rep(1, n.obs)} else if (length(weights) != n.obs){
    stop("weights must have the same length as the estimation sample")
  } else{wi <- weights}
  
  if (model$method=="probit"){link <- pnorm}
  else if (model$method=="logistic"){link <- plogis}
  else if (model$method=="cloglog"){link <- function(x){1-exp(-exp(x))}}
  else {stop("Link function is not supported")}
  
  k <- ncol(model.matrix(terms(model), data=model$model, xlev=model$xlevels))
  n.q <- length(quantiles)
  n.y <- length(levels(model$model[,1]))
  n.z <- length(model$zeta)
  tau <- array(NA, c(n.sims,n.z+2))
  tau[,1] <- -Inf
  tau[,2:(ncol(tau)-1)] <- sims@zeta[,1:n.z]
  tau[,ncol(tau)] <- Inf
  beta <- sims@coef

  if (!is.null(x1name) && !x1name %in% names(model$model)) {
    stop(sprintf("x1name='%s' is not a variable in the model", x1name))
  }
  if (!is.null(x2name) && !x2name %in% names(model$model)) {
    stop(sprintf("x2name='%s' is not a variable in the model", x2name))
  }

  if (!is.null(x1name) && grepl("(", x1name, fixed = TRUE)) {
    stop(sprintf("x1name='%s' appears to use an inline transformation. Define the variable in your data before fitting the model (e.g., data$x <- factor(x)).", x1name))
  }
  if (!is.null(x2name) && grepl("(", x2name, fixed = TRUE)) {
    stop(sprintf("x2name='%s' appears to use an inline transformation. Define the variable in your data before fitting the model (e.g., data$x <- factor(x)).", x2name))
  }

  if (!is.null(x1name) && is.factor(model$model[[x1name]])) {
    x1vals <- as.character(x1vals)
    if (!all(x1vals %in% levels(model$model[[x1name]]))) {
      stop(sprintf("x1vals must be valid levels of factor '%s': %s",
                   x1name, paste(levels(model$model[[x1name]]), collapse = ", ")))
    }
  }
  if (!is.null(x2name) && is.factor(model$model[[x2name]])) {
    x2vals <- as.character(x2vals)
    if (!all(x2vals %in% levels(model$model[[x2name]]))) {
      stop(sprintf("x2vals must be valid levels of factor '%s': %s",
                   x2name, paste(levels(model$model[[x2name]]), collapse = ", ")))
    }
  }

  if (!is.null(holds)) {
    bad <- setdiff(names(holds), names(model$model))
    if (length(bad) > 0) {
      stop(sprintf("holds variable(s) not found in model: %s",
                   paste(bad, collapse = ", ")))
    }
    for (i in seq_along(holds)) {
      nm <- names(holds)[i]
      if (is.factor(model$model[[nm]])) {
        val <- as.character(holds[[i]])
        if (!val %in% levels(model$model[[nm]])) {
          stop(sprintf("holds value '%s' is not a valid level of factor '%s': %s",
                       val, nm, paste(levels(model$model[[nm]]), collapse = ", ")))
        }
      }
    }
  }

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
      X_temp[ , ] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
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
      rownames(l3) <- paste("Y =", levels(model$model[,1]))
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
        X_temp[ , ,i] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.y))
      X <- aperm(X, c(2,1,3))
      for (z in 1:n.y){
        l1[,,,z] <- apply(X, c(2,3), function(x) drop(link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
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
                           paste("Y =", levels(model$model[,1])))
      
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
            for (h in 1:length(holds)){
              newdata[ ,names(holds)[h]] <- holds[[h]]
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
      
      l1 <- array(NA, c(n.sims, n.obs, n.x1, n.x2, n.y))
      for (z in 1:n.y){
        l1[,,,,z] <- apply(X, c(2,3,4), function(x) drop(link(tau[,z+1] - beta %*% x) - link(tau[,z] - beta %*% x)))
      }
      
      l2 <- apply(l1, c(1,3,4,5), function(x) weighted.mean(x, wi))
      l3 <- array(NA, c(n.x1+1, n.q+1, n.x2, n.y))
      for (m in 1:n.y){
        for (j in 1:n.x2){
          for (i in 1:n.x1){
            l3[i,1,j,m] <- mean(l2[,i,j,m])
            l3[i,2:(n.q+1),j,m] <- quantile(l2[,i,j,m], probs=quantiles)
          }
          l3[n.x1+1,1,j,m] <- mean(l2[,n.x1,j,m] - l2[,1,j,m])
          l3[n.x1+1,2:(n.q+1),j,m] <- quantile(l2[,n.x1,j,m] - l2[,1,j,m], probs=quantiles)
        }
      }
      dimnames(l3) <- list(paste(c(rep(paste(x1name," ="),n.x1),paste("\u0394","(",x1vals[1],",",x1vals[length(x1vals)],")")),c(x1vals,"")),
                           c("mean",quantiles),
                           paste(c(rep(paste(x2name,"="),n.x2)),x2vals),
                           paste("Y =", levels(model$model[,1])))
      
      if (is.null(did)){did <- c(x2vals[1],x2vals[n.x2])} else{did <- did}
      l4 <- array(NA, c(n.y,n.q+1))
      for (i in 1:n.y){
        l4[i,1] <- mean((l2[ ,n.x1,match(did[2],x2vals),i] - l2[ ,1,match(did[2],x2vals),i]) -  (l2[ ,n.x1,match(did[1],x2vals),i] - l2[ ,1,match(did[1],x2vals),i]))
        l4[i,2:(n.q+1)] <- quantile((l2[ ,n.x1,match(did[2],x2vals),i] - l2[ ,1,match(did[2],x2vals),i]) -  (l2[ ,n.x1,match(did[1],x2vals),i] - l2[ ,1,match(did[1],x2vals),i]), probs=quantiles)
      }
      dimnames(l4) <- list(paste("Y =", levels(model$model[,1])),c("mean",quantiles))
      
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
      X_temp[ , ] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
      X[ , ] <- X_temp[,-1]
      X <- aperm(X)
      
      l1 <- apply(1 - link(tau[,cut+1] - beta %*% X), 1, function(x) weighted.mean(x, wi))
      l1 <- array(l1, dim = c(length(l1), 1))
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
        X_temp[ , ,i] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
        X[ , ,i] <- X_temp[,-1,i]
      }
      
      X <- aperm(X, c(2,1,3))
      l1 <- apply(apply(X, c(2,3), function(x) drop(1 - link(tau[,cut+1] - beta %*% x))),
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
            for (h in 1:length(holds)){
              newdata[ ,names(holds)[h]] <- holds[[h]]
            }
          }
          newdata[ ,x1name] <- x1vals[i]
          newdata[ ,x2name] <- x2vals[j]
          X_temp[ , ,i,j] <- model.matrix(terms(model), data=newdata, xlev=model$xlevels)
          X[ , ,i,j] <- X_temp[,-1,i,j]
        }
      }
      
      X <- aperm(X, c(2,1,3,4))
      l1 <- apply(apply(X, c(2,3,4), function(x) drop(1 - link(tau[,cut+1] - beta %*% x))), c(1,3,4), function(x) weighted.mean(x, wi))
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

    



