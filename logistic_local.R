library(MASS)
library(Matrix)
library(cfzip)
library(sets)
set.seed(42)
setwd(dir = '/home/liyuelin37/result/sublogit_sim4/n=500_mu=1_MCP')
# generate x、y and mu --------------------------------------------------------------------
P_y <- function(beta,mu,x) {
  x <- matrix(x)
  y_p <- exp(t(beta)%*%x+mu)/(1+exp(t(beta)%*%x+mu))
  y <- rbinom(1,1,prob=as.numeric(y_p))
  return (y)
}

generate_xymu <- function(mu, beta, n){
  Sigma <- diag(length(beta))#协方差矩阵
  x_data <- mvrnorm(n=n, rep(0, length(beta)), Sigma)#生成多元正态分布的variable
  if (n%%2==0) {
    pos_size <- n%/%2
  } else {
    pos_size <- n%/%2+rbinom(1,1,0.5)
  }                                          #确定不同截距项的数量
  neg_size <- n-pos_size
  mu_vec <- c(rep(mu, pos_size), rep(-mu, neg_size))#生成截距项
  y <- 0
  for(i in 1:n) {
    y[i] <- P_y(beta=beta, mu=mu_vec[i],x=x_data[i,])
  }                                         #生成y
  mydata <- list(xydata = cbind(y,x_data), mu = mu_vec)                 #包含x和y
  return(mydata)
}


# 由e_i和e_j定义的n(n-1)/2✖n矩阵 -------------------------------------------------


generate_Del <- function(n){
  Delta <- matrix(NA, n, 1)
  for (i in (n-1):1) {
    sub <- diag(-1, i, i)
    sub_1 <- matrix(1, 1, i)
    sub_0 <- matrix(0, (n-(i+1)), i)
    mat <- rbind(sub_0, sub_1, sub)
    Delta <- cbind(Delta, mat)
  }
  Delta <- Delta[,-1]
  return(Delta)
}

pdmat <- function(n) {
  .Call('_cfzip_pdmat', PACKAGE = 'cfzip', n)
}

# 定义更新w和b的目标函数  -----------------------------------------------------------


#######################################################################################
Obj_mb<- function(par, data, n_beta, delta, lambda, rho, Delta) {
  # if (n_beta+n_mu!=length(par)||nrow(data)!=n_mu){
  #   stop('incorrect input dimension')
  # }
  beta <- par[1:n_beta]
  mu <- par[(n_beta+1):length(par)]
  x <- data[,2:(n_beta+1)]
  y <- data[,1]
  llk <- (t(beta)%*%t(x)+t(mu))%*%y+sum(-log(1+exp(t(beta)%*%t(x)+t(mu))))
  Obj_val <- 0.5*rho*sum((Delta%*%mu-delta+lambda/rho)^2)-llk
  return(Obj_val)
}

# Update beta
update_beta <- function(beta0 = NULL, mu, data) {
  y <- data[,1]
  X <- data[,2:ncol(data)]
  if (is.null(beta0)) beta0 <- rep(0, ncol(data) - 1)
  update_beta_objfun <- function(beta) {
    -sum(y * X %*% beta) + sum(log(1 + exp(X %*% beta + mu)))
  }
  beta_new <- optim(beta0, update_beta_objfun, method = "BFGS")
  beta_new$par
}

# update mu
update_mu <- function(mu0, beta, data, delta, lambda, rho, tol = 1e-5, maxit = 1000) {
  y <- data[,1]
  X <- data[,2:ncol(data)]
  F <- function(mu, y, X, beta, delta, lambda) {
    n <- length(mu)
    delta1 <- lambda1 <- matrix(0, n, n)
    delta1[lower.tri(delta1, diag=FALSE)] <- delta
    delta1 <- t(delta1)
    lambda1[lower.tri(lambda1, diag=FALSE)] <- lambda
    lambda1 <- t(lambda1)
    Fval <- -(y - exp(X %*% beta + mu) / (1 + exp(X %*% beta + mu))) / n + rho * (n * mu - sum(mu) +
                                                                                    colSums(delta1 - lambda1 / rho) - rowSums(delta1 - lambda1 / rho))
    Fval
  }
  
  gF <- function(mu, X, beta, rho) {
    n <- length(mu)
    diag(exp(X %*% beta + mu) / (1 + exp(X %*% beta + mu))^2 / n) + n * rho * (diag(n) - matrix(1/n,n,n))
  }
  
  nitr <- 0
  mu_new <- mu0
  while (abs(F(mu_new, y, X, beta, delta, lambda)) > tol && nitr <= maxit) {
    mu_old <- mu_new
    gradF <- gF(mu_old, X, beta, rho)
    Fval <- F(mu_old, y, X, beta, delta, lambda)
    lugF <- lu(gradF)
    elu <- expand(lugF)
    L <- elu$L
    U <- elu$U
    P <- elu$P
    a <- t(P) %*% Fval
    b <- forwardsolve(L,a)
    d <- backsolve(U,b)
    mu_new <- mu_old - d
    nitr <- nitr + 1
  }
  mu_new
}

update_mb <- function(par, data, n_beta, delta, lambda, rho, Delta, tol = 1e-3, maxit = 1000) {
  x <- data[,2:(n_beta+1)]
  y <- data[,1]
  par_old <- par
  beta_old <- par_old[1:n_beta]
  mu_old <- par_old[(n_beta+1):length(par)]
  beta_new <- update_beta(beta0 = beta_old, mu = mu_old, data = data)
  mu_new <- update_mu(mu0 = mu_old, beta = beta_new, data = data, delta = delta, lambda = lambda, rho = rho)
  par_new <- c(beta_new, mu_new)
  obj_diff <- abs(Obj_mb(par_new,data = data,n_beta = n_beta,delta = delta,lambda = lambda, rho = rho, Delta = Delta) - Obj_mb(par_old,data = data,n_beta = n_beta,delta = delta,lambda = lambda, rho = rho, Delta = Delta))
  nitr <- 0
  while (obj_diff > tol &&nitr <= maxit) {
    par_old <- par_new
    beta_old <- par_old[1:n_beta]
    mu_old <- par_old[(n_beta+1):length(par)]
    beta_new <- update_beta(beta0 = beta_old, mu = mu_old, data = data)
    mu_new <- update_mu(mu0 = mu_old, beta = beta_new, data = data, delta = delta, lambda = lambda, rho = rho)
    par_new <- c(beta_new, mu_new)
    obj_diff <- abs(Obj_mb(par_new,data = data,n_beta = n_beta,delta = delta,lambda = lambda, rho = rho, Delta = Delta) - Obj_mb(par_old,data = data,n_beta = n_beta,delta = delta,lambda = lambda, rho = rho, Delta = Delta))
    nitr <- nitr + 1
  }
  par_new
}




# 定义soft-thresholding operator --------------------------------------------


ST<- function(t, tau) {
  if (abs(t)>tau) {
    st <- sign(t)*(abs(t)-tau)
  }else {
    st <- 0
  }
  return (st)
}



# 定义SCAD惩罚下的δ更新 -----------------------------------------------------------


update_delta_SCAD <- function(u, rho, tau, kappa = 3.7) {
  if (abs(u)<=tau*(1+1/rho)) {
    delta <- ST(u,tau/rho)
  }else if (abs(u)<=kappa*tau) {
    delta <- ST(u,kappa*tau/((kappa-1)*rho))/(1-1/((kappa-1)*rho))
  }else {
    delta <- u
  }
  return (delta)
}


# 定义MCP惩罚下的δ更新 ------------------------------------------------------------


update_delta_MCP <- function (u, rho, tau, kappa = 3.7) {
  if (abs(u)<=kappa*tau) {
    delta <- ST(u,tau/rho)/(1-1/(kappa*rho))
  }else {
    delta <- u
  }
  return (delta)
}


# define BIC --------------------------------------------------------------


bic <- function(para, data, rnd, n_beta, n_mu, d = 1) {
  y <- data[,1]
  x <- data[,2:(n_beta+1)]
  beta <- para[1:n_beta]
  mu <- para[(n_beta+1):(n_beta+n_mu)]
  n <- length(y)
  Khat <- length(unique(round(mu, rnd)))
  p <- n_beta + 1
  llk <- (t(beta)%*%t(x)+t(mu))%*%y+sum(-log(1+exp(t(beta)%*%t(x)+t(mu))))
  -llk + d * log(log(n + p)) * log(n) / n * (Khat + p)
}



# find min BIC ------------------------------------------------------------


find_bic <- function(xydata, result, id, rnd, n_beta, n, n_mu, d = 1) {
  mu.all <- result[(n_beta+1):(n_beta+n),][,-1]
  beta.all <- result[1:n_beta,][,-1]
  para.all <- result[,-1]
  ntau <- ncol(mu.all)
  bic.all <- rep(NA, ntau)
  id.all <- id
  for (i in 1:ntau) {
    bic_num <- bic(para = para.all[,i], data = xydata, rnd = 2, n_beta = n_beta, n_mu = n_mu, d = d)
    bic.all[i] <- bic_num
  }
  bestid <- which.min(bic.all)
  mu.best <- mu.all[, bestid]
  beta.best <- beta.all[, bestid]
  id.best <- id.all[,bestid]
  result_fin <- list(bic.all = bic.all, bestid = bestid, mu.best = mu.best, beta.best = beta.best,
                     id.best = id.best, ngrp.best = length(unique(round(mu.best, rnd))))
  result_fin
}


# prepare for solution paths ----------------------------------------------


id_data <- function(mu, rnd){
  n <- length(mu)
  uniqcoef <- unique(round(mu, rnd))
  id <- rep(NA, n)
  for (i in 1:n) {
    for (j in 1: length(uniqcoef)){
      if (round(mu[i],rnd) == uniqcoef[j]) id[i] <- j
    }
  }
  return(id)
}

plot_sp <- function(id_all, tau_to, tau_by, mu_set){
  id_all <- cbind(1:nrow(id_all),id_all)
  id_new <- id_all[do.call(order,as.data.frame(id_all[,ncol(id_all):1])),]
  id_origin <- id_new[,1]
  n1 <- nrow(id_all)
  n2 <- ncol(id_all)
  id_plot <- id_new
  Jaccard_vec <- rep(NA, ncol(id_new)-1)
  tau_vec <- rep(NA, ncol(id_new)-1)
  for (i in 1:(ncol(id_new)-1)){
    Jaccard_vec[i] <- gset_similarity(gset(mu_set), gset(as.numeric(id_new[,i+1])))
    tau_vec[i] <- paste0('tau:',tau_by*(i-1))
  }
  Jaccard <- rbind(tau_vec,Jaccard_vec)
  save(Jaccard, file = 'Jaccard.RData')
  for (j in n2:1){
    npts <- length(unique(id_new[,j]))
    from <- n1*npts/(npts+1)
    to <- n1/(npts+1)
    id_relab <- id_new[,j]
    id_relab <- rep(seq_along(unique(id_relab)),times = c(as.data.frame(table(id_relab))[unique(id_relab),2]))
    id_plot[,j] <- rep(seq(from,to,length.out = npts), times = c(as.data.frame(table(id_relab))[,2]))
  }
  save(id_new, file = 'id_new.RData')
  matplot(y = t(id_plot), x = seq(0, tau_to+tau_by, tau_by), type = "l", ylab = "ID", xlab = expression(tau))
}


# Find initial parameter -------------------------------------------------

ini_para <- function(mu, beta, n, noise){
  mu_0 <- mu+rnorm(n = n, mean = noise$mean, sd = noise$sd)#生成初始b
  beta_0 <- beta+rnorm(n = length(beta), mean=noise$mean, sd = noise$sd)#生成初始w
  # loglikelihood_0 <- function(par,b_real,n,data) {
  #   w <- par[1:2]
  #   b <- b_real
  #   x <- data[,2:(n+1)]
  #   y <- data[,1]
  #   llk <- 0
  #   llk <- (t(w)%*%t(x)+t(b))%*%y+sum(-log(1+exp(t(w)%*%t(x)+t(b))))
  #   return (-llk)                                           #定义用于确定初值的Loglikelihood
  # }
  # w_0 <- result$par#生成初始w
  #b_0 <- result$par[(length(w)+1)]#生成初始b
  delta <- 0
  vec_n <- 0
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      vec_n <- vec_n+1
      # vec_n <- j-i+sum(seq((n-1), by = -1, length.out = (i-1)))#确定元素在向量中的位置
      delta[vec_n] <- mu_0[i]-mu_0[j]                                    #生成初始delta
    }
  }
  lambda <- rep(0, n*(n-1)/2)#生成初始lambda
  result <- data.frame(ID=c(1:(n+length(beta))))#生成初始result
  id <- as.data.frame(seq(1, n, 1))#生成初始的id
  para <- list(mu = mu_0, beta = beta_0, delta = delta, lambda = lambda, id = id)
}

# core function ------------------------------------------------------
core_fun <- function(tau_from, tau_to, tau_by, epi, beta, mu, delta, lambda, update_del, Delta, mydata){
  id <- data.frame(ini = c(1:length(mu)))
  result <- data.frame(ini = c(1:(length(mu)+length(beta))))
  para <- c(beta, mu)  #设定迭代的初值
  for (tau in seq(tau_from, tau_to, tau_by)) {
    residual_norm <- 1#设置初始残差
    epi <- epi#设置精度
    n_iter <- 0#迭代次数
    while(residual_norm>epi & n_iter<=1000) {
      n_iter <- n_iter+1
      # delta <- Del
      result_update_bm <- update_mb(par = para, data = mydata$xydata, delta = delta, lambda = lambda, rho = 3, Delta = Delta, n_beta = length(beta))
      para <- result_update_bm#update of par
      mu_update <- para[(length(beta)+1):(length(beta)+n)]#更新后的b
      u <- Delta%*%mu_update+lambda/rho
      delta <- sapply(u, FUN = update_del, tau = tau, rho = rho)#更新delta
      lambda <- lambda+rho*(Delta%*%mu_update-delta)#更新lambda
      residual_norm <- sqrt(sum((Delta%*%mu_update-delta)^2))
      words <- paste('tau:',tau,',第', n_iter,'次迭代，残差为', residual_norm)
      print(words)
    }
    col_name <- paste0('tau:',tau)
    id[col_name]<- id_data(mu = mu_update, rnd = 2)
    result[col_name] <- para
  }
  result_rec <- list(para = result, id_all = id)
}



# simulation --------------------------------------------------------------


simulation <- function(times, mu, beta, n, update_method, tau_to, tau_by){
  for(i in 1:times){
    mu_set <- gset(c(rep(1,n/2),rep(2,n/2)))
    mydata <- generate_xymu(mu = mu, beta = beta, n = n)
    Del <- pdmat(n)
    # Del <- generate_Del(n)
    ini_set <- ini_para(mu = mydata$mu, beta = beta, n = n, noise = list(mean = 0, sd = 0.2*mu))
    result_record <- core_fun(tau_from = 0, tau_to = tau_to, tau_by = tau_by, 
                              epi = 0.001, beta = ini_set$beta, mu = ini_set$mu, 
                              delta = ini_set$delta, lambda = ini_set$lambda, 
                              update_del = update_method, Delta = Del, mydata = mydata)
    
    id_plot <- result_record$id_all[,-1]
    filename_img <- paste("",i,"th Solution paths.jpg")
    jpeg(file = filename_img)#保存图片
    plot_sp(id_plot, tau_to, tau_by, mu_set)  ## 画图程序
    dev.off();
    
    filename0 <- paste('result_record_',i,'RData')
    filename1 <- paste('result_para_',i,'.RData')
    filename2 <- paste('mydata',i,'.RData')
    result_para <- result_record$para
    save(result_record, file = filename0)
    save(result_para, file = filename1)
    save(mydata, file = filename2)
  }
}

find_tau <- function(times, c){
  Khat <- rep(NA, times)
  RMSE_mu <- rep(NA, times)
  RMSE_beta <- rep(NA, times)
  tab <- data.frame(Khat, RMSE_mu, RMSE_beta)
  K_mean <- NA
  K_median <- NA
  K_se <- NA
  RMSEmu_se <- NA
  RMSEmu_mean <- NA
  RMSEbeta_se <- NA
  RMSEbeta_mean <- NA
  statistics <- data.frame(K_mean, K_median, K_se, RMSEmu_mean, RMSEmu_se, RMSEbeta_mean, RMSEbeta_se)
  bic_tab <- matrix(NA, nrow = times, ncol = ceiling(tau_to/tau_by)+1)
  for(i in 1:times){
    filename0 <- paste('result_record_',i,'.RData')
    filename1 <- paste('result_para_',i,'.RData')
    filename2 <- paste('mydata',i,'.RData')
    load(filename0)
    load(filename1)
    load(filename2)
    result_fin <- find_bic(xydata = mydata$xydata, result = result_record$para,
                           id = result_record$id_all, n = n,
                           n_beta = length(beta), n_mu = n,rnd = 2 , d = c)
    bic_tab[i,] <- result_fin$bic.all 
    best_tau <- names(result_record$para)[result_fin$bestid+1]
    tab$tau[i] <- best_tau
    tab$Khat[i] <- length(unique(round(result_fin$mu.best, 2)))
    tab$RMSE_mu[i] <- (sum((result_fin$mu.best-mydata$mu)^2)/length(mydata$mu))^0.5
    tab$RMSE_beta[i] <- (sum((result_fin$beta.best-beta)^2)/length(beta))^0.5
  }
  statistics$K_mean <- mean(tab$Khat)
  statistics$K_median <- median(tab$Khat)
  statistics$K_se <- sd(tab$Khat)
  statistics$RMSEmu_mean <- mean(tab$RMSE_mu)
  statistics$RMSEmu_se <- sd(tab$RMSE_mu)
  statistics$RMSEbeta_mean <- mean(tab$RMSE_beta)
  statistics$RMSEbeta_se <- sd(tab$RMSE_beta)
  
  bic_tab <- as.data.frame(bic_tab)
  save(tab, file = 'index.RData')
  save(bic_tab, file = 'bic_Table.RData')
  save(statistics, file = 'statistics.RData')
}

# 变量设置 --------------------------------------------------------------------
mu <- 1#设定截距项
beta <- as.matrix(c(-0.3, 1, -1, 2, 0.5))
kappa <- 3.7 
rho <- 3
n <- 500
e <- rep(0,n)

tau_to <- 0.65
tau_by <- 0.025
iterations_num <- 20

simulation(times = iterations_num, mu = mu, beta = beta, n = n, update_method = update_delta_MCP, 
           tau_to = tau_to, tau_by = tau_by)
# bic_c <-1
# find_tau(times = iterations_num, c = bic_c)
