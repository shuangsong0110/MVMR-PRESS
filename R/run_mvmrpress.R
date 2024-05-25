#' @title Main function
#' @description MVMRPRESS
#'
#' @param summs_exposure GWAS summary statistics for exposure, list, each element should be a data.frame including SNP, A1, A2, BETA, SE
#' @param summs_outcome GWAS summary statistics for outcome, a data.frame including SNP, A1, A2, BETA, SE
#' @param P_hat correlation matrix of exposures
#' @param snp_use rsid of the IVs
#' @param lambda_theta Parameter
#' @param mu_gamma Parameter
#' @param n_boots Number of boostraps
#' @param n_cores Number of cores running in parallel
#' @param n_iter Iterations
#' @param n_iter_boot Iterations in bootstrap
#' @return result
#' @import parallel
#' @import Matrix
#'
#' @export
#'
run_mvmrpress <- function(summs_exposure, summs_outcome, snp_use, P_hat, para_theta=0.01, para_gamma=10,
                          n_boots=100, n_cores=100,n_iter=200,n_iter_boot=100){
  #setwd('/home/songs/mvmrbayes/0406new/share_snp_clumping_format')
  # load(paste0('/home/songs/mvmrbayes/0406new/res/',trait2,'_',
  #             trait1,'.adjust1.newP1.',
  #             'lambda',lambda_theta,'_mu',mu_gamma,'.RData'))
  #summs2 <- summs_outcome
  snp.use <- snp_use
  n_exposure <- length(summs_exposure)
  lambda_theta=para_theta
  mu_gamma=para_gamma
  summs2 <- summs_outcome[match(snp.use, summs_outcome$SNP),]
  summs2$A1 <- toupper(summs2$A1)
  summs2$A2 <- toupper(summs2$A2)
  for(k in 1:length(summs_exposure)){
    summary1 <- summs_exposure[[k]]
    summs1 <- summary1[match(snp.use, summary1$SNP),]
    summs1$A1 <- toupper(summs1$A1)
    summs1$A2 <- toupper(summs1$A2)
    sig <- agtc(summs1$A1,summs1$A2,summs2$A1,summs2$A2)
    if(k==1){
      bx <- summs1$BETA*sig
      bxse <- summs1$SE
    }else{
      bx <- cbind(bx,summs1$BETA*sig)
      bxse <- cbind(bxse,summs1$SE)
    }
  }
  by <- summs2$BETA
  byse <- summs2$SE
  p <- nrow(bx)
  k <- ncol(bx)

  theta_s=rep(0,k)
  gamma_s=rep(0,p)
  res <- MVMR_Bayes(bx,bxse,P_hat,by,byse,lambda_theta=lambda_theta,mu_gamma=mu_gamma,
                    sigma_B=NULL,theta=theta_s, gamma=gamma_s,iter=n_iter,return.se=T)
  theta_s=rep(0,k)
  gamma_s=rep(0,p)
  res.boot <- mclapply(1:n_boots,MVMR_Bayes_boots,bx=bx,bxse=bxse,
                       P_hat=P_hat,by=by,byse=byse,lambda_theta=lambda_theta,mu_gamma=mu_gamma,
                       sigma_B=NULL,theta=theta_s, gamma=gamma_s,iter=n_iter_boot,return.se=F,
                       mc.preschedule=F,mc.cores=n_cores)
  boot.mat <- matrix(unlist(res.boot),nrow=n_exposure)
  res_mvmrpress <- list()
  res_mvmrpress$beta_est <- res$theta
  res_mvmrpress$beta_est_SE <- (res$theta/apply(boot.mat,1,sd)/sqrt(n_exposure)*1.96)
  res_mvmrpress$beta_Z <- res_mvmrpress$beta_est/res_mvmrpress$beta_est_SE
  res_mvmrpress$beta_P <- 2*pnorm(-abs(res_mvmrpress$beta_Z))

  return(res_mvmrpress)

}

#' @import Matrix
MVMR_Bayes=function(bx,bxse,P_hat,by,byse,lambda_theta,mu_gamma,
                    sigma_B=NULL,theta,gamma,iter=200,return.se=F){
  #beta_hat=c(B_hat)
  m=nrow(bx)
  k=ncol(bx)
  temp <- 0
  Omega_hat_list <- Omega_hat_inv_list <- list()
  temp <- rep(0,k)
  for(l in 1:m){
    Omega_hat_list[[l]] <- diag(bxse[l,])%*%P_hat%*% diag(bxse[l,])
    Omega_hat_inv_list[[l]] <- diag(1/bxse[l,])%*%solve(P_hat)%*% diag(1/bxse[l,])
    for(j in 1:k){
      temp[j] <- temp[j]+Omega_hat_list[[l]][j,j]
    }
    #temp <- temp+sum(diag(Omega_hat_list[[l]]))
  }
  if(is.null(sigma_B)){ ## variance
    sigma_B <- c()
    for(l in 1:k){
      sigma_B[l]=(sum(diag(c(bx[,l])%*%t(c(bx[,l]))))-temp[l])/(m)
    }
    #sigma_B=(c(bx1)%*%t(c(bx1))-Omega_hat_x) #(m*k) * (m*k)
  }
  Sigma_inv=diag(1/byse^2) ## diagonal m*m
  #Omega_inv=solve(Omega_hat) ## 1, m+1, 2m+1, ... correlated
  order.mat <- matrix(1:(m*k),m,k)
  order1 <- c(t(order.mat))
  order2 <- order(order1)
  set.seed(0)
  i <- 1
  diff <- 100
  theta0 <- rep(100,k)
  while(diff>5e-4 & i<=iter){
    cat(paste0('Iteration ',i,' / ',iter,' ...','\n'))
    print(theta)
    #for(i in 1:iter){
    #print(i)
    theta.mat <- theta%*%t(theta)
    Sigma_beta_list <- list()
    mu_beta_part1_list <- list() ##omega_inv%*%beta
    for(l in 1:m){
      Sigma_beta_list[[l]] <- solve(Omega_hat_inv_list[[l]]+
                                      theta.mat/byse[l]^2+diag(1/sigma_B))

      mu_beta_part1_list[[l]] <- Omega_hat_inv_list[[l]]%*%bx[l,]
    }
    # Sigma_beta=solve(Omega_inv+kronecker((theta%*%t(theta)),Sigma_inv)
    #                  +1/(sigma_B)*diag(1,(m*k)))
    Sigma_beta0 <- bdiag(Sigma_beta_list)
    Sigma_beta <- Sigma_beta0[order2,order2]
    mu_beta_part1 <- unlist(mu_beta_part1_list)[order2]
    mu_beta_part2 <- kronecker(theta, Sigma_inv%*%(by-gamma))
    mu_beta=Sigma_beta%*%(mu_beta_part1+mu_beta_part2)
    #E step
    mmu_beta=matrix(mu_beta,ncol=k)
    b_theta=t(mmu_beta)%*%Sigma_inv%*%(by-gamma)
    A_theta=matrix(rep(0,k*k),ncol=k)
    system.time(
      for(j in 1:k){
        for(l in 1:k){
          A_theta[j,l]=sum(diag(Sigma_inv%*%Sigma_beta[(m*(j-1)+1):(m*j),(m*(l-1)+1):(m*l)]))
        }
      }
    )
    A_theta=A_theta+t(mmu_beta)%*%Sigma_inv%*%(mmu_beta)
    A_gamma=Sigma_inv
    b_gamma=Sigma_inv%*%(by-mmu_beta%*%theta)

    #M step
   # system.time(
      theta <- PGD(A_theta,b_theta,theta,lambda_theta)
   # )
    # theta
    #print(theta)
    #system.time(
      gamma <- PGD(A_gamma,b_gamma,gamma,mu_gamma)
    #)
    diff <- sum(abs(theta0-theta))
    theta0 <- theta
    # system.time(
    #   gamma1 <- PGD(A_gamma,b_gamma,gamma,mu_gamma)
    # )
    # print(summary(gamma-gamma_true))
    i <- i+1
  }
  if(return.se){
    VAR <- Louis_t(A_theta, mu_beta, Sigma_beta_list,Sigma_beta, order1, m,k, by,byse,
                   theta,gamma)
    r=list(theta=theta,gamma=gamma,VAR=VAR,diff=diff,iter=i)
    # (B_hat,Omega_hat,beta_hat_y,Sigma_hat,lambda_theta,mu_gamma,
    #   sigma_B,theta_r[[1]],theta_r[[2]])
  }else{
    r=list(theta=theta,gamma=gamma,diff=diff,iter=i)
  }
  return(r)

}

#' @import Matrix
MVMR_Bayes_boots=function(bx,bxse,P_hat,by,byse,lambda_theta,mu_gamma,
                          sigma_B=NULL,theta,gamma,iter=50,return.se=F,seed=0){
  set.seed(seed)
  #beta_hat=c(B_hat)
  m=nrow(bx)
  k=ncol(bx)
  boot.ind <- sample(m,m,replace = T)
  bx <- bx[boot.ind,]
  by <- by[boot.ind]
  bxse <- bxse[boot.ind,]
  byse <- byse[boot.ind]
  temp <- 0
  Omega_hat_list <- Omega_hat_inv_list <- list()
  temp <- rep(0,k)
  for(l in 1:m){
    Omega_hat_list[[l]] <- diag(bxse[l,])%*%P_hat%*% diag(bxse[l,])
    Omega_hat_inv_list[[l]] <- diag(1/bxse[l,])%*%solve(P_hat)%*% diag(1/bxse[l,])
    for(j in 1:k){
      temp[j] <- temp[j]+Omega_hat_list[[l]][j,j]
    }
    #temp <- temp+sum(diag(Omega_hat_list[[l]]))
  }
  if(is.null(sigma_B)){ ## variance
    sigma_B <- c()
    for(l in 1:k){
      sigma_B[l]=(sum(diag(c(bx[,l])%*%t(c(bx[,l]))))-temp[l])/(m)
    }
    #sigma_B=(c(bx1)%*%t(c(bx1))-Omega_hat_x) #(m*k) * (m*k)
    ##?
  }
  #
  # for(l in 1:m){
  #   Omega_hat_list[[l]] <- diag(bxse[l,])%*%P_hat%*% diag(bxse[l,])
  #   Omega_hat_inv_list[[l]] <- diag(1/bxse[l,])%*%solve(P_hat)%*% diag(1/bxse[l,])
  #   temp <- temp+sum(diag(Omega_hat_list[[l]]))
  # }
  # if(is.null(sigma_B)){
  #   #sigma_B=(c(bx1)%*%t(c(bx1))-Omega_hat_x) #(m*k) * (m*k)
  #   sigma_B=(sum(diag(c(bx)%*%t(c(bx))))-temp)/(m*k) ##?
  # }
  #
  Sigma_inv=diag(1/byse^2) ## diagonal m*m
  #Omega_inv=solve(Omega_hat) ## 1, m+1, 2m+1, ... correlated
  order.mat <- matrix(1:(m*k),m,k)
  order1 <- c(t(order.mat))
  order2 <- order(order1)
  set.seed(0)
  i <- 1
  diff <- 100
  theta0 <- rep(100,k)
  while(diff>1e-3 & i<=iter){
    #for(i in 1:iter){
    #print(i)
    theta.mat <- theta%*%t(theta)
    Sigma_beta_list <- list()
    mu_beta_part1_list <- list() ##omega_inv%*%beta
    for(l in 1:m){
      Sigma_beta_list[[l]] <- solve(Omega_hat_inv_list[[l]]+
                                      theta.mat/byse[l]^2+diag(1/sigma_B))

      mu_beta_part1_list[[l]] <- Omega_hat_inv_list[[l]]%*%bx[l,]
    }
    # Sigma_beta=solve(Omega_inv+kronecker((theta%*%t(theta)),Sigma_inv)
    #                  +1/(sigma_B)*diag(1,(m*k)))
    Sigma_beta0 <- bdiag(Sigma_beta_list)
    Sigma_beta <- Sigma_beta0[order2,order2]
    mu_beta_part1 <- unlist(mu_beta_part1_list)[order2]
    mu_beta_part2 <- kronecker(theta, Sigma_inv%*%(by-gamma))
    mu_beta=Sigma_beta%*%(mu_beta_part1+mu_beta_part2)
    #E step
    mmu_beta=matrix(mu_beta,ncol=k)
    b_theta=t(mmu_beta)%*%Sigma_inv%*%(by-gamma)
    A_theta=matrix(rep(0,k*k),ncol=k)
    system.time(
      for(j in 1:k){
        for(l in 1:k){
          A_theta[j,l]=sum(diag(Sigma_inv%*%Sigma_beta[(m*(j-1)+1):(m*j),(m*(l-1)+1):(m*l)]))
        }
      }
    )
    A_theta=A_theta+t(mmu_beta)%*%Sigma_inv%*%(mmu_beta)
    A_gamma=Sigma_inv
    b_gamma=Sigma_inv%*%(by-mmu_beta%*%theta)

    #M step
    #system.time(
      theta <- PGD(A_theta,b_theta,theta,lambda_theta)
   # )
    # theta
    #print(theta)
    #system.time(
      gamma <- PGD(A_gamma,b_gamma,gamma,mu_gamma)
   # )
    diff <- sum(abs(theta0-theta))
    theta0 <- theta
    # system.time(
    #   gamma1 <- PGD(A_gamma,b_gamma,gamma,mu_gamma)
    # )
    # print(summary(gamma-gamma_true))
    i <- i+1
  }
  #print(theta)
  if(return.se){
    SE <- Louis_t(A_theta, mu_beta, Sigma_beta_list,Sigma_beta, order1, m,k, by,byse,
                  theta,gamma)
    r=list(theta=theta,gamma=gamma,SE=SE,diff=diff,iter=i)
    # (B_hat,Omega_hat,beta_hat_y,Sigma_hat,lambda_theta,mu_gamma,
    #   sigma_B,theta_r[[1]],theta_r[[2]])
  }else{
    r=list(theta=theta)
  }
  return(r)

}


#' @import Matrix
Louis_t=function(A_theta, mu_beta, Sigma_beta_list, Sigma_beta,order1, m,k, by,byse,
                 theta,gamma){
  # Sigma_beta_list: by column
  ix=A_theta
  iz=matrix(0,k,k)
  Sigma_inv=diag(1/byse^2) ## diagonal m*m
  for(j in 1:k){
    for(l in j:k){
      A1=matrix(0,k,k)
      A2=matrix(0,k,k)
      A1[,j]=theta
      A1[j,]=A1[j,]+t(theta)
      A2[,l]=theta # k*k
      A2[l,]=A2[l,]+t(theta)
      temp0 <- 0
      #temp1 <- 0
      for(q in 1:m){
        temp0 <- temp0+sum(diag(A1%*%Sigma_beta_list[[q]]%*%A2%*%Sigma_beta_list[[q]]
        ))/byse[q]^4
      }
      A1=kronecker(A1, Sigma_inv)
      A2=kronecker(A2, Sigma_inv)
      iz[j,l]=temp0/2+as.numeric(t(mu_beta)%*%A1%*%Sigma_beta%*%A2%*%mu_beta)
      t1=rep(0,k)
      t2=rep(0,k)
      t1[j]=1
      t2[l]=1
      iz[j,l]=iz[j,l]+as.numeric(kronecker(t(t2), t(by-gamma)%*%Sigma_inv)%*%Sigma_beta%*%
                                   kronecker(t1, Sigma_inv%*%(by-gamma)))
      iz[j,l]=iz[j,l]-as.numeric(kronecker(t(t2), t(by-gamma)%*%Sigma_inv)%*%
                                   Sigma_beta%*%A1%*%mu_beta-
                                   kronecker(t(t1), t(by-gamma)%*%Sigma_inv)%*%Sigma_beta%*%A2%*%mu_beta)
      iz[l,j]=iz[j,l]
    }
  }
  iy=ix-iz
  vartheta=solve(iy)
  return(vartheta)
}

#' @import Matrix
lihre=function(bx,bxse,by,byse,P_hat, lambda,mu,theta,gamma,G, sigma_B=NULL){
  #calculate the basic item to calculate the likelihood
  m=nrow(bx)
  k=ncol(bx)
  order.mat <- matrix(1:(m*k),m,k)
  order1 <- c(t(order.mat[G,]))
  order2 <- order(order1)
  m <- length(G)
  #print(length(G))
  bx1=bx[G,]
  by1=by[G]
  bxse1 <- bxse[G,]
  byse1 <- byse[G]
  Sx_x=diag(c(bxse1))
  Sigma_hat_y=diag(byse1^2)
  Sigma_inv=diag(1/byse1^2) ## diagonal m*m

  Omega_hat_list <- Omega_hat_inv_list <- list()


  temp <- rep(0,k)
  for(l in 1:m){
    Omega_hat_list[[l]] <- diag(bxse[l,])%*%P_hat%*% diag(bxse[l,])
    Omega_hat_inv_list[[l]] <- diag(1/bxse[l,])%*%solve(P_hat)%*% diag(1/bxse[l,])
    for(j in 1:k){
      temp[j] <- temp[j]+Omega_hat_list[[l]][j,j]
    }
    #temp <- temp+sum(diag(Omega_hat_list[[l]]))
  }
  if(is.null(sigma_B)){ ## variance
    sigma_B <- c()
    for(l in 1:k){
      sigma_B[l]=(sum(diag(c(bx[,l])%*%t(c(bx[,l]))))-temp[l])/(m)
    }
    #sigma_B=(c(bx1)%*%t(c(bx1))-Omega_hat_x) #(m*k) * (m*k)
    ##?
  }
  theta.mat <- theta%*%t(theta)
  Sigma_beta_list <- list()
  mu_beta_part1_list <- list() ##omega_inv%*%beta
  for(l in 1:length(G)){
    Sigma_beta_list[[l]] <- solve(Omega_hat_inv_list[[l]]+
                                    theta.mat/byse1[l]^2+diag(1/sigma_B))

    mu_beta_part1_list[[l]] <- Omega_hat_inv_list[[l]]%*%bx1[l,]
  }
  gamma=rep(0,length(G))
  # Sigma_beta=solve(Omega_inv+kronecker((theta%*%t(theta)),Sigma_inv)
  #                  +1/(sigma_B)*diag(1,(m*k)))
  Sigma_beta0 <- bdiag(Sigma_beta_list)
  Omega_hat_inv0 <- bdiag(Omega_hat_inv_list)
  Omega_hat_inv <- Omega_hat_inv0[order2,order2]
  Sigma_beta <- Sigma_beta0[order2,order2] ## by column
  mu_beta_part1 <- unlist(mu_beta_part1_list)[order2]
  mu_beta_part2 <- kronecker(theta, Sigma_inv%*%(by1-gamma))
  mu_beta=Sigma_beta%*%(mu_beta_part1+mu_beta_part2)
  #E step
  mmu_beta=matrix(mu_beta,ncol=k)
  if(T){
    b_theta=t(mmu_beta)%*%Sigma_inv%*%(by1-gamma)
    A_theta=matrix(rep(0,k*k),ncol=k)
    system.time(
      for(j in 1:k){
        for(l in 1:k){
          A_theta[j,l]=sum(diag(Sigma_inv%*%Sigma_beta[(m*(j-1)+1):(m*j),(m*(l-1)+1):(m*l)]))
        }
      }
    )
    A_theta=A_theta+t(mmu_beta)%*%Sigma_inv%*%(mmu_beta)
    A_gamma=Sigma_inv
  }
  #b_gamma=Sigma_inv%*%(by-mmu_beta%*%theta)

  #claculate likelihood
  # likeli <- -0.5*t(mu_beta-c(bx1))%*%Omega_hat_inv%*%(mu_beta-c(bx1))-
  #   0.5*t(by1-mmu_beta%*%theta)%*%Sigma_inv%*%(by1-mmu_beta%*%theta)
  #  -lambda*sum(abs(theta))
  #-sum(c(bx1)^2)/2/sigma_B
  # likeli=t(by1)%*%Sigma_inv%*%gamma+
  #   t(theta)%*%b_theta-
  #   0.5*t(theta)%*%A_theta%*%theta-
  #  0.5*t(gamma)%*%Sigma_inv%*%gamma
  # -lambda*sum(abs(theta))-mu*sum(abs(gamma))
  likeli=t(theta)%*%b_theta-
    0.5*t(theta)%*%A_theta%*%theta-
    lambda*sum(abs(theta))
  return(likeli)
}


#' @import Matrix
cvre=function(bx,bxse,by,byse,P_hat, lambda_theta,mu_gamma,G,iter=20,sigma_B=NULL){
  #G=t
  m=nrow(bx)
  order.mat <- matrix(1:(m*k),m,k)
  order1 <- c(t(order.mat[G,]))
  order2 <- order(order1)
  m <- length(G)
  k=ncol(bx)
  #print(length(G))
  bx1=bx[G,]
  by1=by[G]
  bxse1 <- bxse[G,]
  byse1 <- byse[G]
  Omega_hat_list <- Omega_hat_inv_list <- list()
  # Sx_x=diag(c(bxse1))
  # Sigma_hat_y=diag(byse1^2)
  # Sigma_inv=diag(1/byse1^2) ## diagonal m*m
  temp <- rep(0,k)
  for(l in 1:m){
    Omega_hat_list[[l]] <- diag(bxse[l,])%*%P_hat%*% diag(bxse[l,])
    Omega_hat_inv_list[[l]] <- diag(1/bxse[l,])%*%solve(P_hat)%*% diag(1/bxse[l,])
    for(j in 1:k){
      temp[j] <- temp[j]+Omega_hat_list[[l]][j,j]
    }
    #temp <- temp+sum(diag(Omega_hat_list[[l]]))
  }
  if(is.null(sigma_B)){ ## variance
    sigma_B <- c()
    for(l in 1:k){
      sigma_B[l]=(sum(diag(c(bx[,l])%*%t(c(bx[,l]))))-temp[l])/(m)
    }
    #sigma_B=(c(bx1)%*%t(c(bx1))-Omega_hat_x) #(m*k) * (m*k)
    ##?
  }

  #print(length(G))
  #Sigma_hat_y=diag(byse1^2)
  #Bx=B_hat[G,]
  #betay=beta_hat_y[G]
  #Sx_x=diag(c(B_SE[G,]))
  # Omega_hat_x=Sx_x%*%(kronecker(P_hat,diag(1,nrow = length(G))))%*%Sx_x
  # sigma_B_x=(c(Bx)%*%t(c(Bx))-Omega_hat_x)
  # sigma_B_x=sum(diag(sigma_B_x))/(length(G)*k)
  #
  theta_s=rep(0,k)
  gamma_s=rep(0,length(G))

  r=MVMR_Bayes(bx1,bxse1,P_hat,by1,byse1,lambda_theta,mu_gamma,
               sigma_B,theta_s,gamma_s,iter=iter)
  # r=MVMR_Bayes(Bx,Omega_hat_x,betay,Sigma_hat_y,lambda_theta,mu_gamma,sigma_B_x,theta_s
  #              ,gamma_s,iter)
  return(r)
}

