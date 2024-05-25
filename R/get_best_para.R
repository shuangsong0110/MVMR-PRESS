#' @title Tune parameters
#' @description A unified Bayesian framework to calculate PRS based on GWAS summary statistics
#'
#' @param summs_exposure GWAS summary statistics for exposure, list, each element should be a data.frame including SNP, A1, A2, BETA, SE
#' @param summs_outcome GWAS summary statistics for outcome, a data.frame including SNP, A1, A2, BETA, SE
#' @param P_hat correlation matrix of exposures
#' @param snp.use rsid of the IVs
#' @param n_cores Number of cores running in parallel
#' @param n_iter_cv Iterations in CV
#' @return best_parameters
#' @import parallel
#'
#' @export
#'
get_best_para <- function(summs_exposure, summs_outcome, P_hat, snp.use, ncv=3,ncores=100,n_iter_cv=30){
  set.seed(0)
  # snp.use <- fread(paste0('/home/songs/mvmrbayes/0406new/share_snp_clumping_format/',
  #                         trait1,'_adjust_0.005.clumped'))$SNP
  ## strict: --clump-kb 10000 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001
  ## 2:  --clump-kb 1000 --clump-p1 5e-8 --clump-p2 0.01 --clump-r2 0.1
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
  GC=1:p
  folds=rep(1,p)
  rou=round(p/ncv,0)
  for(i in 2:ncv){
    sam=sample(which(folds==1),rou)
    folds[sam]=i
  }
  likemax=-100000
  lambda_c=NULL
  mu_c=NULL
  #cho=c(100,10,0.1,0.01)
  cho=c(10,3,1,0.3,0.1)
  lambda_all <- mu_all <- ifold_all <- c()
  for(lambda in cho){
    for (mu in cho){
      for(ifold in 1:ncv){
        lambda_all <- c(lambda_all,lambda)
        mu_all <- c(mu_all,mu)
        ifold_all <- c(ifold_all,ifold)
      }
    }}
  para <- data.frame(lambda=lambda_all,mu=mu_all,ifold=ifold_all)
  #library(parallel)
  res_para <- mclapply(1:nrow(para),get.cv.res.mc,bx=bx,bxse=bxse,
                       by=by,byse=byse,P_hat=P_hat,para=para,folds=folds,GC=GC,
                       iter=n_iter_cv,
                       mc.cores=ncores)
  ##find the largest llh parameter settings
  llh <- lambda_all <- mu_all <- c()
  temp <- unlist(lapply(res_para,'[','lih'))
  for(lambda in cho){
    for (mu in cho){
      #for(ifold in 1:5){
      ind <- which(para$lambda==lambda & para$mu==mu)
      lambda_all <- c(lambda_all,lambda)
      mu_all <- c(mu_all,mu)
      llh <- c(llh,mean((as.numeric(temp[ind]))))
      # }
    }
  }
  tt <- which.max(llh)
  cbind(lambda_all,mu_all,llh)
  return(list(lambda=lambda_all[tt],
              mu=mu_all[tt]))

}


get.cv.res.mc <- function(bx,bxse,by,byse,P_hat,para,k,folds, GC=1:p,iter=20){
  #print(k)
  lambda <- para$lambda[k]
  mu <- para$mu[k]
  ifold <- para$ifold[k]
  return(get.cv.res(bx,bxse,by,byse,P_hat,lambda=lambda,
                    mu=mu,
                    ifold=ifold,
                    folds=folds,
                    GC=GC,iter=iter))
}

get.cv.res <- function(bx,bxse,by,byse,P_hat,lambda,mu,ifold, folds, GC=1:p,iter=20){
  print(paste("lambda",lambda))
  print(paste("mu",mu))
  lih=NULL
  te=GC[folds==ifold]#test data
  tr=GC[folds!=ifold]#training data
  #r=cvre(lambda,mu,tr,30)#calculate theta and gamma
  #bx,bxse,by,byse,P_hat, lambda_theta,mu_gamma,G,iter=20,sigma_B=NULL
  r=cvre(bx,bxse,by,byse,P_hat,
         lambda_theta=lambda,mu_gamma=mu,
         G=tr,iter=iter,sigma_B=NULL)
  #print('haha')
  theta_r=r[[1]]
  gamma_r=r[[2]]
  lih=lihre(bx,bxse,by,byse,P_hat, lambda=lambda,mu=mu,theta=theta_r,gamma=gamma_r,
            G=te, sigma_B=NULL)
  lih ## log likelihood
  #lih=lihre(lambda,mu,theta_r,gamma_r,te)#calculate likelihood
  #print(paste("theta_r",theta_r))
  return(list(lih=as.numeric(lih),theta=theta_r))
}

