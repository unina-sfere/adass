## Data generation ---------------------------------------------------------


simulate_data<-function(case,n_obs=3000) {

  length_tot<-500
  grid_s<-grid_t<-seq(0,1,length.out = length_tot)



  # generate X --------------------------------------------------------------

  domain<-c(0,1)
  n_basis_x<-32     #random chosen between 10 and 50
  X_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_x)
  X_coef<-matrix(rnorm(n_obs*n_basis_x),nrow = n_basis_x,ncol = n_obs )
  X_fd<-fd(X_coef,X_basis)
  X<-eval.fd(grid_s,X_fd)

  # Generate ERROR ----------------------------------------------------------


  n_basis_eps<-20 #random chosen between 10 and 50
  eps_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_eps)
  eps_coef<-matrix(rnorm(n_obs*n_basis_eps),nrow = n_basis_eps,ncol = n_obs )
  eps_fd<-fd(eps_coef,eps_basis)
  Eps<-t(eval.fd(grid_t,eps_fd))


  # Define beta -----------------------------------------------------------

  if(case=="Scenario I"){
    cat("Scenario I")
    beta<-function(s,t){

      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')

      z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
      z_matrix}
      else{
        z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
        z_matrix}
    }



  }
  if(case=="Scenario II"){
    cat("Scenario II")
    beta<-function(s,t){
      a=0.15
      b=0.15
      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')
      z<- -(((d$s-0.5)/a)^2 + ((d$t-0.5)/b)^2) +3
      z[z<0]<-0
      z_matrix<-matrix(z,nrow=length(s))
      z_matrix}
      else{
        z<- -(((s)/a)^2 + ((d)/b)^2) + 1
        z[z<0]<-0
        z+10}
    }}
  if(case=="Scenario III"){
    cat("Scenario III")

    beta<-function(s,t){
      a<-0.05
      b<-0.5
      f_1<-function(s,t){b*(1-s)*sin(10*pi*(s-a-1+sqrt(1-(t-0.5)^2)))}
      f_2<-function(s,t){b*sin(10*pi*(s+a+1-sqrt(1-(t-0.5)^2)))}

      z<-matrix(0,length_tot,length_tot)
      for (ii in 1:length(grid_t)) {
        t<-grid_t[ii]
        s_0_1<-grid_s[grid_s>( a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_0_2<-grid_s[grid_s<( -a+sqrt(1-(t-0.5)^2))&grid_s>=0.5]
        s_n0_1<-grid_s[grid_s<=(a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_n0_2<-grid_s[grid_s>=(-a+sqrt(1-(t-0.5)^2))&grid_s>0.5]
        z_i<-c(f_1(s_n0_1,t),rep(0,length(s_0_1)),rep(0,length(s_0_2)),f_2(s_n0_2,t))
        z[ii,]=z_i
      }
      return(t(z))
    }
  }

  if(case=="Scenario IV"){
    cat("Scenario IV")
    beta<-function(s,t){
      a=0.5
      b=0.5
      c=0.5
      d=0.5
      f_1<-function(s,t){((t-0.5)/c)^3+((s-0.5)/d)^3+((t-0.5)/b)^2 + ((s-0.5)/a)^2+5}
      z<- outer(s,t,f_1)
      z
    }
  }
  if(case=="Historical"){
    cat("Historical")
    beta<-function(grid_s,grid_t){

      f_1<-function(s,t){0.5*sin(30*s)+0.5*cos(30*t)}
      grid<-expand.grid(s=grid_s,t=grid_t)
      z<-f_1(grid$s,grid$t)
      ind_0<-which(grid$s>grid$t)
      z[ind_0]<-0
      matrix(z,length(grid_s),length(grid_t))

    }
  }
  if(case=="Concurrent"){
    cat("Concurrent")
    beta<-function(grid_s,grid_t){

      f_1<-function(s,t){0.5*sin(30*s)+0.5*cos(30*t)}
      grid<-expand.grid(s=grid_s,t=grid_t)
      z<-f_1(grid$s,grid$t)
      ind_0<-which(grid$s!=grid$t)
      z[ind_0]<-0
      matrix(z,length(grid_s),length(grid_t))

    }
  }
  if(case=="Scenario HAT"){
    cat("Scenario HAT")
    beta<-function(grid_s,grid_t){
      eval_point<-expand.grid(grid_s,grid_t)
      eval_pdf<-dmvnorm(eval_point,mean = c(0.6,0.6), sigma = diag(c(0.001,0.001)))
      eval<--1+as.matrix(eval_point)%*%c(1.5,1.5)+0.05*eval_pdf
      Z<-matrix(eval,length(grid_s),length(grid_t))
    }
  }


  if(case=="Scenario DAMP"){
    cat("Scenario DAMP")
    beta<-function(grid_s,grid_t){
      eval_point<-as.matrix(expand.grid(grid_s,grid_t))
      eval<-5*exp(-5*eval_point%*%c(1,1))*cos(10*pi*eval_point%*%c(1,0))+5*exp(-5*eval_point%*%c(1,1))*cos(10*pi*eval_point%*%c(0,1))+1

      Z<-matrix(eval,length(grid_s),length(grid_t))
    }
  }
  if(case=="Scenario RCHANGE"){
    cat("Scenario RCHANGE")
    beta<-function(grid_s,grid_t){
      eval_point<-as.matrix(expand.grid(grid_s,grid_t))
      eval<-1-5/(1+exp(10*(eval_point%*%c(1,1)-0.2)))+5/(1+exp(75*(eval_point%*%c(1,1)-0.8)))

      Z<-matrix(eval,length(grid_s),length(grid_t))
    }
  }
  if(case=="Concurrent"){G<-t(eval.basis(grid_s,X_basis))%*%beta(grid_s,grid_t)}
  else{
    G<-(1/length(grid_s))*t(eval.basis(grid_s,X_basis))%*%beta(grid_s,grid_t)
  }
  Y_parz<-t(X_coef)%*%G

  signal_to_noise_ratio<-4
  ##for each observation SN=(sigma^2_signal/sigma^2_error)
  if(case=="Scenario I"){Y = Y_parz + Eps%*%diag(colVars(Eps)^(-1/2))}
  else{
    k <- sum(sqrt(colVars(Y_parz)))/sum((signal_to_noise_ratio*colVars(Eps)))

    Y = Y_parz + Eps%*%diag(rep(k,length(Y_parz[1,])))
  }

  out<-list(X=X,
            Y=t(Y),
            X_fd=X_fd,
            Eps=Eps,
            beta_matrix=beta(grid_s,grid_t))

  return(out)
}
## B splines penalized regression Ramsey 2005 -------------------------------


fr.usc.cv<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=10^seq(5,15,by=1),lambdas_t=10^seq(5,15,by=1)){

  n_obs<-length(X_fd$coefs[1,])


  inpr_vec<-list()
  comb_list<-expand.grid(lambdas_s,lambdas_t)



  parr_smooth<-function(ii){

    lambdas<-as.numeric(comb_list[ii,])
    lambda_s<-lambdas[1]
    lambda_t<-lambdas[2]
    ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
    split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))

    for(ll in 1:K){

      Y_i<-Y_fd[split_vec[[ll]]]
      X_i<-X_fd[split_vec[[ll]],]
      X_minus<-X_fd[-split_vec[[ll]],]
      Y_minus<-Y_fd[-split_vec[[ll]],]
      mod<-fregre.basis.fr(X_minus,Y_minus,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s,lambda.t = lambda_t)
      Y_hat<-predict(mod,X_i)
      inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
    }
    mean<-mean(unlist(inpr_vec))
    sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
    out<- list(mena=mean,
               sd=sd)

    return(out)

  }
  cores <- detectCores()
  vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_smooth,mc.cores = cores)
  par<-sapply(vec_par,"[[",1)
  sd<-sapply(vec_par,"[[",2)

  l_opt<-as.numeric(comb_list[max(which(par==min(par))),])


  lambda_s_opt<-l_opt[1]
  lambda_t_opt<-l_opt[2]
  mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s_opt,lambda.t = lambda_t_opt)

  B<-mod$beta.estbifd$coefs
  Beta_hat_fd<-mod$beta.estbifd

  out<-list(B=B,
            Beta_hat_fd=Beta_hat_fd,
            lambda_s_opt=lambda_s_opt,
            lambda_t_opt=lambda_t_opt,
            CV=par,
            CV_sd=sd,
            comb_list=comb_list,
            type="SMOOTH")
  return(out)
}
fr.usc<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=0,lambdas_t=0){

  n_obs<-length(X_fd$coefs[1,])
  mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambdas_s,lambda.t = lambdas_t)

  B<-mod$beta.estbifd$coefs
  Beta_hat_fd<-mod$beta.estbifd
  l_opt<-c(lambdas_s,lambdas_t)
  cat("SMOOTH:",l_opt,"     ")
  out<-list(B=B,
            Beta_hat_fd=Beta_hat_fd,
            lambda_s=lambdas_s,
            lambda_t=lambdas_t,
            mod=mod,
            Y=Y_fd,
            X=X_fd,
            type="SMOOTH")
  return(out)
}



## Smoothing adaptive spline -------------------------------




adass.fr_eaass <-
  function(Y_fd,
           X_fd,
           basis_s,
           basis_t,
           weight_fun_s = NULL,
           weight_fun_t = NULL,
           X_fd_test = NULL,
           Y_fd_test = NULL,
           rand_search_par = list(c(-4, 4), c(-4, 4), c(0,1,5,10, 15), c(0,1,2,3,4), c(0,1,5,10, 15), c(0,1,2,3,4)),
           popul_size = 12,
           iter_num = 10,
           grid_evals_der=NULL,
           grid_evalt_der=NULL,
           progress=TRUE,
           r=0.2,
           pert_vec=c(0.8, 1.2),
           ncores=1,
           ...) {

    inf_lam_s<-rand_search_par[[1]][1]
    sup_lam_s<-rand_search_par[[1]][2]
    inf_lam_t<-rand_search_par[[2]][1]
    sup_lam_t<-rand_search_par[[2]][2]
    inf_sum_s<-rand_search_par[[3]][1]
    sup_sum_s<-rand_search_par[[3]][2]
    inf_ind_s<-rand_search_par[[4]][1]
    sup_ind_s<-rand_search_par[[4]][2]
    inf_sum_t<-rand_search_par[[5]][1]
    sup_sum_t<-rand_search_par[[5]][2]
    inf_ind_t<-rand_search_par[[6]][1]
    sup_ind_t<-rand_search_par[[6]][2]


    list_l<-lapply(1:popul_size,function(ii){c(10^runif(1,inf_lam_s,sup_lam_s),10^runif(1,inf_lam_t,sup_lam_t),
                                          runif(1,inf_sum_s,sup_sum_s),runif(1,inf_ind_s,sup_ind_s),
                                          runif(1,inf_sum_t,sup_sum_t),runif(1,inf_ind_t,sup_ind_t))})

    ## Evaluation ----
    comb_list_eval<-do.call("rbind",list_l)
    mean_par<-NULL
    pb = txtProgressBar(min = 0, max = iter_num, initial = 0,style = 3)
    if(progress)
      setTxtProgressBar(pb,0)
    for (iter in 1:iter_num) {
      parr_eval <- function(ii) {
        lambda_s_i <- comb_list_eval[ii, 1]
        lambda_t_i <- comb_list_eval[ii, 2]
        sum_s_i <- comb_list_eval[ii, 3]
        ind_s_i <- comb_list_eval[ii, 4]
        sum_t_i <- comb_list_eval[ii, 5]
        ind_t_i <- comb_list_eval[ii, 6]

        mod_i <-
          adass.fr(
            Y_fd,
            X_fd,
            basis_s = basis_s,
            basis_t = basis_t,
            tun_par = c(lambda_s_i,lambda_t_i,sum_s_i,ind_s_i,sum_t_i,ind_t_i),
            weight_fun_s = weight_fun_s,
            weight_fun_t = weight_fun_t,
            X_fd_test = X_fd_test,
            Y_fd_test = Y_fd_test,
            grid_evals_der=grid_evals_der,
            grid_evalt_der=grid_evalt_der,
            CV = TRUE,...
          )

        a<-mod_i$CV
        b<-mod_i$CV_sd
        out <- list(a, b)
        return(out)
      }
      vec_par <-
        mclapply(seq(1, length(comb_list_eval[, 1])), parr_eval, mc.cores = ncores)
      mean_eval <- sapply(vec_par, "[[", 1)
      sd_eval <- sapply(vec_par, "[[", 2)

      ##Exploit (Truncation method) ----
      if (is.null(mean_par)) {
        mean_tot <- mean_eval
        sd_tot <- sd_eval
        comb_list_tot<-comb_list_eval
      }
      else{
        mean_tot <- c(mean_par, mean_eval)
        sd_tot <- c(sd_par, sd_eval)
      }
      if(progress&iter==iter_num){
        setTxtProgressBar(pb,iter)
        break
      }
      rank_vec <- rank(mean_tot)
      treshold <- length(mean_tot) * (1-r)
      in_vec <- which(rank_vec <= treshold)
      out_vec <- which(rank_vec > treshold)
      in_comb <- comb_list_tot[in_vec, ]
      out_comb <- comb_list_tot[out_vec, ]
      ns_out <- length(out_vec)
      ind_new_comb <- sample(1:length(in_vec), size = ns_out, replace = TRUE)
      new_comb <- in_comb[ind_new_comb, ]
      mean_par <- mean_tot[in_vec]
      sd_par <- sd_tot[in_vec]

      mean_tot <- c(mean_par, mean_par[ind_new_comb])

      ## Exploration (perturbation) ----
      fac_vec <- pert_vec
      ind_fac <- sample(1:2, length(new_comb), replace = T)
      fac_comb <- sapply(ind_fac, function(ii)
        fac_vec[ii])

      fac_mat <- matrix(fac_comb, ns_out, 6)

      new_comb <- new_comb * fac_mat

      comb_list_eval <- new_comb
      comb_list_tot <- rbind(in_comb, new_comb)
      if(progress)
        setTxtProgressBar(pb,iter)
    }

    ind_opt=which(mean_tot==min(mean_tot))
    lambda_s_opt <- comb_list_tot[ind_opt, 1]
    lambda_t_opt <- comb_list_tot[ind_opt, 2]
    sum_s_opt <- comb_list_tot[ind_opt, 3]
    ind_s_opt <- comb_list_tot[ind_opt, 4]
    sum_t_opt <- comb_list_tot[ind_opt, 5]
    ind_t_opt <- comb_list_tot[ind_opt, 6]


    tun_par_opt<-c(lambda_s_opt,lambda_t_opt,sum_s_opt,ind_s_opt,sum_t_opt,ind_t_opt)
    out<-list(tun_par_opt=tun_par_opt,
              CV=mean_tot,
              CV_sd=sd_tot,
              comb_list_tot=comb_list_tot,
              X=X_fd,
              Y=Y_fd)
    class(out)<-"adass_eaass"
    return(out)

  }

adass.fr<-function(Y_fd,X_fd,basis_s,basis_t,
                  weight_fun_s=NULL,weight_fun_t=NULL,grid_evals_der=NULL,
                  grid_evalt_der=NULL,tun_par=c(lambda_s=10^4,lambda_t=10^4,delta_s=0,gamma_s=1,delta_t=0,delta_t=1),
                  CV=FALSE,K=10,X_fd_test=NULL,Y_fd_test=NULL){


  n_obs<-dim(X_fd$coefs)[2]
  domain_s<-basis_s$rangeval
  domain_t<-basis_t$rangeval
  n_basis_s<-basis_s$nbasis
  n_basis_t<-basis_t$nbasis
  orders<-n_basis_s-length(basis_s$params)
  ordert<-n_basis_t-length(basis_t$params)
  breaks_s<-basis_s$params
  breaks_t<-basis_t$params
  ext_break_s<-c(rep(domain_s[1],orders),breaks_s,rep(domain_s[2],orders))
  ext_break_t<-c(rep(domain_s[1],ordert),breaks_t,rep(domain_s[2],ordert))

  ##Centering data
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  X_fd_cen<-center.fd(X_fd)
  Y_fd_cen<-center.fd(Y_fd)

  X_new<-inprod(X_fd_cen,basis_s)
  Y_new<-inprod(Y_fd_cen,basis_t)

  ##Eval weight function
  if(is.null(weight_fun_s)){
    beta_der_eval_s<-matrix(1,1,1)
    tun_par[3]=0
    tun_par[4]=0
    grid_evals_der=domain_s
  }
  else if(class(weight_fun_s)[1]=="bifd" ){
    if(is.null(grid_evals_der)|is.null(grid_evalt_der))
      stop("Please set the grid to evaluate the partial derivatives!")
    beta_der_eval_s<-eval.bifd(grid_evals_der,grid_evalt_der,weight_fun_s)
  }
  else{
    beta_der_eval_s<-weight_fun_s
    if(is.null(grid_evals_der)|is.null(grid_evalt_der)) stop("Please provide the grids to evaluate the partial derivatives!")
    if(sum(range(grid_evals_der)!=domain_s)|sum(range(grid_evalt_der)!=domain_t))
      stop("Grid domains to evaluate the partial derivatives different from the functional data domains!")
    if((length(grid_evals_der)!=dim(beta_der_eval_s)[1])|(length(grid_evalt_der)!=dim(beta_der_eval_s)[2]))
      stop("Wrong dimensions between grids and partial derivatives evaluation!")
  }
  if(is.null(weight_fun_t)){
    beta_der_eval_t<-matrix(1,1,1)
    tun_par[5]=0
    tun_par[6]=0
    grid_evalt_der=domain_t
  }
  else if(class(weight_fun_t)[1]=="bifd" ){
    if(is.null(grid_evals_der)|is.null(grid_evalt_der))
      stop("Please provide the grid to evaluate the partial derivatives!")
    beta_der_eval_t<-eval.bifd(grid_evals_der,grid_evalt_der,weight_fun_t)
  }
  else{
    beta_der_eval_t<-weight_fun_t
    if(is.null(grid_evals_der)|is.null(grid_evalt_der)) stop("Please provide the grids to evaluate the partial derivatives!")
    if(sum(range(grid_evals_der)!=domain_s)|sum(range(grid_evalt_der)!=domain_t))
      stop("Grid domains to evaluate the partial derivatives different from the functional data domains!")
    if((length(grid_evals_der)!=dim(beta_der_eval_t)[1])|(length(grid_evalt_der)!=dim(beta_der_eval_t)[2]))
      stop("Wrong dimensions between grids and partial derivatives evaluation!")
  }

  W_si<-R_si<-W_ti<-R_ti<-list()
  for (ii in 1:(length(grid_evals_der)-1)) {
    W_si[[ii]]<-eval.penalty(basis_s,rng =c(grid_evals_der[ii],grid_evals_der[ii+1]) )
    R_si[[ii]]<-eval.penalty(basis_s,2,rng =c(grid_evals_der[ii],grid_evals_der[ii+1]))
  }

  for (ii in 1:(length(grid_evalt_der)-1)) {
    W_ti[[ii]]<-eval.penalty(basis_t,rng =c(grid_evalt_der[ii],grid_evalt_der[ii+1]) )
    R_ti[[ii]]<-eval.penalty(basis_t,2,rng =c(grid_evalt_der[ii],grid_evalt_der[ii+1]) )
  }

  ## Calculation C1 and C2
  comb_list_mat<-expand.grid(1:(length(grid_evals_der)-1),1:(length(grid_evalt_der)-1))
  ##Cs ------
  par_Cs <- function(ll) {
    ii <- comb_list_mat[ll, 1]
    jj <- comb_list_mat[ll, 2]
    as.matrix.csr(W_ti[[jj]])%x% as.matrix.csr(R_si[[ii]])*(1 / (abs(beta_der_eval_s[ii, jj]) +tun_par[3]*max(abs(beta_der_eval_s))))^tun_par[4]
  }
  Cs_list<-lapply(seq(1,length(comb_list_mat[,1])),par_Cs)
  Cs<-Reduce("+",Cs_list)
  ##Ct ------
  par_Ct <- function(ll) {
    ii <- comb_list_mat[ll, 1]
    jj <- comb_list_mat[ll, 2]
    as.matrix.csr(R_ti[[jj]])%x% as.matrix.csr(W_si[[ii]])*(1/(abs(beta_der_eval_t[ii, jj]) +tun_par[5]*max(abs(beta_der_eval_t))))^tun_par[6]
  }
  Ct_list<-lapply(seq(1,length(comb_list_mat[,1])),par_Ct)
  Ct<-Reduce("+",Ct_list)


  if(CV){
    if(is.null(X_fd_test)){
      fun<-function(ii){
        ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
        split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
        mean_vec<-numeric()
        for(ll in 1:K){
          Y_i_fd<-Y_fd[split_vec[[ll]]]
          X_i<-X_new[split_vec[[ll]],]
          X_minus<-X_new[-split_vec[[ll]],]
          Y_minus<-Y_new[-split_vec[[ll]],]
          A<-W_t%x%(t(X_minus)%*%X_minus)+tun_par[1]*Cs+tun_par[2]*Ct
          B_s<-vec(t(X_minus)%*%Y_minus)
          b_vec<-SparseM::solve(A,B_s)
          B<-matrix(b_vec,basis_s$nbasis,basis_t$nbasis)
          Beta_hat_fd<-bifd(B,basis_s,basis_t)
          Y_hat<-fd(t(X_i%*%B),basis_t)
          mean_vec[ll]<-mean(diag(inprod(Y_i_fd-Y_hat,Y_i_fd-Y_hat)))
        }
        mean<-mean(mean_vec)
        sd<-sd(mean_vec)/sqrt(K)
        out<- list(mean=mean,
                   sd=sd)
        return(out)
      }
    }
    else{
      fun<-function(ii){
        A<-W_t%x%(t(X_new)%*%X_new)+tun_par[1]*Cs+tun_par[2]*Ct
        B_s<-vec(t(X_new)%*%Y_new)
        b_vec<-SparseM::solve(A,B_s)
        B<-matrix(b_vec,basis_s$nbasis,basis_t$nbasis)
        Beta_hat_fd<-bifd(B,basis_s,basis_t)
        mean<-get_PMSE(Y_fd_test,X_fd_test,Beta_hat_fd)
        out<- list(mean=mean,
                   sd=0)
        return(out)
      }
    }
    vec_par<-lapply(1,fun)
    CV<-sapply(vec_par,"[[",1)
    CV_sd<-sapply(vec_par,"[[",2)

  }
  A<-W_t%x%(t(X_new)%*%X_new)+tun_par[1]*Cs+tun_par[2]*Ct
  B_s<-vec(t(X_new)%*%Y_new)
  b_vec<-SparseM::solve(A,B_s)
  B<-matrix(b_vec,basis_s$nbasis,basis_t$nbasis)
  Beta_hat_fd<-bifd(B,basis_s,basis_t)

  if(CV==FALSE){
    CV=CV_sd=NA
  }
  out<-list(B=B,
            Beta_hat_fd=Beta_hat_fd,
            tun_par=tun_par,
            CV=CV,
            CV_sd=CV_sd,
            X=X_fd,
            Y=Y_fd)
  class(out)<-"adass"
  return(out)
}
get_PMSE<-function(Y_fd_test,X_fd_test,Beta){

  length_grid<-1000
  grid<-seq(0,1,length.out = length_grid)
  delta<-1/length_grid
  X_fd_eval<-t(eval.fd(grid,X_fd_test))
  Y_fd_eval<-t(eval.fd(grid,Y_fd_test))
  Beta_mat<-eval.bifd(grid,grid,Beta)

  Y_hat<-delta*X_fd_eval%*%Beta_mat

  PMSE<-mean(delta*rowSums((Y_fd_eval-Y_hat)^2))

  return(PMSE)
}
mean_rep<-function (fdobj,nobs)
{
  coef <- as.array(fdobj$coefs)
  coefd <- dim(coef)
  ndim <- length(coefd)
  basisobj <- fdobj$basis
  nbasis <- basisobj$nbasis
  coefmean <- apply(coef, 1, mean)
  mean_rep <- fd(matrix(rep(coefmean,nobs),coefd[1],nobs), basisobj)
  return(mean_rep)
}

