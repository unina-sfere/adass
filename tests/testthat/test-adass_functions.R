

case<-"Scenario HAT"
n_obs_tra<-10
n_obs_test<-4000
n_obs<-n_obs_tra+n_obs_test

data<-simulate_data(case,n_obs=n_obs)
Beta_vero_mat<-data$beta_matrix
X<-data$X
Y<-data$Y

domain<-c(0,1)
length_grid<-500
norder<-4
grid<-seq(0,1,length.out = length_grid)


# Smoothing
n_basis_x<-min(80,length_grid)
n_basis_y<-min(80,length_grid)
breaks_x<-seq(0,1,length.out = (n_basis_x-2))
breaks_y<-seq(0,1,length.out = (n_basis_y-2))
basis_x <- create.bspline.basis(domain,breaks=breaks_x)
basis_y <- create.bspline.basis(domain,breaks=breaks_y)
X_fd <- smooth.basis(grid,X,basis_x)$fd
Y_fd <- smooth.basis(grid,Y,basis_y)$fd
inter_basis<-create.bspline.basis(domain,nbasis = length(grid),norder = 3)


Beta_vero_fd<-bifd(Beta_vero_mat,inter_basis,inter_basis)

# Basis s and t
n_basis_s<-min(10,length_grid)
n_basis_t<-min(10,length_grid)
breaks_s<-seq(0,1,length.out = (n_basis_s-2))
breaks_t<-seq(0,1,length.out = (n_basis_t-2))
basis_s <- create.bspline.basis(domain,breaks=breaks_s,norder = 4)
basis_t <- create.bspline.basis(domain,breaks=breaks_t,norder = 4)

# Matrices
W_s<-eval.penalty(basis_s)
W_s_sp<-sparseMatrix(which(W_s>0,arr.ind=T)[,1],which(W_s>0,arr.ind=T)[,2],x=W_s[which(W_s>0,arr.ind=T)])
W_t<-eval.penalty(basis_t)
W_t_sp<-sparseMatrix(which(W_t>0,arr.ind=T)[,1],which(W_t>0,arr.ind=T)[,2],x=W_t[which(W_t>0,arr.ind=T)])
W_st<-inprod(basis_s,basis_t)
W_st_sp<-sparseMatrix(which(W_st>0,arr.ind=T)[,1],which(W_st>0,arr.ind=T)[,2],x=W_st[which(W_st>0,arr.ind=T)])
R_s<-eval.penalty(basis_x,2)
R_s_sp<-sparseMatrix(which(R_s!=0,arr.ind=T)[,1],which(R_s!=0,arr.ind=T)[,2],x=R_s[which(R_s!=0,arr.ind=T)])
R_t<-eval.penalty(basis_t,2)
R_t_sp<-sparseMatrix(which(R_t!=0,arr.ind=T)[,1],which(R_t!=0,arr.ind=T)[,2],x=R_t[which(R_t!=0,arr.ind=T)])

# Training and Test set
X_fd_tra_nc<-X_fd[1:n_obs_tra]
Y_fd_tra_nc<-Y_fd[1:n_obs_tra]
X_fd_tra<-center.fd(X_fd_tra_nc)
Y_fd_tra<-center.fd(Y_fd_tra_nc)
X_fd_test<-X_fd[(n_obs_tra+1):(n_obs)]-mean_rep(X_fd_tra_nc,n_obs_test)
Y_fd_test<-Y_fd[(n_obs_tra+1):(n_obs)]-mean_rep(Y_fd_tra_nc,n_obs_test)
X_coef_tra<-t(X_fd_tra$coefs)
Y_coef_tra<-t(Y_fd_tra$coefs)
X_coef_test<-t(X_fd_test$coefs)
Y_coef_test<-t(Y_fd_test$coefs)


# Smooth regression (fda) --------------------------------------------------------
mod_smooth_cv<-fr.usc.cv(Y_fd_tra,X_fd_tra,basis_s,basis_t,K=10,lambdas_s = 10^seq(-8,-2),lambdas_t = 10^seq(-8,-2))
mod_smooth<-fr.usc(Y_fd_tra,X_fd_tra,basis_s,basis_t,lambdas_s=10^-2,lambdas_t =10^-2)
ISE_smooth<-get_ISE(mod_smooth$Beta_hat_fd,Beta_vero_fd,case)
PMSE_smooth<-get_PMSE(Y_fd_test,X_fd_test,mod_smooth$Beta_hat_fd)


# Smooth adaptive ---------------------------------------------------------

grid_s<-seq(0,1,length.out = 10)
grid_t<-seq(0,1,length.out = 10)
inter_basis2<-create.bspline.basis(domain,nbasis = length(grid_s),norder = 3)
beta_der_vero_s<-eval.bifd(grid_s,grid_t,Beta_vero_fd,sLfdobj = 2)
beta_der_vero_t<-eval.bifd(grid_s,grid_t,Beta_vero_fd,tLfdobj = 2)
beta_der_vero_s_fd<-bifd(beta_der_vero_s,inter_basis2,inter_basis2)
beta_der_vero_t_fd<-bifd(beta_der_vero_t,inter_basis2,inter_basis2)
beta_der_eval_s<-eval.bifd(grid_s,grid_t,mod_smooth$Beta_hat_fd,sLfdobj = 2)
beta_der_eval_t<-eval.bifd(grid_s,grid_t,mod_smooth$Beta_hat_fd,tLfdobj = 2)


ISE_der_s<-get_ISE(bifd(beta_der_eval_s,inter_basis2,inter_basis2),beta_der_vero_s_fd,case)
ISE_der_t<-get_ISE(bifd(beta_der_eval_t,inter_basis2,inter_basis2),beta_der_vero_t_fd,case)
mod_adsm<-adass.fr_eaass(Y_fd_tra,X_fd_tra,basis_s,basis_t,
                      weight_fun_s=beta_der_eval_s, weight_fun_t=beta_der_eval_t,X_fd_test=X_fd_tra[1:10],Y_fd_test=X_fd_tra[1:10],
                      rand_search_par=list(c(-8,4),c(-8,4),c(0,0.1),c(0,4),c(0,0.1),c(0,4)), grid_evals_der=grid_s,
                      grid_evalt_der=grid_t,
                      popul_size = 2,ncores=1,iter_num=1)


mod_opt <-
  adass.fr(
    Y_fd_tra,
    X_fd_tra,
    basis_s = basis_s,
    basis_t = basis_t,
    tun_par=mod_adsm$tun_par_opt,
    weight_fun_s = beta_der_eval_s,
    weight_fun_t = beta_der_eval_t,
    grid_evals_der=grid_s,
    grid_evalt_der=grid_t,
    CV = FALSE
  )
