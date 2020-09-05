source("C:/R/mcmc.R")


# Define the distribution functions first: 

# probability density (can be calculated for an array of values):
fgev=function(x,theta)
{
    mu=theta[1]
    sigma=theta[2]
    chi=theta[3]

    if(sigma<0)      
      return(0)
    else
     {
      tx=1+(x-mu)/sigma*chi
      ind=1*(tx>0)
      tx=ind*tx+1*(1-ind)
      return(ind/sigma*tx^(-1/chi-1)*exp(-tx^(-1/chi)))
     }
}	

# Cumulative distribution:
Fgev=function(x,theta)
{
    mu=theta[1]
    sigma=theta[2]
    chi=theta[3]

    if(sigma<0)
       return(0)
    else
     {
       tx=1+(x-mu)/sigma*chi
       ind=1*(tx>0)
       tx=ind*tx+1*(1-ind)
       return(ind*exp(-tx^(-1/chi))+(chi<0)*(1-ind))
     }
}

# Quantile function (inverse cumulative function):
Fgev.inv=function(x,theta)
{
    mu=theta[1]
    sigma=theta[2]
    chi=theta[3]

    return(mu+sigma/chi*((-log(x))^(-chi)-1))
}

gev.ksidist.result=function(data.yearmax, # array of yearly maximums 
        T10.95.lower=3, T10.95.upper=600, # 95% cred. band for T10
	T100_T10.95.lower=3, T100_T10.95.upper=600, # ditto for T100-T10
	ksi.lower=0.106, ksi.upper=0.127, # 95% cred. band for ksi
        mcmc.numsamples=1000, # MCMC specification
        mcmc.burnin=10000, mcmc.spacing=10, mcmc.numtemp=1)
{
  y=data.yearmax
  
  # Normal distribution hyperparameters for ksi:
  ksi.mu=(ksi.lower+ksi.upper)/2
  ksi.sd=(ksi.upper-ksi.lower)/2/1.96 
      
  # A functiond fining how wrong the log-transformed hyperparameters 
  # for the T10 is:
  gamma1.error=function(hyper)
    {
	alpha=exp(hyper[1])
        beta=exp(hyper[2])
        ret=(qgamma(0.025,shape=alpha, scale=beta)-T10.95.lower)^2+
	    (qgamma(0.975,shape=alpha, scale=beta)-T10.95.upper)^2
        return(ret)
    }
  
  # Define the gamme distribution hyperparameters for T10:
  v1=optim(c(0,0),gamma1.error)
  alpha1=exp(v1$par[1])
  beta1=exp(v1$par[2])
  
  # Same for T100-T10:
  gamma2.error=function(hyper)
    {
	alpha=exp(hyper[1])
        beta=exp(hyper[2])
        ret=(qgamma(0.025,shape=alpha, scale=beta)-T100_T10.95.lower)^2+
	    (qgamma(0.975,shape=alpha, scale=beta)-T100_T10.95.upper)^2
        return(ret)
    }
  v2=optim(c(0,0),gamma2.error)
  alpha2=exp(v2$par[1])
  beta2=exp(v2$par[2])
  

  # Define the hyperparameters
  hyper.gev=list(alpha.T10=alpha1, beta.T10=beta1, 
	         alpha.T100_T10=alpha2, beta.T100_T10=beta2,
	         ksi.mu=ksi.mu, ksi.sd=ksi.sd)
	
  # Define the parameter structure. (Uses ls=log(sigma))
  par.gev=list(mu=0,ls=0,ksi=ksi.mu)
  	
  # Define transformation formulas from parameters to Q10(q1) and Q100-Q10(q2)
  q1=function(p1,mu,logsigma,ksi)
    mu+exp(logsigma)/ksi*((-log(1-p1))^(-ksi) -1)
  q2=function(p1,p2,mu,logsigma,ksi)
    exp(logsigma)/ksi*((-log(1-p2))^(-ksi) -1)-
    exp(logsigma)/ksi*((-log(1-p1))^(-ksi) -1)
 
  logprior.gev=function(par,hyp)
  {
    # Do a numerical derivation of the transformation formala from 
    # (mu,log(sigma)) to (T10, T100-T10)
    theta.der=array(0,c(2, 2))
    dx=0.0001
    p1=1/10
    p2=1/100
    q1.curr=q1(p1,par$mu,par$ls,par$ksi)
    q2.curr=q2(p1,p2,par$mu,par$ls,par$ksi)
    theta.der[1,1]=(q1.curr-q1(p1,par$mu-dx,par$ls,par$ksi))/dx
    theta.der[1,2]=(q1.curr-q1(p1,par$mu,par$ls-dx,par$ksi))/dx
    theta.der[2,1]=(q2.curr-q2(p1,p2,par$mu-dx,par$ls,par$ksi))/dx
    theta.der[2,2]=(q2.curr-q2(p1,p2,par$mu,par$ls-dx,par$ksi))/dx
    logdet=log(det(theta.der))
    
    if(is.na(logdet))
      return(-1e+200)
    if(!(logdet> -1e+200 & logdet<1e+200))
      return(-1e+200)
    
    if(is.na( dgamma(q1.curr,shape=hyp$alpha.T10,
	       scale=hyp$beta.T10,log=T)))
     {
       show(c(q1.curr, hyp$alpha.T10, hyp$beta.T10))
       return(-1e+200)
     }

    if(is.na( dgamma(q2.curr,shape=hyp$alpha.T100_T10,
	             scale=hyp$beta.T100_T10,log=T)))
     {
       show(c(q2.curr, hyp$alpha.T100_T10, hyp$beta.T100_T10))
       return(-1e+200)
     }

    ret=logdet + 
	dgamma(q1.curr,shape=hyp$alpha.T10,
	       scale=hyp$beta.T10,log=T) +
	dgamma(q2.curr,shape=hyp$alpha.T100_T10,
	       scale=hyp$beta.T100_T10,log=T) +
        dnorm(par$ksi,hyp$ksi.mu, hyp$ksi.sd, log=T)

    if(is.na(ret))
      return(-1e+200)
    if(!(ret> -1e+200 & ret<1e+200))
      return(-1e+200)

    return(ret) 
  }

  loglik.gev=function(x,co,par)
  {
    if((par$ksi>0 & sum(x<par$mu-exp(par$ls)/par$ksi)>0) | (par$ksi<0 & sum(x>par$mu-exp(par$ls)/par$ksi)>0))
    return(-1e+200)
  
    ret=sum(log(fgev(x,c(par$mu,exp(par$ls),par$ksi))))

    if(is.na(ret))
      return(-1e+200)

    if(ret> -1e+200 & ret< 1e+200)
      return(ret)

    return(-1e+200)
  }
  
  res.gev=sample.mcmc(loglik.gev,logprior.gev,y,NA,hyper.gev, par.gev,
    mcmc.numsamples, # number of MCMC samples retrieved
    mcmc.burnin,     # burnin
    mcmc.spacing,    # spacing of samples retrieved
    mcmc.numtemp)    # number of tempering chains

  mcmc.mu=as.numeric(res.gev[1,])
  mcmc.sigma=exp(as.numeric(res.gev[2,]))
  mcmc.ksi=as.numeric(res.gev[3,])

  par.median=list(mu=median(mcmc.mu),sigma=median(mcmc.sigma),ksi=median(mcmc.ksi))
 
  par.mean=list(mu=mean(mcmc.mu), sigma=mean(mcmc.sigma),ksi=mean(mcmc.ksi))
 
  par.95cred=list(mu=c(quantile(mcmc.mu,0.025),quantile(mcmc.mu,0.975)),
                  sigma=c(quantile(mcmc.sigma,0.025),quantile(mcmc.sigma,0.975)),
                  ksi=c(quantile(mcmc.ksi,0.025),quantile(mcmc.ksi,0.975)))

  ret=list(mcmc.mu=mcmc.mu, mcmc.sigma=mcmc.sigma, mcmc.ksi=mcmc.ksi,
           par.median=par.median, par.mean=par.mean, par.95cred=par.95cred)
}

gev.returnvalue.medianpar=function(gev.ksidist.res, return.time)
  Fgev.inv(1-1/return.time , c(gev.ksidist.res$par.median$mu,
                             gev.ksidist.res$par.median$sigma,
                             gev.ksidist.res$par.median$ksi))


gev.returnvalue.meanpar=function(gev.ksidist.res, return.time)
  Fgev.inv(1-1/return.time , c(gev.ksidist.res$par.mean$mu,
                             gev.ksidist.res$par.mean$sigma,
                             gev.ksidist.res$par.mean$ksi))

gev.returnvalue.quantile=function(gev.ksidist.res, return.time,quantile.value)
  {
    quantiles=array(NA,c(length(gev.ksidist.res$mcmc.ksi),length(return.time)))
    for(i in 1:length(gev.ksidist.res$mcmc.ksi))
     quantiles[i,]=Fgev.inv(1-1/return.time , 
                    c(gev.ksidist.res$mcmc.mu[i],gev.ksidist.res$mcmc.sigma[i],gev.ksidist.res$mcmc.ksi[i]))
    qq=function(x) quantile(x,quantile.value)
        
    apply(quantiles,2,qq)
  }

gev.returnvalue.mean=function(gev.ksidist.res, return.time)
  {
    quantiles=array(NA,c(length(gev.ksidist.res$mcmc.ksi),length(return.time)))
    for(i in 1:length(gev.ksidist.res$mcmc.ksi))
     quantiles[i,]=Fgev.inv(1-1/return.time , 
                    c(gev.ksidist.res$mcmc.mu[i],gev.ksidist.res$mcmc.sigma[i],gev.ksidist.res$mcmc.ksi[i]))
        
    apply(quantiles,2,mean)
  }