
##################################################
# Adaptive RW Metropolis with parallel tempering 
##################################################

# Logarithmic proability density=log(likelihood*prior):
logprob=function(loglik,logprior, response, covariate, params, hyperparam, T)
{
  p=logprior(params, hyperparam)
  if(p<=-1e+199 || is.na(p))
    {
      return(-1e+200)
    }

  l=loglik(response, covariate, params)
  if(l<=-1e+199 || is.na(l))
    {
      return(-1e+200)
    }

  ret=(l+p)/T

  return(ret)
}


# Single RW Metropolis sampling. Uses tempering algorithm also
# params should be a collection of parameters, one for eahc value
# in the T array. Returns an object containing new parameter
# values for each T, an array of the new logprobs and updated
# acceptance numbers (for adaptive usages)
# rw is an array of random walk standard deviations for each parameter
# Acceptances and number of tests (number of possible acceptances)
# for each parameter, are also only for for the first chain
# (the chain which is the object of interest)
# Input: likelihood
sample.once=function(loglik,logprior, response, covariate, params, rw, logprobs, acceptances, num.tests, hyperparam, swaps, T)
  {
    numtemp=length(T) # fetch the number of tempering chains
    #i630=FALSE

    
    if(numtemp>1 && runif(1)<0.5) # with probability 0.1, do parallell tempering
      {
        index=sample(1:(numtemp-1),1)
        new.logprob1=logprob(loglik,logprior, response, covariate, params[,index], hyperparam, T[index+1])
        new.logprob2=logprob(loglik,logprior, response, covariate, params[,index+1], hyperparam, T[index])

        # accept switching of tempering chains?
	if(is.finite(new.logprob1) & is.finite(new.logprob2))
         if(log(runif(1))<new.logprob1+new.logprob2-logprobs[index]-logprobs[index+1])
          {
            swaps[index]=swaps[index]+1
            
            param.buffer=params[,index]
            params[,index]=params[,index+1]
            params[,index+1]=param.buffer

            logprobs[index]=new.logprob2
            logprobs[index+1]=new.logprob1
          }
      }
    else
      {
        for(t in 1:numtemp) # traverse the parallell tempering chains
          {
            rw=as.numeric(rw)

            pars=unlist(params[,t])
            for(i in 1:length(pars)) # traverse the parameters
              {
                #if(i630 & t==1 && i==631)
		#{
		#  show("i=631")
		#  show(unlist(params[,1])[630])
	        #}
		
                # suggest a new parameter value:
                new.param=pars
                new.param[i]=as.numeric(new.param[i])+sqrt(T[t])*rw[i]*rnorm(1)
                new.param=relist(new.param,params[,t])

                # fetch the logprob:
                new.logprob=logprob(loglik,logprior, response, covariate,
                  new.param, hyperparam, T[t])

                # accept/reject:
                if(is.finite(new.logprob))
                 if(log(runif(1))<new.logprob-logprobs[t])
                  {
                    logprobs[t]=new.logprob
                    params[,t]=new.param
		    pars=unlist(params[,t])
                    if(t==1)
                      {
                        acceptances[i]=acceptances[i]+1
			#if(i==630)
			# {
 			#   show(c(unlist(new.param)[i],unlist(params[,t])[i]))
			#   i630=TRUE
 			# }
                      }
                  }
                
                # update the number of tests if this is the original chain:
                if(t==1)
                  {
                    num.tests[i]=num.tests[i]+1
                  }
              }
          }
      }

    #show(params[1:10,1])
    #show(loglik(response,covariate,params[,1]))
    
    # return structure:
    ret=list(logprobs=logprobs, params=params,
      acceptances=acceptances, num.tests=num.tests, swaps=swaps)

    #if(i630)
     # show(c(1,unlist(params[,1])[630], unlist(ret$params[,1])[630]))
    return(ret)
  }








# Repeated sampling
#
# INPUT:
# loglik - logarithmic likelihood, should take response,
#          covariates and parameters as input
# logprior - logarithmic prior, should take parameres and
#            hyperparameters as input
# response - response variable
# covariates - if applicable
# hyperparam - hyper parameters, used in the prior
# param.template - gives and indication of how the parameter
#                  vector should look like
# numsample - number of samples to be returned
# burnin - burnin period length
# indep - spacing between samples to be included in the output
# numtemp - number of "temperatures" (minimum 1)
#
# OUTPUT: 
# Returns a numsamples x numparams array containing parameter samples 
sample.mcmc=function(loglik, logprior, response, 
                     covariate, hyperparam, param.template,
  	                numsamples, burnin=1000, 
                     indep=10, numtemp=1, 
                     base.temp=1.25, talkative=F)
  {
    numparam=length(unlist(param.template))
    
    # initialization of tempering chains:
    T=base.temp^(0:(numtemp-1))
    params=cbind(param.template)
    if(numtemp>=2)
      {
        for(t in 2:numtemp)
          {
            params=cbind(params,param.template)
          }
      }
    
    # initialization of random walk step lengths:
    rw=unlist(param.template)
    rw[1:length(rw)]=rep(1,length(rw))
    
    # Initialization of number of acceptances and tests:
    acc=rep(0,numparam)
    num=rep(0,numparam)
    swaps=rep(0,numtemp)
    
    # Initialization of logarithmic probability densities:
    lp=rep(0,numtemp)
    for(t in 1:numtemp)
      {
        lp[t]=logprob(loglik,logprior, response,covariate,params[,t],hyperparam,T[t])
      }

    numburn=max(2,ceiling(burnin/5000))
    # burnin:
    for(i in 1:numburn)
      {
         if(talkative)
           show(c("Burnin - non-adaptive rw steps",i,"of",numburn))
        # non-adaptive part:
	if(ceiling(burnin/numburn/2)>=1)
         for(j in 1:ceiling(burnin/numburn/2))
          {
            res=sample.once(loglik,logprior,response, covariate,
	                    params, rw, lp, acc, num,
                            hyperparam, swaps, T)
            params=res$params
	    #if(numparam>630)
	    #  if(acc[630]!=res$acceptances[630])
	    #    show(unlist(params[,1])[630])
            acc=res$acceptances
            num=res$num.tests
            lp=res$logprobs
          }

        if(talkative)
          show(c("Burnin - adaptive rw steps",i,"of",numburn))
        # adaptive part:
	if(ceiling(burnin/numburn/200)>=1)
         for(k in 1:ceiling(burnin/numburn/200))
          {
            # Run the chain 100 times, recoding the number of acceptances
            # of each parameter:
            acc=rep(0,numparam)
            num=rep(0,numparam)
            for(j in 1:100)
              {
                res=sample.once(loglik,logprior,response, covariate, params,
                  rw, lp, acc, num, hyperparam, swaps, T)
                params=res$params
                acc=res$acceptances
                num=res$num.tests
                lp=res$logprobs
              }

            # Update the random walk step lengths:
            for(i in 1:length(rw))
              {
                rw[i]=as.numeric(rw[i])*exp((acc[i]/num[i]-0.33)*2)
              }

            # debug:
            if(talkative)
            {
              show("rw (debug):")
              show(as.numeric(rw))
              show("acc (debug):")
              show(as.numeric(acc))
	      show(params[,1])
            }
          }
      }

    if(talkative)
      show("start sampling...")
    
    # sampling:
    for(i in 1:numsamples)
      {
        # run the chain 'indep' number of times:
        for(j in 1:indep)
          {
            res=sample.once(loglik,logprior,response,
	      covariate, params, rw, lp, acc,
              num, hyperparam, swaps, T)
            params=res$params
            acc=res$acceptances
            num=res$num.tests
            lp=res$logprobs
            swaps=res$swaps
          }

        # add the current parameter to the set of samples:
        if(i==1)
          {
            ret=cbind(params[,1])
          }
        else
          {
            ret=cbind(ret,params[,1])
          }

        if(i%%10==0 & talkative)
          {
            show(i)
            show(params[,1])
          }
      }

    if(numtemp>1 & talkative)
      {
        show("Swaps:")
        show(swaps[1:(numtemp-1)])
      }
    
    return(ret)
  }




library(MASS)
# Next is to calculate model likelihoods:
model.loglik=function(loglik, logprior, response, covariate,
  hyperparam, param.template, samples, num.iterations)
{
  numparam=dim(samples)[1]
  numsamp=dim(samples)[2]

  samp=array(NA,c(numparam,numsamp))
  for(i in 1:numparam)
    samp[i,]=as.numeric(samples[i,])

  m=array(NA,numparam)
  for(i in 1:numparam)
    m[i]=mean(samp[i,])
  m=c(m)
  S=var(t(samp))
  Sinv=solve(S)
  logSdet=as.numeric(determinant(S,logarithm=TRUE)$modulus)
  
  param=param.template
  lx0=-0.5*numparam*log(2*pi)-0.5*logSdet
  param[1:length(param)]=m
  lp0=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)

  mm=max(10,ceiling(num.iterations/100))
  xp=mvrnorm(mm,m,S)

  for(i in 1:mm)
    {
      lx=-0.5*numparam*log(2*pi)-0.5*logSdet-0.5*(xp[i,]-m)%*%Sinv%*%(xp[i,]-m)
      param=param.template
      param[1:length(param)]=xp[i,]
      lp=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)
      if(!is.infinite(lx) & !is.infinite(lp) & lp+lx > lp0+lx0)
        {
          lp0=lp
          lx0=lx
        }
    }

  m.lik=0
  xp=mvrnorm(num.iterations,m,S)
  for(i in 1:num.iterations)
    {
      lx=-0.5*numparam*log(2*pi)-0.5*logSdet-0.5*(xp[i,]-m)%*%Sinv%*%(xp[i,]-m)
      param=param.template
      param[1:length(param)]=xp[i,]
      lp=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)
      w=(lp-lp0)-(lx-lx0)
      m.lik=m.lik+exp(w)
    }
  m.lik=m.lik/num.iterations

  log(m.lik)+lp0-lx0
}



# Next is to calculate model likelihoods:
model.loglik.indep=function(loglik, logprior, response, covariate,
  hyperparam, param.template, samples, num.iterations)
{
  numparam=dim(samples)[1]
  numsamp=dim(samples)[2]

  samp=array(NA,c(numparam,numsamp))
  for(i in 1:numparam)
    samp[i,]=as.numeric(samples[i,])

  m=array(NA,numparam)
  s=array(NA,numparam)
  for(i in 1:numparam)
  {
    m[i]=mean(samp[i,])
    s[i]=sd(samp[i,])
  }
  
  param=param.template
  lx0=sum(dnorm(m,m,s,log=T))
  param[1:length(param)]=m
  lp0=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)

  mm=max(1000,ceiling(num.iterations/5))

  
  m.use=rep(1,mm)%*%t(m)
  s.use=rep(1,mm)%*%t(s)
  xp=matrix(rnorm(mm*numparam),nrow=mm)
  xp=m.use+s.use*xp

  #xp=mvrnorm(mm,c(m),diag(s))

  for(i in 1:mm)
    {
      lx=sum(dnorm(xp[i,],m,s,log=T))
      param=param.template
      param[1:length(param)]=xp[i,]
      lp=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)
      if(!is.infinite(lx) & !is.infinite(lp) & lp-lx > lp0-lx0)
        {
          lp0=lp
          lx0=lx
	  #show(i)
        }
    }

  m.lik=0
  #xp=mvrnorm(num.iterations,c(m),diag(s))
  m.use=rep(1,num.iterations)%*%t(m)
  s.use=rep(1,num.iterations)%*%t(s)
  xp=matrix(rnorm(num.iterations*numparam),nrow=num.iterations)
  xp=m.use+s.use*xp

  for(i in 1:num.iterations)
    {
      lx=sum(dnorm(xp[i,],m,s,log=T))
      param=param.template
      param[1:length(param)]=xp[i,]
      lp=logprob(loglik,logprior, response, covariate, param, hyperparam, 1)
      w=(lp-lp0)-(lx-lx0)
      m.lik=m.lik+exp(w)
      #show(c(m.lik,w))
    }
  m.lik=m.lik/num.iterations

  log(m.lik)+lp0-lx0
}

