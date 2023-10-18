#' Title Simulation data generation with mixture_noise.
#' @param lambda The probability of the sample() function process.
#' @param spots The number of spots.
#' @param gns the number of SVGs，
#' @param type A character.'ZINB' or 'ZIP','Pois', or "NB".
#' @param mix the mixture proportion
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @return a data frame
#' @export
mix_zi=function(lambda = lambda, mix=0.1,spots, gns, type, ptn,png_dir) {
  n=spots;lambda=lambda
  ## type='logN','Pois','ZIP','ZINB','NB'
  if (type=='logN') {
    zi.p=0.03;FC=2
    mu1=1.5;cut=5;sd1=1;sd2=1
    df=HS.1004.newsimu::noise.norm(noise=0,cut=cut,spots = n,se=gns,ns=gns,se.mu=mu1*FC,ns.mu = mu1,
                                   lambda=lambda,se.p = zi.p/3,ns.p=zi.p,se.sd=sd1,ns.sd = sd2,
                                   ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='Pois'){
    zi.p=0.03
    high=1.5;low=0.5
    df=noise.zi(noise=0,spots = n,se=gns,ns=gns,type='ZIP',lambda=lambda,
                se.mu=high,ns.mu =low,#se.mu=mu1*FC,ns.mu = mu1,
                se.p = zi.p/3,ns.p=zi.p, ptn=ptn,png_dir =  ptn_dir)
  } else if(type=='ZIP'){
    zi.p=0.6; high=6;low=2
    df=noise.zi(noise=0,spots = n,se=gns,ns=gns,type='ZIP',lambda=lambda,
                se.mu=high,ns.mu =low,#se.mu=mu1*FC,ns.mu = mu1,
                se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='ZINB'){
    zi.p=0.8;size1=0.5;mu1=0.5;FC=2
    df=noise.zi(noise=0,spots = n,se=gns,ns=gns,type='ZINB',lambda=lambda,
                se.mu=mu1*FC,ns.mu = mu1, se.size = size1*FC,ns.size = size1,
                se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='NB'){
    zi.p=0.03;size1=1.5;mu1=0.5;FC=3
    df=noise.zi(noise=0,spots = n,se=gns,ns=gns,type='ZINB',lambda=lambda,
                se.mu=mu1*FC,ns.mu = mu1, se.size = size1*FC,ns.size = size1,
                se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)

  }

  # df=new.zi(lambda = lambda, spots=spots, se=gns,ns=gns, type=type,se.p=se.p,se.size=se.size,
  #             se.mu=se.mu, ns.p=ns.p, ns.size = ns.size, ns.mu=ns.mu, ptn=ptn,png_dir=png_dir)
  d_1=df[,3:(gns+2)]
  d_2=df[,(gns+3):(2*gns+2)]
  # new_d=data.frame()
  mix_s=seq(0.1,mix,by=0.1)
  mix_s=c(mix_s,1-mix_s)
  for (i in mix_s) {
    dd=d_1*(1-i)+d_2*i
    a=ifelse(i<=mix,'tt','ff')
    cat('i======',i,'\n')
    cat('a======',a,'\n')
    colnames(dd)<-paste(a,i*10,1:ncol(dd),sep='_')
    # colnames(dd)<-paste(a,1:ncol(dd),sep='_')

    df=cbind(df,dd)
  }
  end_d=df
}


#' Title For irregular pattern, Simulation data generation with mixture_noise.
#' Title Simulation data generation with mixture_noise.
#' @param lambda The probability of the sample() function process.
#' @param spots The number of spots.
#' @param gns the number of SVGs，
#' @param type A character.'ZINB' or 'ZIP','Pois', or "NB".
#' @param mix the mixture proportion
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @return a data frame
#' @export

mix_ire=function(lambda = lambda, mix=0.1,spots, gns, type, ptn,png_dir) {
  n=spots;lambda=lambda
  ## type='logN','Pois','ZIP','ZINB','NB'
  if (type=='logN') {
    zi.p=0.03;FC=2
    mu1=1.5;cut=5;sd1=1;sd2=1
    df=HS.1004.newsimu::ire.norm(noise=0,cut=cut,spots = n,se=gns,ns=gns,se.mu=mu1*FC,ns.mu = mu1,
                                 lambda=lambda,se.p = zi.p/3,ns.p=zi.p,se.sd=sd1,ns.sd = sd2,
                                 ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='Pois'){
    zi.p=0.03
    high=1.5;low=0.5
    df=ire.noise(noise=0,spots = n,se=gns,ns=gns,type='ZIP',lambda=lambda,
                 se.mu=high,ns.mu =low,#se.mu=mu1*FC,ns.mu = mu1,
                 se.p = zi.p/3,ns.p=zi.p, ptn=ptn,png_dir =  ptn_dir)
  } else if(type=='ZIP'){
    zi.p=0.6; high=6;low=2
    df=ire.noise(noise=0,spots = n,se=gns,ns=gns,type='ZIP',lambda=lambda,
                 se.mu=high,ns.mu =low,#se.mu=mu1*FC,ns.mu = mu1,
                 se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='ZINB'){
    zi.p=0.8;size1=0.5;mu1=0.5;FC=2
    df=ire.noise(noise=0,spots = n,se=gns,ns=gns,type='ZINB',lambda=lambda,
                 se.mu=mu1*FC,ns.mu = mu1, se.size = size1*FC,ns.size = size1,
                 se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)
  } else if (type=='NB'){
    zi.p=0.03;size1=1.5;mu1=0.5;FC=3
    df=ire.noise(noise=0,spots = n,se=gns,ns=gns,type='ZINB',lambda=lambda,
                 se.mu=mu1*FC,ns.mu = mu1, se.size = size1*FC,ns.size = size1,
                 se.p = zi.p/3,ns.p=zi.p,ptn=ptn,png_dir =  ptn_dir)

  }

  # df=new.zi(lambda = lambda, spots=spots, se=gns,ns=gns, type=type,se.p=se.p,se.size=se.size,
  #             se.mu=se.mu, ns.p=ns.p, ns.size = ns.size, ns.mu=ns.mu, ptn=ptn,png_dir=png_dir)
  d_1=df[,3:(gns+2)]
  d_2=df[,(gns+3):(2*gns+2)]
  # new_d=data.frame()
  mix_s=seq(0.1,mix,by=0.1)
  mix_s=c(mix_s,1-mix_s)
  for (i in mix_s) {
    dd=d_1*(1-i)+d_2*i
    a=ifelse(i<=mix,'tt','ff')
    cat('i======',i,'\n')
    cat('a======',a,'\n')
    colnames(dd)<-paste(a,i*10,1:ncol(dd),sep='_')
    # colnames(dd)<-paste(a,1:ncol(dd),sep='_')

    df=cbind(df,dd)
  }
  end_d=df
}

