# Ising attitude

ham=function(x,n,t,w) -sum(t*x)-sum(w*x%*%t(x)/2) # hamiltonian 

learnw=function(x,w,d,e) # hebb rule for w
{ 
  for(i in 1:n)
    for(j in 1:n)
      w[i,j]=(1-e)*w[i,j]+d*x[i]*x[j]
    w
}

glauber=function(n,beta,t,w,x,iterations,plot=F,learning=F)
{
if(plot) xx=ww=hh=rep(0,iterations)  
for (j in 1:iterations)
  {
    i = sample(1:n,size=1) # consider neibor case
    x2=x;x2[i]=x2[i]*-1
    p=1/(1+exp(beta*(ham(x2,n,t,w)-ham(x,n,t,w))))  # metropolis
    if(runif(1)<p) x=x2# update state
    if(learning) w=learnw(x,w,.005,.005)
    if(plot) {xx[j]=sum(x);ww[j]=sum(w);hh[j]=ham(x,n,t,w)}
    
}
 if(plot) {layout(1:3);plot(xx);plot(ww);plot(hh)}
  list(x,w)
}

#pdf('figures/figureA1.pdf')
#dpi=300
tiff("figures/figureA1.tiff",height=7,width=5,res=dpi,units='in')


# Simulation

layout(matrix(c(1:18),6,3,byrow=T),heights = c(1,1,1,1,1,1.4))
n=40 # 40 nodes
iterations=500 # 500 # interations
m=200 #200
############## Pitchfork

for(pl in  1:6)
{
  if(pl==1)
  {p_missing_link=0
  mw=.1
  mt=0
  sdw=0
  sdt=0
  sym=T}
  
  if(pl==2)
  {p_missing_link=0
  mw=.1
  mt=0
  sdw=0.1
  sdt=0
  sym=T}
  
  
  if(pl==3)
  {p_missing_link=0.5
  mw=.1
  mt=0
  sdw=0.1
  sdt=0
  sym=T}
  
  if(pl==4)
  {p_missing_link=0.5
  mw=.1
  mt=0
  sdw=0.1
  sdt=0
  sym=F}
   
  if(pl==5)
  {p_missing_link=0
  mw=0
  mt=0
  sdw=0.1
  sdt=0
  sym=F}
  
  
  if(pl==6)
  {p_missing_link=0.95
  mw=.1
  mt=0
  sdw=0.0
  sdt=0.0
  sym=T
  
  }

w=matrix(rnorm(n*n,mw,sdw),n,n) # biased weights
w=(1-diag(n))*w
if(sym) w[lower.tri(w)] = t(w)[lower.tri(w)]
wm=matrix(sample(0:1,n*n,T,prob=c(p_missing_link,1-p_missing_link)),n,n)
w=wm*w


t=matrix(rnorm(n*n,mt,sdt),n,n) #tresholds

ww=w
ww[1-diag(n)==0]=NA
minw=min(ww-.1,-.5,na.rm=T)
mint=min(t-.1,-.5,na.rm=T)
maxw=max(ww+.1,.5,na.rm=T)
maxt=max(t+.1,.5,na.rm=T)

if(pl<6) par(mar=c(1.5,2.5,1.5,2)) else par(mar=c(5,2.5,1.5,2))


if(pl==6) xlab=expression(omega) else xlab=''
hist(ww,breaks=seq(minw,maxw,.05),xlim=c(minw,maxw),col='grey',main='',xlab=xlab,ylab='',freq=F,axes=F,cex.lab=1.2)
axis(1)
# if(pl==6) xlab=expression(tau) else xlab=''
# hist(t,breaks=seq(mint,maxt,.05),col='grey',xlim=c(mint,maxt),main='',ylab='',xlab=xlab,freq=F,axes=F)
# axis(1)

data=rep(0,m)
betas=seq(0,2,length=m)
i=0
for(beta in betas)
{
i=i+1
x=sample(c(-1,1),n,T)
data[i]=sum(glauber(n,beta,t,w,x,iterations)[[1]])
}
if(pl==6) xlab='A' else xlab=''
plot(data~betas,type='p',xlab=xlab,ylab=expression(sum(X)),bty='n',pch=20,cex.lab=1.2,cex.axis=.9)

data=rep(0,m)
i=0
beta=4
tau=seq(-.2,.2,length=.5*m)
tau=c(tau,rev(tau))

for(t_common in tau)
{
  i=i+1
  if(i%%2==0) ph=.2 else ph=.8
  x=sample(c(-1,1),prob=c(ph,1-ph),n,T)
  data[i]=sum(glauber(n,beta,t+t_common,w,x,iterations)[[1]])
}

if(pl==6) xlab="I" else xlab=''
plot(data~tau,type='p',xlab=xlab,ylab=expression(sum(X)),bty='n',pch=20,cex.lab=1.2)
}
dev.off()

