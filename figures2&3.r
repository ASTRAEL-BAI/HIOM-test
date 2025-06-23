#set.seed(1)
stoch_cusp=function(x,b,a,sdW,dt)
{
#  dt=.1
  x=x-dt*(x^3-b*x-a)+rnorm(1,0,sdW)
  x[is.nan(x)]=0
  x
}

stoch_cusp_run=function(x0,b,a,s,n,dt)
{
o=rep(0,n)
o[1]=x0
for(i in 2:n)
  o[i]=stoch_cusp(o[i-1],b,a,s,dt)
o
}

#pdf("figures/figure2.pdf",h=6,w=12)
dpi=300
tiff("figures/figure2.tiff",h=4,w=10,res=dpi,units='in')
layout(matrix(1:2,1,2))
set.seed(10)
dt=.15
s=.1
a=0;b=1
n=2000
x0=0
I_min=-.5
# illustion of stochastic cusp
plot(stoch_cusp_run(x0,b+I_min,a,s,n,dt)[-1:-50],type='l',bty='n',ylab='Opinion',xlab="time",ylim=c(-1.5,1.5))

n=400;dt=.2;bb=c(1.5,1,.5,0)
a=seq(-.5,.5,le=n)
a=c(a,rev(a))
bi=0
so=.01
for(b in bb)
{
  bi=bi+1
  o=rep(0,n)
  runs=1000
  for(i in 2:(2*n))
    o[i]=stoch_cusp_run(o[i-1],b+I_min,a[i],s=so,n=runs,dt)[runs] # o[i-1]-dt*(o[i-1]^3-b*o[i-1]-a[i])
  if(b==bb[1])plot(o[-1:-10]~a[-1:-10],xlim=c(-.65,.65),type='l',bty='n',ylab='Opinion',xlab="Information (I)",main="") else lines(o[-1:-10]~a[-1:-10],col=bi,lty=bi)
}
legend("bottomright",lty=1:4,col=1:4,legend=paste("A = ",bb),cex = 0.9,bty='n')

dev.off()



#pdf('figures/figure3.pdf',h=4,w=8)  # Dynamics of involvement
dpi=300
tiff("figures/figure3.tiff",h=4,w=8,res=dpi,units='in')

set.seed(6)
layout(matrix(1:2,1,2))
n=200

x=rep(0,n)
x[1]=.5
p=.05
d=.2
x_star=1

u=sample(0:1,n,T,prob=c(1-p,p))
u=rep(0,n)
u[c(10,50,90,110,115,130,155)]=1
for(i in 2:n)
  x[i]=x[i-1]+d*u[i]*(2*x_star-x[i-1])-2*d*x[i-1]*p
plot(x,type='l',bty='n',main="",xlab='time',ylab="Attention",ylim=c(0,1))
lines(u/15,col='red',lty=3)
legend(120,.3,lty=c(1,3),col=c(1,2),legend=c('Attention','Interaction'),bty='n',cex=.8)

rmin=.1
p=2
f=function(x,rmin,p) rmin + (1-rmin)/(1+exp(-p*x))
curve(f(x,rmin=0,p=1),-3,3,bty='n',xlab=expression('A'['i']-'A'['j']),ylab='r',ylim=c(0,1),lty=1)
curve(f(x,rmin=0,p=2),-3,3,bty='n',xlab=expression('A'['i']-'A'['j']),ylab='r',add=T,col=2,lty=2)
curve(f(x,rmin=.4,p=1),-3,3,bty='n',xlab=expression('A'['i']-'A'['j']),ylab='r',add=T,col=4,lty=3)
legend('bottomright',lty=c(1,2,3),col=c(1,2,4),legend=c('              ','             ','              '),bty='n')
text(1,.24,bquote("p = 1,"~r[min]~"= 0"),pos=4,cex=.8)
text(1,.15,bquote("p = 2,"~r[min]~"= 0"),pos=4,cex=.8)
text(1,0.05,bquote("p = 1,"~r[min]~"= .4"),pos=4,cex=.8)
dev.off()


