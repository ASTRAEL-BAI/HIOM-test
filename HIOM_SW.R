rm(list=ls())

library(igraph)
library(colorspace)
library(ggplot2)
library(diptest)
library(png)
library(qgraph)
library(tidyverse)

hartigan_d= function(x)
{
  dip=dip.test(x)
  dip.p=dip$p.value
  dip=dip$statistic
  dipstar=4-as.numeric(cut(dip.p,c(-Inf,.001,.01,.05,Inf)))
  if(dipstar!=0) dipstar=paste(rep('*',dipstar),sep='',collapse="") else dipstar=""
  return(list(dip,dipstar))
}

sample = function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}

information_update=function(i1,i2,a1,a2,persuasion,r_min,info_update=T)
{
  r=r_min+(1-r_min)/(1+exp(-1*persuasion*(a1-a2))) # resistance
  inf=r*i1+(1-r)*i2  # compute weighted man of information
  if(!info_update) inf=i1
  return(inf)
}

stoch_cusp=function(N,x,b,a,s_O,maxwell_convention=F) 
{
  dt=.01
  x=x-dt*(x^3-b*x-a)+rnorm(N,0,s_O) # updates opinion according to cusp
  if(maxwell_convention) x[x*a<0]=-x[x*a<0]
  x[is.nan(x)]=0
  x
}

make_network = function(model='SBM',N=400,clusters=10,p_within=.01,
                        p_between=.001,rewiring=.02)
{
  if(model=='WS')
  {
    g <- gg <- sample_smallworld(1, N,  1,rewiring); # The Watts-Strogatz small-world model
    l <-layout_nicely(g)
    adj <- get.adjacency(g)
  }
  
  if(model=='SBM')
  {
    pm=matrix(p_between,clusters,clusters);
    diag(pm)=p_within
    g <- gg <- sample_sbm(N, pref.matrix=pm, block.sizes=rep(N/clusters,clusters))
    l <-layout_nicely(g)
    adj <- get.adjacency(g)
  }
  
  if(model=='lattice')
  {
    g=gg=make_lattice(c(sqrt(N), sqrt(N), 1),nei=1);
    l = layout_on_grid(g);
    adj <- get.adjacency(g)
  }
  
  if(model=='SBM_4')
  {  
    pm=matrix(p_between,clusters,clusters);
    diag(pm)=p_within
    g <- gg <- sample_sbm(.5*N, pref.matrix=pm, block.sizes=rep((.5*N)/clusters,clusters))
    l <-layout_nicely(g)
    g=g %du% g
    l1=l
    l1[,1]=l1[,1]
    l1[,2]=l1[,2]+36 
    l=rbind(l,l1)
    adj <- get.adjacency(g)
  }
  
  return(list(g,l,adj,gg))
}

plot.graph = function(adj,l,opin,inform,atten,CD,title="",shape='circle',c.title=1.5)
{
  vertex.c=diverging_hsv(256,s=2,v=1)[as.numeric(cut(opin,seq(min.o,max.o,length=255)))]
  vertex.border.c=rep('red',length(opin))
  vertex.border.c[inform<0]='blue'
  vertex.border.c[opin*inform>=0 | CD>0]='black'
  border.w=rep(1,length(opin))
  border.w[opin*inform>=0| CD>0]=.5
  par(mar=c(0,0,0,0))
  if(scenario %in% c(2,3)) par(mar=c(1,2,3,4))
  plot(range(l[,1]),range(l[,2]),type = "n", bty = "n", axes = FALSE)
  
  qgraph(adj, layout = l, vsize = 0.4 * (2+3*atten), labels = FALSE, color=vertex.c,
         shape = shape, esize = 0.5, edge.color = rgb(0.5,0.5,0.5, 0.5),border.color=vertex.border.c,
         border.width = border.w, plot = F, rescale = F,title=title,title.cex=c.title)
  par(mar=c(5,4,4,2));
}

plot.histo= function(x,min,max,xlab='')
{
  hist.c=rev(diverging_hsv(40,s=2,v=1))
  par(mar=c(3,0,0,0));par(mgp = c(1, 0, 0))
  hist(x,main='',xlab=xlab,ylab='',breaks=seq(min,max,le=40),col=hist.c,cex.lab=1.1,axes=F,freq=F);
  axis(1,at=0,labels = '')
  par(mar=c(5,4,4,2));par(mgp = c(3, 1, 0))
}

##############

scenario=4  # set scenario

# scenario 1: removed from manuscript
# scenario 2 : figure 4
# scenario 3 : figure 5
# scenario 31 : figure 6 (warning: takes hours)
# scenario 4 : figure 7
# scenario 41 : figure 8 (warning: takes hours)

pdfplot=T # if TRUE to plot to pdf files used for manuscript figures
plottype='tiff'

PNG=F # save png for gif
if(PNG) unlink(paste0("figures/pngplots_",scenario,"/*"))

#layout for plots
if(scenario %in% c(1)) layout.m=layout.basis=layout(1)
if(scenario %in% c(2,3,31)) layout.m=layout.basis=matrix(c(rep(1,8),0,2,2,3,3,4,4,0),2,8,byrow=T)
if(scenario %in% c(4,41)) layout.m=layout.basis=matrix(c(rep(1,16),2,3,4,5),4,5)
if(scenario %in% c(1)) heights=heights.basis=1
if(scenario %in% c(2,3,31)) heights=heights.basis=c(4,1)
if(scenario %in% c(4,41)) heights=heights.basis=c(1,1)

if(scenario>30) simulation=T else simulation=F # simulation runbs oer conbination of parameter settings

if(simulation) {plotting=F} else{plotting=T;sim_values=sim_values2=1;sim_values3='SBM'}  

#simulation settings
if(scenario==31) {sim_values=rep(seq(0,.4,by=.05),each=5);sim_var=c('d_A');
sim_values2=rep(seq(.0,1,by=.1),each=5);sim_var2=c('r_min')
sim_values3=c(1,2);sim_var3=c('network')}
if(scenario==41) {sim_values=rep(seq(0,.002,by=.0005),each=5);sim_var=c('p(pertubation)')
sim_values2=rep(seq(0,.5,by=.1),each=5);sim_var2=c('t_d')
sim_values3=c(1,2);sim_var3=c('network')}


datasim=matrix(NA,length(sim_values)*length(sim_values2)*length(sim_values3),7) # data storage simulations
sim_i=0 #counter in simulation loop

## simulation loop ; not used in scenario 1,2,3,4
for(sim_value in sim_values) 
{
  for(sim_value2 in sim_values2)
  {
    for(sim_value3 in sim_values3)
    {
      sim_i=sim_i+1
      print(sim_i)
      
      if(scenario==1) #the cusp of cusp; removed from manuscript
      {
        Ni= 5000 # max iterations
        if(pdfplot)
        { 
          set.seed(1);
          pdfname=paste('figures/cuspofcusps',scenario,'.pdf',collapse='',sep='')
          pdf(pdfname,h=8,w=11);
          plot_iteration=c(1,(1/10)*Ni,Ni)
          layout.m=matrix(1:3,3,1)
        }
        
        n1=200;n2=200  # dimensions of lattice
        N=n1*n2
        lattice=T
        g=make_lattice(c(n1, n2, 1),nei=1);l = layout_on_grid(g); # lattice
        
        sd_noise_information=.0
        info_update=F # if false information is not updated
        persuasion=1
        r_min=0
        
        s_O=.00001
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0 # increase in attention when interacting
        decay_attention=delta_attention/(attention_star*(N^2))
        deffuant_c=Inf
        
        # initial values for A and I
        k_init=.2
        i_init_min=.5
        i_init_max=1
        information=rep(seq(-1*k_init,k_init,length=n1),each=n2)
        attention=rep(seq(i_init_min,i_init_max,length=n1),n2)
        opinion=rnorm(N,0,.3)
      }
      
      if(scenario==2)
      {
        Ni=15000 # max iterations
        if(pdfplot)
        { 
          set.seed(1);
          if(plottype=='pdf'){
            pdfname='figures/figure4.pdf'
            pdf(pdfname,h=8,w=11);
          }
          if(plottype=='tiff'){
            pdfname='figures/figure4.tiff'
            dpi=300
            tiff(pdfname,h=8,w=11,res=dpi,units = "in");
          }
      
          layout.m=cbind(layout.basis,layout.basis+4*(layout.basis>0));
          layout.m=rbind(layout.m,layout.m+8*(layout.m>0))
          heights=c(3,1,3,1)
          plot_iteration=c(1,(1/3)*Ni,(2/3)*Ni,Ni)
        }
        
        N=400
        network=make_network(model='WS',N=N,clusters=10,p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        sd_noise_information=.0005
        info_update=T
        persuasion=1
        r_min=0
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5
        delta_attention=0
        decay_attention=delta_attention/(attention_star*(N^2))
        deffuant_c=Inf
        
        # all slighty positive low attention
        information = rnorm(N, 0, .3)
        attention = runif(N, 1, 1)
        opinion=rnorm(N,0,.2)
        
      }
      
      if(scenario==3)
      {
        Ni= 15000 #25000 # max iterations
        if(pdfplot)  
        {
          set.seed(12);
          if(plottype=='pdf'){
            pdfname='figures/figure5.pdf'
            pdf(pdfname,paper='a4r',h=8,w=12);
          }
          if(plottype=='tiff'){
            pdfname='figures/figure5.tiff'
            dpi=300
            tiff(pdfname,h=8,w=11,res=dpi,units = "in");
          }
          
          layout.m=cbind(layout.basis,layout.basis+4*(layout.basis>0));
          layout.m=rbind(layout.m,layout.m+8*(layout.m>0))
          layout.m[3:4,9:16]=13
          heights=c(3,1,3,1)
          plot_iteration=c(600,3000,Ni)
        }
        
        N=400
        social_network='SBM'
       # social_network='lattice'
        
        network=make_network(social_network,N=N,clusters=10,p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        info_update=T
        sd_noise_information=.005
        persuasion=2
        r_min=0.1
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0.1
        decay_attention=delta_attention/(attention_star*(N^2))
        deffuant_c=Inf
        
        # all slighty positive low attention
        information = rnorm(N, .1, 0)
        attention = runif(N, .0, .0)
        opinion=rnorm(N,0,.01)
        for( i in 1:500) opinion=stoch_cusp(N,opinion,attention+min_attention,
                                            information,s_O,maxwell_convention)
        
        # except for some negative high attention persons (activists)
        m=rep(0,N)
        m[seq(1,N,N/3)[-1]-1]=1
        information_activists=-.5
        attention_activists=1
        opinion_activists=-.5
        activism_on_at_iteration=300
        
      }
      
      if(scenario %in% c(31))
      {
        N=400
        network_sim=c('SBM','lattice')
        network=make_network(network_sim[sim_value3],N=N,clusters=10,
                             p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        info_update=T
        Ni= 25000# 5000 # max iterations
        sd_noise_information=.0005
        persuasion=2
        r_min=0
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0.1
        decay_attention=delta_attention/(attention_star*(N^2))
        deffuant_c=Inf
        
        delta_attention= sim_value
        r_min= sim_value2
        
        # all slighty positive low attention
        information = rnorm(N, .1, 0)
        attention = runif(N, .0, .0)
        opinion=rnorm(N,0,.01)
        for( i in 1:500) opinion=stoch_cusp(N,opinion,attention+min_attention,
                                            information,s_O,maxwell_convention)
        # except for some negative high attention persons (activists)
        m=rep(0,N)
        m[seq(1,N,N/3)[-1]-1]=1
        information_activists=-.5
        attention_activists=1
        opinion_activists=-.5
        activism_on_at_iteration=300    
        
        plot_iteration=c(1,Ni)
        plot_iteration=c(1,Ni)=NA
      }
      
      if(scenario==4) # Meat eating vegies
      {
        Ni= 30000 #150000 # max iterations
        if(pdfplot) 
        {
          set.seed(3);
          if(plottype=='pdf'){
           pdfname='figures/figure7.pdf'
           pdf(pdfname,h=8,w=11,paper='a4r');
          }
          if(plottype=='tiff'){
            pdfname='figures/figure7.tiff'
            dpi=300
            tiff(pdfname,h=8,w=11,res=dpi,units = "in");
          }
          layout.m=layout.m=cbind(layout.m,layout.m+5)
          plot_iteration=c(1,Ni) 
        }
        
        N=800
        network=make_network('SBM_4',N=N,clusters=10,p_within=.2,p_between=.001)
        g=network[[1]];l=network[[2]];adj=network[[3]];gg=network[[4]]
        
        lattice=F
        info_update=T
        
        sd_noise_information=.0001
        persuasion=2
        r_min=0
        init_opinion_random = FALSE
        s_O=.01
        attention_star=1
        min_attention=-0.5 # to include continuous change in O as function of K
        delta_attention=0.01
        deffuant_c=0.2 # only interactions when |o1-o2| < .2
        p_pertubation=0.0005 # probability of inserting meat-eating vegatarians
        maxwell_convention=F
        
        information=sample(c(.1,-.4),N,T,prob=c(.8,.2))
        attention=rep(1,N)
        attention[information>0]=0.1
        opinion=rnorm(N,0,.01)
        
        information[(.5*N+1):N]=information[1:(.5*N)]
        attention[(.5*N+1):N]=attention[1:(.5*N)]
        opinion[(.5*N+1):N]=opinion[1:(.5*N)]
        for( i in 1:500) opinion=stoch_cusp(N,opinion,attention+min_attention,
                                            information,s_O,maxwell_convention)
      }
      
      if(scenario==41) # Meat eating veggies
      {
        N=400
        network=make_network('SBM',N=N,clusters=10,p_within=.2,p_between=.001)
        # network=make_network('lattice',N=N,clusters=10,p_within=.2,p_between=.001)
        g=network[[1]];l=network[[2]];adj=network[[3]];gg=network[[4]]
        
        lattice=F
        info_update=T
        Ni= 30000 #150000 # max iterations
        sd_noise_information=.0001
        persuasion=2
        r_min=0
        
        s_O=.01
        min_attention=-.5 # to include continuous change in O as function of K
        attention_star=1
        delta_attention=.01 # increase in attention when interacting
        deffuant_c=0.2 # only interactions when |o1-o2| < .2
        p_pertubation=0.001 # probability of inserting meat-eating vegatarians
        maxwell_convention=F
        
        p_pertubation= sim_value
        deffuant_c= sim_value2
        
        information=sample(c(.1,-.4),N,T,prob=c(.8,.2))
        attention=rep(1,N)
        attention[information>0]=0.1
        opinion=rnorm(N,0,.01)
        for( i in 1:500) opinion=stoch_cusp(N,opinion,attention+min_attention,
                                            information,s_O,maxwell_convention)
      }
      
      # plot and write settings
      if(!pdfplot & !simulation & !PNG){layout(matrix(1:1,1,1));plot_iteration=seq(0,Ni,2000)}
      if(PNG &scenario %in% c(1,2)) plot_iteration=c(seq(1,Ni,500))
      if(PNG &scenario %in% c(3,4)) plot_iteration=c(seq(1,1000,200),
                                                       seq(1000,.5*Ni,500),seq(.5*Ni,Ni,1000))
      
      
      #remember initial values
      opinion_init=opinion
      attention_init=attention
      information_init=information
      
      # make neighbor matrix for fast sampling neighbors in loop
      max_neighb=max(degree(g))
      m_neigh=matrix(NA,N,max_neighb)
      for(i in 1:N)
      {
        nb=as.numeric(neighbors(g, i, mode = c("all"))[[]])
        lengt_nb=length(nb)
        nb=c(nb,rep(NA,max_neighb-lengt_nb))
        m_neigh[i,]=nb
      }
      
      # collect data in
          freq_data_store=1
          if(!simulation) data=matrix(NA,Ni/freq_data_store,25) else data=matrix(NA,20,25) 
       
      layout(layout.m,heights=heights)
      
      ######################################
      #############  Main loop #############
      ######################################
      
      for(iteration in 1:Ni)
      {
        if(iteration%%1000==0 & !simulation) print(iteration)
        
        # introduce activista in scenario 3 and 31 after activism_on_at_iteration interations
        if(scenario %in% c(3,31)) if(iteration == activism_on_at_iteration) 
        {information[m == 1] = information_activists
        attention[m == 1] = attention_activists
        opinion[m == 1] = opinion_activists
        }
        
        #sample agent
        if(max(attention)==0) agent=0 else agent=sample(1:N,1,prob=attention)
        
        #sample partner
        if(agent!=0) 
        {
          neighb=m_neigh[agent,][!is.na(m_neigh[agent,])]
          if(length(neighb>0))  partner=sample(neighb,1) else partner=0
        } else partner=0
        
        if(partner!=0&agent!=0 ) 
        {
          I1=information[agent];A1=attention[agent];  O1=opinion[agent];
          I2=information[partner];A2=attention[partner];O2=opinion[partner]
        }
        
        if(iteration>1) # exclude first iteration from updating to get a plot of initial state
        {
          ### information update for interacting agents
          if(partner!=0&agent!=0 ) 
          {
            
            if( abs(O1-O2) < deffuant_c) # bounded confidence check
            {
              information[agent]=information_update(I1,I2,A1,A2,persuasion,r_min,info_update)
              information[partner]=information_update(I2,I1,A2,A1,persuasion,r_min,info_update)
            }
            
            ### attention update for interacting agents
            if( abs(O1-O2) < deffuant_c )
            {
              attention[agent]=attention[agent]+delta_attention*(2*attention_star-attention[agent]) 
              attention[partner]=attention[partner]+delta_attention*(2*attention_star-attention[partner])
            }
          }
          
          ### attention decay for all agents
          information=information+rnorm(N,0,sd_noise_information)
          attention=attention-2*delta_attention*attention/N   # correction 2 times if interaction
          
          # scenario 2 schrinking
          if(scenario==2) if(iteration>(1/3)*Ni) {information=.999*information;
                                          sd_noise_information=0}  # schrinking I in scenario 1
          if(scenario==2) if(iteration>(2/3)*Ni) {attention=.999*attention;
                                          delta_attention=0} # schrinking A in scenario 1
          
          ###  behavior update
          opinion=stoch_cusp(N,opinion,attention+min_attention,information,s_O,maxwell_convention)
          
          # introducing meat eating vegatations in scenario 4
          if(scenario==4) 
          {
            perturb_vector=sample(c(F,T),N,T,prob=c(1-p_pertubation,p_pertubation))
            opinion[opinion<0 & perturb_vector & rep(c(T,F),each=.5*N)] = -1*opinion[opinion<0 & 
                                                            perturb_vector & rep(c(T,F),each=.5*N)]
          }
          
          if(scenario==41) # introducing meat eating vegatations in scenario 4
          {
            perturb_vector=sample(c(F,T),N,T,prob=c(1-p_pertubation,p_pertubation))
            opinion[opinion<0 & perturb_vector] = -1*opinion[opinion<0 & perturb_vector]
          }
        }
        
        # reports
        if(partner!=0&agent!=0 & iteration %% freq_data_store==0 & !simulation)  
        {
          CD=27*information^2-4*(attention+min_attention)^3 # Cardan's discriminant
          CD_p=sum(CD<0)/N
          freq_pos_opinion=sum(opinion>0)/N  # p(O>0)
          assortativity_g=assortativity(g,opinion)
          ambivalence_r=cor(opinion,information)
          dip=dip(opinion)
          
          data[iteration/freq_data_store,]=c(agent,partner,O1,O2,I1,I2,A1,A2,opinion[agent],opinion[partner],
                                 information[agent], information[partner],attention[agent],
                                 attention[partner],
                                 CD_p,ambivalence_r,freq_pos_opinion,
                                 assortativity_g,dip,mean(opinion),sd(opinion),mean(information),
                                 sd(information),mean(attention),sd(attention))
        }
        
        if(partner!=0&agent!=0 & iteration > (Ni-21) & simulation)  
        {
          CD=27*information^2-4*(attention+min_attention)^3
          CD_p=sum(CD<0)/N
          freq_pos_opinion=sum(opinion>0)/N
          assortativity_g=assortativity(g,opinion)
          ambivalence_r=cor(opinion,information)
          dip=dip(opinion)
          
          data[iteration-(Ni-20),]=c(agent,partner,O1,O2,I1,I2,A1,A2,opinion[agent],
                                     opinion[partner],information[agent],
                                     information[partner],attention[agent],attention[partner],
                                     CD_p,ambivalence_r,freq_pos_opinion,
                                     assortativity_g,dip,mean(opinion),sd(opinion),mean(information),
                                     sd(information),mean(attention),sd(attention))
        }
        
        ###### Plotting
        if(plotting)
        {
          if(PNG) {
            png(paste0("figures/pngplots_",scenario,"/",scenario,"_",100000+iteration,".png"),
                            width=8,height=5,units="in",res=200)
            layout(layout.m,heights=heights)
          }
          
          if(iteration %in% plot_iteration)
          {
            min.o=-1.5;max.o=1.5;min.i=-1;max.i=1;min.a=0;max.a=2
            opin=sapply(-opinion, function(y) min(max(y,min.o),max.o))
            inform=sapply(-information, function(y) min(max(y,min.i),max.i))
            atten=attention;atten[atten>2]=2
            CD=27*information^2-4*(attention+min_attention)^3
            assortativity_g=assortativity(g,opinion)
            
            shape=rep("circle",N)
            if(pdfplot)c.title=1.5  else c.title=1.5
            title=sub_H=sub_A=""
            if(scenario==2 )
            {
              h_d=hartigan_d(opinion)
              if (iteration >= 1) title = "1: Initial state"
              if (iteration >= (1/3)*Ni) title = "2: A = 1, I varies"
              if (iteration >= (2/3)*Ni) title= "3: A = 1, I = 0"
              if (iteration == Ni) title = "4: A = 0, I = 0"
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity_g,2))
              sub=paste(sub_H,sub_A)
            }
            
            if(scenario==3 )
            {
              h_d=hartigan_d(opinion)
              title=paste('p(O>0) =',round(sum(opinion>0)/N,2),'\nsd(O) =',round(sd(opinion),2))
              shape[m==1]='star'
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity_g,2))
              sub=paste(title,sub_H,sub_A)
              title=""
              
            }
            
            if(scenario %in% c(2,3)) 
            {
              plot.graph(adj,l,opin,inform,atten,CD,title,shape,c.title)
              if(scenario==3) text(9,-16,sub,cex=1.2,pos=4,xpd=NA)
              
              plot.histo(-opin,min.o,max.o,'O')
              plot.histo(-inform,min.i,max.i,'I')
              plot.histo(atten, min.a,max.a,'A')
            }
            if(scenario==2) text(.9,32,sub,cex=1.2,pos=4,xpd=NA)
            
            if(scenario==4){
              
              h_d=hartigan_d(opinion[(1+(N/2)):N])
              title=paste('p(O>0) =',round(sum(opinion[(1+(N/2)):N]>0)/(.5*N),2),
                          '\nsd(O) =',round(sd(opinion[(1+(N/2)):N]),2))
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity( gg,opinion[(1+(N/2)):N]),2))
              sub1=paste(title,sub_H,sub_A)
              
              h_d=hartigan_d(opinion[1:(N/2)])
              sub=paste('p(O>0) =',round(sum(opinion[1:(N/2)]>0)/(.5*N),2),
                        '\nsd(O) =',round(sd(opinion[1:(N/2)]),2))
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity(gg,opinion[1:(N/2)]),2))
              sub=paste(sub,sub_H,sub_A)
              
              plot.graph(adj,l,opin,inform,atten,CD,title="",shape,c.title)
              par(mar=c(0,0,0,0))
              plot(c(0,1),c(0,1),type='n',axes=FALSE,ann=FALSE)
              text(-.22,.3,sub1,cex=1.2,pos=4,xpd=NA)
              plot.histo(-opin[(1+(N/2)):N],min.o,max.o,'O')
              par(mar=c(0,0,0,0))
              plot(c(0,1),c(0,1),type='n',axes=FALSE,ann=FALSE)
              text(-.22,.3,sub,cex=1.2,pos=4,xpd=NA)
              plot.histo(-opin[1:(N/2)],min.o,max.o,'O')
            }  
            
            if(scenario==1 )
            {
              opin_m=matrix(opin,n1,n2)
              image(opin_m,axes=F,col = diverging_hsv(256,s=1.6,v=1),
                    breaks=seq(-1.5,1.5,le=257),xlab='A', ylab='I(init)')
              axis(1,seq(0,1,length=6),seq(i_init_min,i_init_max,by=.1))
              axis(2,seq(0,1,by=.25),seq(-1*k_init,k_init,length=5))
              scale_factor= (i_init_max-i_init_min)/(2*k_init)
              bf2=function(b) .5+scale_factor*(2 * (b+min_attention+i_init_min)^(3/2))/(3 *3^(1/2))
              plot(bf2,0,1,ylab='b',xlab='t',add=T)
              bf2=function(b) .5-1*scale_factor*(2 * (b+min_attention+i_init_min)^(3/2))/(3 *3^(1/2))
              plot(bf2,0,1,ylab='b',add=T)
            }
          }      
          if(PNG) dev.off()
        }
      }
      ####################
      ########end#########
      ####################
      
      # cusp fig plot in scenario 3
      if(scenario==3 & pdfplot) 
      {cuspfig <- readPNG('figures/figure7.png')
      plot(c(0,1),c(0,1), type='n', main="", xlab="", ylab="",axes=F)
      par(mar=c(0,0,0,0))
      rasterImage(cuspfig, .05,-.05, .9,1.05)
      }
      
      if(pdfplot) dev.off()
      
      datad=as.data.frame(data)
      names(datad)=c('agent','partner','o1o','o2o','i1o','i2o','a1o','a2o','o1n',
                     'o2n','i1n','i2n','a1n','a2n','CD',
                     'ambivalence_r','freq_pos_opinion','assortativity','dip','mean_o','sd_o',
                     'mean_i','sd_i','mean_a','sd_a')
      
      datasim[sim_i,]=c(sim_value,sim_value2,sim_value3,mean(datad$freq_pos_opinion),mean(datad$sd_o),
                        mean(datad$dip),mean(datad$assortativity))
    }
  }
} # end of simulation loops

if(scenario==31)
{
#  save(datasim,file='figures/datasim_sc31')
  load('figures/datasim_sc31_org')
  sim_var=c('d_A');sim_var2=c('r_min'); sim_var3=c('network')
  mean_sim=as.matrix(aggregate(datasim[,4:7],list(datasim[,1],datasim[,2],datasim[,3]),mean,na.rm=T))
  colnames(mean_sim)=c(sim_var,sim_var2,sim_var3,c('p(O>0)','SD(O)','Hartigan D','Assortativity'))
  mean_sim2=as_tibble(mean_sim)
  mean_sim2$r_min=as.factor(mean_sim2$r_min)
  mean_sim2$network[mean_sim2$network==1]='Stoch. block network'
  mean_sim2$network[mean_sim2$network==2]='Lattice network'
  mean_sim2$network=factor(mean_sim2$network, levels=c('Stoch. block network','Lattice network'))
  mean_sim3 = mean_sim2 %>%
    gather('p(O>0)','Hartigan D', key = measure, value = value)
  mean_sim3$measure=factor(mean_sim3$measure, levels=c('p(O>0)','Hartigan D'))
  
  if(plottype=='pdf')pdf('figures/figure6.pdf')
  if(plottype=='tiff'){
  dpi=300
  tiff('figures/figure6.tiff',h=6,w=6,res=dpi,units = "in")
  }
  print(ggplot(data = mean_sim3) + 
          geom_point(mapping = aes(x = d_A, y = value,color = r_min),size=1.2)+
          geom_line(mapping = aes(x = d_A, y = value,color = r_min),size=.5)+
          facet_grid(measure~network,scales="free",switch='y')+
          scale_colour_hue(h = c(90, 360)) + theme_grey()+ xlab(expression('d'['A']))+ylab('')+
          guides(color = guide_legend(reverse = TRUE))+
          theme(strip.text.x = element_text( face = "bold", size=12))+
          theme(strip.text.y = element_text( size=11),strip.background = element_blank(),strip.placement = "outside")+
          labs(color = expression('r'['min']))+
          theme(panel.spacing = unit(1, "cm")))
  
  dev.off()
}

if(scenario==41)
{
  #save(datasim,file='figures/datasim_sc41')
  load('figures/datasim_sc41_org')
  sim_var=c('p(pertubation)');sim_var2=c('t_d');sim_var3=c('network')
  mean_sim=as.matrix(aggregate(datasim[,4:7],list(datasim[,1],datasim[,2],datasim[,3]),mean,na.rm=T))
  colnames(mean_sim)=c('p_perturbation',sim_var2,sim_var3,c('p(O>0)','SD(O)','Hartigan D','Assortativity'))
  mean_sim2=as_tibble(mean_sim)
  mean_sim2$t_d=as.factor(mean_sim2$t_d)
  mean_sim2$network[mean_sim2$network==1]='Stoch. block network'
  mean_sim2$network[mean_sim2$network==2]='Lattice network'
  mean_sim2$network=factor(mean_sim2$network, levels=c('Stoch. block network','Lattice network'))
  
  mean_sim3 = mean_sim2 %>%
    gather('p(O>0)','Hartigan D', key = measure, value = value)
  mean_sim3$measure=factor(mean_sim3$measure, levels=c('p(O>0)','Hartigan D'))

  if(plottype=='pdf') pdf('figures/figure8.pdf')
  if(plottype=='tiff'){
    dpi=300
    tiff('figures/figure8.tiff',units="in",h=6,w=6,res=dpi)
  }
  
  scaleFUN <- function(x) sprintf("%.4f", x)

  print(ggplot(data = mean_sim3,aes(x = p_perturbation, y = value,color = t_d)) + 
    geom_point(size=1.2)+
    geom_line(size=.5)+
    facet_grid(measure~network,scales="free",switch='y')+
    scale_colour_hue(h = c(90, 360)) + theme_grey()+ xlab('p(pertubation)')+ylab('')+
    guides(color = guide_legend(reverse = TRUE))+
    theme(strip.text.x = element_text( face = "bold",size=12))+
    theme(strip.text.y = element_text( size=11),strip.background = element_blank(),strip.placement = "outside")+
    labs(color = expression('t'['O']))+
    theme(panel.spacing = unit(1, "cm")))
  
  dev.off()
}

layout(1:3)
hist(opinion,20,col='grey');hist(information,20,col='grey');hist(attention,20,col='grey')

layout(1)
matplot(datad[,c('CD','freq_pos_opinion','assortativity','ambivalence_r')],
        type='l',bty='n',xlab='x 100',ylab="value",axes=F,col=1:4,lty=1:4)
legend('topright',col=1:7,lty=1:7,legend=c('in bifurcation set','ambivalence',
        'freq_pos_opinion','assortativity','ambivalence_r'),bty='n')
axis(2)
axis(1)

# make gif if PGN = TRUE
if(PNG)
{
  library(dplyr)
  library(purrr) 
  library(magick)
  
  list.files(path=paste0("figures/pngplots_",scenario,"/"), pattern = '*.png', full.names = TRUE) %>% 
    image_read() %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=2,loop=1) %>% # animates, can opt for number of loops
    image_write(paste0("figures/Anim_",scenario,".gif")) # write to current dir
  
}


