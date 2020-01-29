rm(list=ls(all=TRUE))
#e0: parâmetro da priori dirichilet
#c0,C0: parâmetro da priori IG
#b0,B0: parâmetro da priori N
#p0: priori S
#e0: priori epsilon
require(circular)



#parâmetros dos dados gerados
T=200 #tamanho da série temporal
r=100 #número de réplicas
mu1=pi/2
mu2=3*pi/2
sigma21=0.1
sigma22=0.2

#inicializando as dimensões
Sreal=matrix(0,T,r)
dad=matrix(0,T,r)
e=array(0,dim=c(2,2,r))

#Gerando os dados
for(j in 1:r){
  Sreal[1,j]=rbinom(1,1,0.3)
  Sreal[1,j]=Sreal[1,j]+1
  if(Sreal[1,j]==1){
    dad[1,j]=rwrappednormal(1,mu1,sd=sqrt(sigma21))
  }
  if(Sreal[1,j]==2){
    dad[1,j]=rwrappednormal(1,mu2,sd=sqrt(sigma22))
  }
  e[1,1,j]=0.8
  e[1,2,j]=1-e[1,1,j]
  e[2,2,j]=0.8
  e[2,1,j]=1-e[2,2,j]
  
  for(i in 2:T){
    u=runif(1,0,1)
    if(Sreal[i-1,j]==1){
      if(u<0.8){
        Sreal[i,j]=1
        dad[i,j]=rwrappednormal(1,mu1,sd=sqrt(sigma21))
      }
      else{
        Sreal[i,j]=2
        dad[i,j]=rwrappednormal(1,mu2,sd=sqrt(sigma22))
      }
    }
    if(Sreal[i-1,j]==2){
      if(u<0.2){
        Sreal[i,j]=1
        dad[i,j]=rwrappednormal(1,mu1,sd=sqrt(sigma21))
      }
      else{
        Sreal[i,j]=2
        dad[i,j]=rwrappednormal(1,mu2,sd=sqrt(sigma22))
      }
    }
  }
}





##############################################################################################
#                                                                                            #
#                                                                                            #
#                                          FUNÇÃO                                            #
#                                                                                            #
#                                                                                            #
##############################################################################################



mist=function(dad=dad,e0,c0,C0,b0,B0,k,muini,sigma2ini,Sini,p0,n){
  #Definindo a dimensão dos parâmetros
  T=length(dad)
  e=array(0,dim=c(k,k,n))
  p=matrix(0,k,T)
  probfilt=matrix(0,T,k)
  probcond=matrix(0,T,k)
  c=matrix(data=0,k,n)
  C=c
  sigma2=c
  S=matrix(data=0,T,n)
  desvq=matrix(data=0,k,n) #soma dos desvios (y-mu)^2 em cada grupo 
  N=c
  b=c
  B=c
  media=c
  w=c
  aly=S
  mu=c
  sigma2=c
  auxi1=matrix(0,2,2)
  auxi2=matrix(0,2,2)
  y=matrix(0,T, n)
  K=matrix(0,T, n)
  prob=NA
  auxiprob=NA
  mu[,1]=muini
  sigma2[,1]=sigma2ini
  y[,1]=dad
  S[,1]=Sini
  
  
  for(i in 2:n){
    #########################################################INÍCIO DA PARTE WRAPPED
    ##############estimação de K
    for (q in 1:T){
      prob=auxiprob
      #########Estimação de K para o grupo 1
      if(S[q,i-1]==1){
        men=1+trunc(3*sqrt(sigma2[1,i-1])/(2*pi))
        for (j in 1:(2*men+1)){
          mauxik=j-(men+1)
          pt=(dad[q]+2*pi*mauxik-mu[1,i-1])/sqrt(sigma2[1,i-1])
          prob[j]=dnorm(pt,0,1)
        }
        prob=prob/(sum(prob))
        p=runif(1,0,1)
        if(0<p & (p<=prob[1])){
          K[q,i]=-men
        }
        con=prob[1]
        for (j in 2:(2*men)){
          mauxik=j-(men+1)
          pauxi=prob[j]+con
          if(con<p & p<=pauxi){ 
            K[q,i]=mauxik 
          }
          con=con+prob[j]
        }
        if(con<p & p<=1){ 
          K[q,i]=men
        }
      }
      ########Estimação de k para o grupo 2
      if(S[q,i-1]==2){
        men=1+trunc(3*sqrt(sigma2[2,i-1])/(2*pi))
        for (j in 1:(2*men+1)){
          mauxik=j-(men+1)
          pt=(dad[q]+2*pi*mauxik-mu[2,i-1])/sqrt(sigma2[2,i-1])
          prob[j]=dnorm(pt,0,1)
        }
        prob=prob/(sum(prob))
        p=runif(1,0,1)
        if(0<p & (p<=prob[1])){
          K[q,i]=-men
        }
        con=prob[1]
        for (j in 2:(2*men)){
          mauxik=j-(men+1)
          pauxi=prob[j]+con
          if(con<p & p<=pauxi){ 
            K[q,i]=mauxik 
          }
          con=con+prob[j]
        }
        if(con<p & p<=1){ 
          K[q,i]=men
        }
      }
      
    }
    y[,i]=K[,i]*2*pi+dad
    #########################################################FIM DA PARTE WRAPPED
    
    #Calculando Nt e N
    Nt=matrix(0,k,k)
    for(j in 1:(T-1)){
      for(o in 1:k){
        for(m in 1:k){
          if(S[j,(i-1)]==o){
            if(S[(j+1),(i-1)]==m)
              Nt[o,m]=Nt[o,m]+1
          }
        }
      }
    }
    for (j in 1:k){
      N[j,i]=sum(S[,(i-1)]==j)
    }
    
    #Amostrando os pesos
    
    e11=e0[1,1]+Nt[1,1]
    e12=e0[1,2]+Nt[1,2]
    e[1,1,i]=rbeta(1,e11,e12)
    e[1,2,i]=1-e[1,1,i]
    e21=e0[2,1]+Nt[2,1]
    e22=e0[2,2]+Nt[2,2]
    e[2,2,i]=rbeta(1,e21,e22)
    e[2,1,i]=1-e[2,2,i]
    
    #Amostrando as variâncias
    
    com=1
    fin=0
    for (j in 1:k){
      ind=0
      yax=0
      nulo=sum(which(S[,(i-1)]==j))
      if(nulo==0){
        c[j,i]=c0+0.5*N[j,i]
        C[j,i]=C0
        sigma2[j,i]=1/rgamma(1,c[j,i],scale=1/C[j,i])
      }
      else{
        ind=which(S[,(i-1)]==j)
        yax=y[ind,i]
        fin=com+N[j,i]-1
        aly[com:fin,i]=yax
        media[j,i]=mean(yax)
        soma=NA
        for(o in 1:N[j,i]){
          soma[o]=(yax[o]-mu[j,i-1])^2
        }
        desvq[j,i]=sum(soma)
        c[j,i]=c0+0.5*N[j,i]
        C[j,i]=C0+0.5*desvq[j,i]
        sigma2[j,i]=1/rgamma(1,c[j,i],scale=1/C[j,i])
        com=fin+1
      }
    }
    
    #Amostrando os muk's
    for (j in 1:k){
      B[j,i]=1/(1/B0+N[j,i]/sigma2[j,i])
      b[j,i]=B[j,i]*(N[j,i]*media[j,i]/sigma2[j,i]+b0/B0)
      mu[j,i]=rnorm(1,b[j,i],B[j,i])
    }
    
    #Amostrando as alocações
    
    predS=matrix(0,T,k)
    somaalo=matrix(0,k,k)
    for (l in 1:k){
      for (o in 1:k){
        somaalo[o,l]=e[o,l,i]*p0[o]
      }
      predS[1,l]=sum(somaalo[,l])
    }
    for(l in 1:k){
      probfilt[1,l]=dnorm(y[1,i],mu[l,i],sigma2[l,i])*predS[1,l]
    }
    probfilt[1,]=probfilt[1,]/sum(probfilt[1,])
    
    for(r in 2:T){
      
      #one-step ahead prediction of St
      somaalo=matrix(0,k,k)
      for (l in 1:k){
        for (o in 1:k){
          somaalo[o,l]=e[o,l,i]*probfilt[(r-1),o]
        }
        predS[r,l]=sum(somaalo[,l])
      }
      
      #filtering for St
      
      for(l in 1:k){
        probfilt[r,l]=dnorm(y[r,i],mu[l,i],sigma2[l,i])*predS[r,l]
      }
      probfilt[r,]=probfilt[r,]/sum(probfilt[r,])
    }
    p=runif(1,0,1)
    con=0
    for (z in 1:(k-1)){
      pauxi=probfilt[z,l]+con
      if(con<p & p<=pauxi){ 
        S[T,i]=z 
      }
      con=con+probfilt[z,l]
    }
    if(con<p & p<=1){ 
      S[T,i]=k
    }
    
    
    #sampling from the conditional distributions
    
    
    probcond[T,]=probfilt[T,]
    for(z in 1:(T-1)){
      inv=T-z
      for(u in 1:k){
        for (o in 1:k){
          auxi1[u,o]=e[u,o,i]*probfilt[inv,u]*probcond[(inv+1),o] 
          for (h in 1:k){
            auxi2[h,o]=e[h,o,i]*probfilt[inv,h] #eh possivel melhorar o cod
          }
          auxi1[u,o]=auxi1[u,o]/sum(auxi2[,o])
        }
        probcond[inv,u]=sum(auxi1[u,])
      }
      
      
      p=runif(1,0,1)
      con=0
      for (d in 1:(k-1)){
        pauxi=probcond[inv,d]+con
        if(con<p & p<=pauxi){ 
          S[inv,i]=d 
        }
        con=con+probcond[inv,d]
      }
      if(con<p & p<=1){ 
        S[inv,i]=k
      }
    }
    
    #permutação
    moeda=runif(1,0,1)
    
    if(moeda>0.5){
      ap=mu[1,i]
      mu[1,i]=mu[2,i]
      mu[2,i]=ap
      bp=sigma2[1,i]
      sigma2[1,i]=sigma2[2,i]
      sigma2[2,i]=bp
      
      #elementos da matriz de transição
      
      ap=e[1,1,i]
      bp=e[1,2,i]
      cp=e[2,1,i]
      dp=e[2,2,i]
      
      e[1,1,i]=dp
      e[1,2,i]=cp
      e[2,1,i]=bp
      e[2,2,i]=ap
      
      for(t in 1:T){
        if(S[t,i]==1){ 
          S[t,i]=2
        }
        else{
          S[t,i]=1
        }
      }
      
    }
    ########################
    
    print(i)
  }
  
  lista=list(mu,sigma2,S,K,e)
  names(lista)=c("mu", "sigma2","S","K","e")
  return(lista)
  
}



##############################################################################################
#                                                                                            #
#                                                                                            #
#                                      Rodando o MCMC                                        #
#                                                                                            #
#                                                                                            #
##############################################################################################



setwd("H:\\Raniere\\Dissertação\\Exercicio simulado\\Wrapped Normal\\Modelos de Misturas\\MMM\\Estudo Simulado")

k=2
p0=c(0.5,0.5)
e0=matrix(0,2,2) #colocar o pimeiro parametro de cada linha maior faz com que a probabilidadede permanencia seja maior
e0[1,1]=1.5
e0[1,2]=1
e0[2,1]=1.5
e0[2,2]=1
muini=c(3,3)
sigma2ini=c(1,1)
Sini=rbinom(T,1,0.5)

n=100000

t1=proc.time()
for(i in 1:100){
  cadeia=mist(dad=dad[,i],e0=e0,c0=2.0225,C0=0.153375,b0=3,B0=1,k=2,muini=muini,sigma2ini=sigma2ini,Sini=Sini, p0=p0,n=n)
  #priori inv gama com media 0.15 e var 1
  mu=cadeia$mu
  sigma2=cadeia$sigma2
  S_mcmc=cadeia$S
  K_mcmc=cadeia$K
  e_mcmc=cadeia$e
  dados=dad[,i]
  aloc=Sreal[,i]
  save(dados,aloc,mu,sigma2,S_mcmc,K_mcmc,e_mcmc,file=paste("réplica_",toString(i),".Rdata"))
}


t2=proc.time()


##############################################################################################
#                                                                                            #
#                                                                                            #
#                                         Resultados                                         #
#                                                                                            #
#                                                                                            #
##############################################################################################



####Cadeias para os mu's

par(mfrow=c(2,2))
ts.plot(cadeia$mu[1,])
abline(h=mu1,col=2)
ts.plot(cadeia$mu[2,])
abline(h=mu2,col=2)
ts.plot(cadeia$sigma2[1,])
abline(h=sigma22,col=2)
ts.plot(cadeia$sigma2[2,])
abline(h=sigma21,col=2)


mu1mcmc=NA
mu2mcmc=NA
sigma21mcmc=NA
sigma22mcmc=NA


b=40000

for(j in 1:(n-b)){
  mu1mcmc[j]=min(cadeia$mu[1,j+b],cadeia$mu[2,j+b])
  mu2mcmc[j]=max(cadeia$mu[1,j+b],cadeia$mu[2,j+b])
  
  sigma21mcmc[j]=min(cadeia$sigma2[1,j+b],cadeia$sigma2[2,j+b])
  sigma22mcmc[j]=max(cadeia$sigma2[1,j+b],cadeia$sigma2[2,j+b])
}

x11()

par(mfrow=c(2,2))
ts.plot(mu1mcmc)
abline(h=mu1,col=2)
ts.plot(mu2mcmc)
abline(h=mu2,col=2)
ts.plot(sigma21mcmc)
abline(h=sigma22,col=2)
ts.plot(sigma22mcmc)
abline(h=sigma21,col=2)