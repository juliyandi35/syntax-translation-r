n<-100000
par(mfrow = c(1,2))
g1<-function(x){ #g1 est la première densité marginale de g, M1*g étant la densité qui majore  f1 tilde.
  return(dnorm(x,mean=0,sd=2))
}

f1x<-function(x){ #f1x est la première marginale de f1.
  return(dnorm(x,mean=0,sd=2))
}

#on commence par les algorithmes de rejet des marginales de f1.

rejet.f1x<-function(n,M){  #Méthode de rejet pour simuler la 1ère marginale de f1 qui est une N(0,4)
  G<-numeric(n)
  for (j in (1:n)){
    u<-runif(1)
    y<-2*rnorm(1)
    while (u>f1x(y)/g1(y)){
      u<-runif(1)
      y<-2*rnorm(1)
    }
    G[j]<-y
  }
  return(G)
}
#ici on teste si notre rejet fonctionne bien en comparant les sorties 
#aux densités théoriques
X11<-rejet.f1x(1000,M1)
hist(X11,freq = F,breaks="FD",col="skyblue", main = "Histogramme de f1x")
t<- seq(-max(X11), max(X11), 0.01)
lines(t,dnorm(t,mean=0,sd=2),col=53,lwd=2)
legend("topright", "f1x(x)", box.lty = 0, bg = "gray90", col = "skyblue", lty = 1, lwd = 3, inset = 0.05)

# simulation par rejet de la 2ème marginale de f1.
f1y <- function(y) {
  return( exp(-(0.5)*y^2)*(abs(y)<= 1)/(sqrt(2*pi)*(pnorm(1)-pnorm(-1))))
}
g1y <- function(y) {
  return( exp(-(0.5)*y^2)/(sqrt(2*pi)))
}
rejet.basique.f1y <- function(n) { #méthode de rejet basique pour simuler f1y
  x<- numeric(n)
  for (i in 1:n) {
    x[i] <- rnorm(1)
    u <- runif(1)
    while (u >(abs(x[i])<=1)) {
      x[i] <- rnorm(1)
      u <- runif(1)
    }
  }
  return(x)
}

#ici on teste si notre rejet fonctionne bien en comparant les sorties 
#aux densités théoriques
X12<-rejet.basique.f1y(n)
hist(X12, freq = FALSE, breaks = "FD", main = "Histogramme de f1y", ylab = "Fréquences")
t <- seq(-1, 1, 0.01)
lines(t, f1y(t), col = "mediumseagreen", lwd = 3)
legend("topright", "f1y(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)

#Faisons l'algorithme de rejet pour f2.

f2 <- function(x,y) {
  return( (cos(x)^2+0.5*(sin(3*y)^2)*(cos(x)^4))*exp((-0.5)*((x^2)/4+(y^2))))
}
g <- function(x,y) { #loi normale N(0,1)*N(0,4).
  return( exp(-(0.5)*(y^2+(x^2)/4))/(4*pi))
}
coef_f2<-function(x,y) {
  return ((cos(x)^2+0.5*(sin(3*y)^2)*(cos(x)^4)))
}
rejet.basique.f2 <- function(n) { #méthode de rejet basique pour simuler f2.
  x<- numeric(n)
  y<- numeric(n)
  for (i in 1:n) {
    z1 <- 2*rnorm(1)
    z2 <- rnorm(1)
    u <- runif(1)
    while (u > 2/3*coef_f2(z1,z2) ){
      z1 <- rnorm(1, mean = 0, sd = 2)
      z2 <- rnorm(1)
      u <- runif(1)
    }
    x[i]<-z1
    y[i]<-z2
  }
  v<-rbind(x,y)
  return(v)
}

vf2<-rejet.basique.f2(n)
plot(vf2[1,], vf2[2,])


######################################################
######################################################
#Faisons l'algorithme de rejet pour f1. (On ne se restreint pas juste aux marginales)

f1 <- function(x,y) {
  return( (abs(y) <= 1)*exp((-0.5)*((x^2)/4+(y^2))))
}

rejet.basique.f1 <- function(n) { #méthode de rejet basique pour simuler f1
  x<- numeric(n)
  y<- numeric(n)
  for (i in 1:n) {
    x[i] <- 2*rnorm(1)
    y[i] <- rnorm(1)
    u <- runif(1)
    while (u > (abs(y[i]) <= 1)) {
      x[i] <- 2*rnorm(1)
      y[i] <- rnorm(1)
      u <- runif(1)
    }
  }
  v<-rbind(x,y)
  return(v)
}

#ici on teste si notre algorithme de rejet pour f1 fonctionne bien en comparant les sorties 
#aux densités théoriques
vf1<-rejet.basique.f1(n)
plot(vf1[1,], vf1[2,])

hist(vf1[2,], freq = FALSE, breaks = 25, main = "Histogramme de f1y", ylab = "Fréquences")
t <- seq(-1, 1, 0.01)
lines(t, f1y(t), col = "mediumseagreen", lwd = 3)
legend("topright", "f1y(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)
##############
hist(vf1[1,], freq = FALSE, breaks = 30, main = "Histogramme de f1x", ylab = "Fréquences")
s <- seq(-max(vf1[1,]), max(vf1[1,]), 0.01)
lines(s, f1x(s), col = "mediumseagreen", lwd = 3)
legend("topright", "f1x(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)



#PARTIE 2
######
#I
######
#
#on va maintenant essayer d'estimer par différentes méthodes de Monte Carlo la probabilité:
#P(exp(x11)+exp(X12)>5)
par(mfrow = c(3,2))

evol_estim<-function(y){ #fonction qui donne un vecteur de la somme de ses éléments divisé par le nombre d'éléments sommés
  return(cumsum(y)/(1:length(y)))
}

evol_IC<-function(y,level){ #fonction qui donne l'intervalle de confiance d'un estimateur au niveau "level"
  mu <- evol_estim(y)
  V= (cumsum((y-mu)^2)/(0:(length(y)-1)))[-1]
  mu<-mu[-1]
  q=qnorm((1+level)/2)*sqrt(V/(2:length(y)))
  Dsup=mu+q
  Dinf=mu-q
  return(data.frame(min = Dinf, max = Dsup))
}

#estimateur de Monte Carlo classique
pre_mc1<-exp(X11)+exp(X12)
pre_mc1<-pre_mc1>5
mc1<-evol_estim(pre_mc1)
mc1_IC<-evol_IC(pre_mc1, 0.95)
plot(1:n, mc1,col='mediumseagreen',type='l',ylim = c(0.23, 0.38), main = "Monte Carlo Classique", xlab="n iter", ylab ="estimation")
lines(2:n, mc1_IC$min,col='blue')
lines(2:n, mc1_IC$max,col='red')
legend("topleft", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
  col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur de Monte carlo par variable antithétique
#X11 et -X11 suivent la même loi car X11 gaussienne
#X12 et -X12 suivent la même loi (vérifié par méthode de la fonction test)
pre_mc1.ant<-exp(-X11)+exp(-X12)
pre_mc1.ant<-pre_mc1.ant>5
pre_mc1.ant<-(pre_mc1.ant+pre_mc1)/2
mc1.ant<-evol_estim(pre_mc1.ant)
mc1.ant_IC<-evol_IC(pre_mc1.ant, 0.95)
plot(1:n, mc1.ant,col='mediumseagreen',type='l',ylim = c(0.23, 0.38) ,main = "Monte Carlo Va Antithétiques", xlab="n iter", ylab ="estimation")
lines(2:n, mc1.ant_IC$min,col='blue')
lines(2:n, mc1.ant_IC$max,col='red')
legend("topleft", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur avec une variable de controle
#on va noter Z la fonction indicatrice(X11+X12>log(25)) 
#on note g la fonction indicatrice(exp(X11)+exp(X12)>5)
p<-n/50
Z<-X11+X12>log(25)
esp_Z<-sum(Z[1:p])/(p)
esp_g<-sum(pre_mc1[1:p])/(p)
cov_gZ<-sum((pre_mc1[1:p]-esp_g)*(Z[1:p]-esp_Z))
var_Z<-sum((Z[1:p]-esp_Z)*(Z[1:p]-esp_Z))
b<-cov_gZ/var_Z
Ztronc<-Z[p+1:n]
length(Ztronc)<-n-p #changement manuel de la taille de Ztronc
pre_mc1.tronc<-pre_mc1[p+1:n]
length(pre_mc1.tronc)<-n-p #changement manuel de la taille de pre_mc.tronc
pre_mc1.controle<-pre_mc1.tronc-b*(Ztronc-esp_Z)
mc1.controle<-evol_estim(pre_mc1.controle)
mc1.controle_IC<-evol_IC(pre_mc1.controle, 0.95)
m<-p+1
plot(m:n, mc1.controle,col='mediumseagreen',type='l',ylim = c(0.23, 0.38) ,main = "Monté Carlo avec controle", xlab="n iter", ylab ="estimation")
lines((m+1):n, mc1.controle_IC$min,col='blue')
lines((m+1):n, mc1.controle_IC$max,col='red')
legend("topleft", legend=c("MC.controle", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)


#nouveau controle
p<-n/50
x1c<-2*rnorm(n)
x2c<-rnorm(n)
Z1<-(x1c+x2c)>=log(4)
pre_mc.cont.1<-exp(x1c)+exp(x2c)
pre_mc.cont.1<-(pre_mc.cont.1>5)*(abs(x2c)<=1)/(pnorm(1)-pnorm(-1))
esp_Z1<-(1-pnorm(log(4), mean=0, sd=sqrt(5)))

Z1av<-Z1[1:p]
length(Z1av)<-p
g1av<-pre_mc.cont.1[1:p]
length(g1av)<-p

esp_g1<-sum(g1av)/p
cov_gZ1<-sum((g1av-esp_g1)*(Z1av-esp_Z1))
var_Z1<-sum((Z1av-esp_Z1)*(Z1av-esp_Z1))
b1<-cov_gZ1/var_Z1

Z1tronc<-Z1[p+1:n]
length(Z1tronc)<-n-p
pre_mc1.tronc1<-pre_mc.cont.1[p+1:n]
length(pre_mc1.tronc1)<-n-p

pre_mc1.controle1<-pre_mc1.tronc1-b1*(Z1tronc-esp_Z1)
mc1.controle1<-evol_estim(pre_mc1.controle1)
mc1.controle1_IC<-evol_IC(pre_mc1.controle1, 0.95)
m<-p+1
plot(m:n, mc1.controle1,col='mediumseagreen',type='l' ,ylim = c(0.23, 0.38) ,main = "Monté Carlo avec controle", xlab="n iter", ylab ="estimation")
lines((m+1):n, mc1.controle1_IC$min,col='blue')
lines((m+1):n, mc1.controle1_IC$max,col='red')
legend("topleft", legend=c("MC.controle", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)


#estimation par stratification
MC_estim.strat<-function(n, mu, sigma, L){ #car P(X E Dl)=1/L et nl=n/L donc toutes les strates font la meme taille
  Y=numeric(n)   #la probabilité d'appartenir à Dl est touijours la même
  nv <- c(n%/%L, diff((n * (1:L))%/%L))
  for (l in (1:L)){
    U<-runif(nv[l])
    Y[(((l-1)*n/L)+1):(l*n/L)]<-qnorm(((l-1)/L+U/L),mean=mu,sd=sigma)
    
  }
  return(Y)
}
#stratification<-MC_estim.strat(n, 0, 2, 10)
#print(stratification)
pre_mc1.strat<-MC_estim.strat(n, 0, 2, 10)
pre_mc1.strat<-exp(pre_mc1.strat)+exp(X12)
pre_mc1.strat<-pre_mc1.strat>5
mc1.strat<-evol_estim(pre_mc1.strat)
mc1.strat_IC<-evol_IC(pre_mc1.strat, 0.95)
plot(1:n, mc1.strat,col='mediumseagreen',type='l', main = "Monte Carlo Strates", xlab="n iter", ylab ="estimation")
lines(2:n, mc1.strat_IC$min,col='blue')
lines(2:n, mc1.strat_IC$max,col='red')
legend("topleft", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)


#Résumé
mc1[n]
mc1.ant[n]
mc1.controle[n-p]
mc1.controle1[n-p]
mc1.strat[n]

ic_mc<-c(mc1_IC$min[n-1],mc1_IC$max[n-1])
ic_mc_ant<-c(mc1.ant_IC$min[n-1],mc1.ant_IC$max[n-1])
ic_cont<-c(mc1.controle_IC$min[n-p-1],mc1.controle_IC$max[n-p-1])
ic_cont1<-c(mc1.controle1_IC$min[n-p-1],mc1.controle1_IC$max[n-p-1])
ic_strat<-c(mc1.strat_IC$min[n-1],mc1.strat_IC$max[n-1])

ic_mc
ic_mc_ant
ic_cont
ic_cont1
ic_strat

var(pre_mc1)
var(pre_mc1.ant)
var(pre_mc1.controle)
var(pre_mc1.controle1)
var(pre_mc1.strat)

######
#II
######
#on va maintenant essayer d'estimer par différentes méthodes de Monte Carlo la valeur suivante :
#E(cos(X21*X22)*sin(X21)*exp(sin(X21+X22)))
par(mfrow = c(1,2))
X21<-vf2[1,]
X22<-vf2[2,]
#estimateur de Monte Carlo classique
pre_mc2<-cos(X21*X22)*sin(X21)*exp(sin(X21+X22))
mc2<-evol_estim(pre_mc2)
mc2_IC<-evol_IC(pre_mc2, 0.95)
plot(1:n, mc2,col='mediumseagreen',type='l',ylim = c(0, 0.2), main = "Monte Carlo Classique", xlab="n iter", ylab ="estimation")
lines(2:n, mc2_IC$min,col='blue')
lines(2:n, mc2_IC$max,col='red')
legend("topright", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur de Monte Carlo par variable antithétique
#Soit (X,Y) qui suit la loi de densité f2
#on remarque en étudiant E(h(-X,-y)) pour toute fonction h mesurable bornée que E(h(-X,-Y))=E(h(X,y)).
#Ceci est dû à la parité de cos et l'imparité de sin
Xbis21<-(-vf2[1,])
Xbis22<-(-vf2[2,])
pre_mc2.ant<-cos(Xbis21*Xbis22)*sin(Xbis21)*exp(sin(Xbis21+Xbis22))
pre_mc2.ant<-(pre_mc2+pre_mc2.ant)/2
mc2.ant<-evol_estim(pre_mc2.ant)
mc2.ant_IC<-evol_IC(pre_mc2.ant,0.95)
plot(1:n, mc2.ant,col='mediumseagreen',type='l',ylim = c(0, 0.2), main = "Monté Carlo Va Antithétiques", xlab="n iter", ylab ="estimation")
lines(2:n, mc2.ant_IC$min,col='blue')
lines(2:n, mc2.ant_IC$max,col='red')
legend("topright", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)


#
#
# PARTIE 3
#
#
#
rejet.recycle.f1 <- function(n) { #méthode de rejet recyclage pour simuler f1
  x1<- numeric(n)
  x2<- numeric(n)
  z1<- numeric(n)
  z2<- numeric(n)
  for (i in 1:n) {
    x1[i] <- 2*rnorm(1)
    x2[i] <- 2*rnorm(1)
    u <- runif(1)
    if (u > exp(-(x2[i]^2)/2+(x2[i]^2)/8)*(abs(x2[i]) <= 1)/2) { #on a multiplié M par 2 pour avoir f=M*g sur un espace de mesure nulle
      z1[i]<-x1[i]
      z2[i]<-x2[i]
      x1[i]<-0
      x2[i]<-0
    }
  }
  return(data.frame(accept1 = x1,accept2 = x2, rejet1 = z1, rejet2 = z2))
}
#Afin de n'avoir que des vecteurs pleins, on enlèvera les 0 des vecteurs créés.
#En effet, on a commencé avec des vecteurs nuls de taille n dans lesquels nous avons placé,
#à place i, soit la valeur acceptée soit celle rejettée, soit rien.
#Enlever les 0 n'a donc pas d'influence car on travaille avec des variables continues (i.e. 0 est tiré avec probabilité 0).
wf2<-rejet.recycle.f1(n)
acc1<-wf2$accept1
acc2<-wf2$accept2
rej1<-wf2$rejet1
rej2<-wf2$rejet2
acc1<-acc1[acc1!=0]
acc2<-acc2[acc2!=0]
rej1<-rej1[rej1!=0]
rej2<-rej2[rej2!=0]

#ici on teste si notre algorithme de rejet recyclé fonctionne bien en comparant les sorties 
#aux densités théoriques
hist(acc1,freq = F,breaks="FD",col="skyblue", main = "Histogramme de f1x")
t<- seq(-max(acc1), max(acc1), 0.01)
lines(t,dnorm(t,mean=0,sd=2),col=53,lwd=2)
legend("topright", "f1x(x)", box.lty = 0, bg = "gray90", col = "skyblue", lty = 1, lwd = 3, inset = 0.05)
hist(acc2, freq = FALSE, breaks = "FD", main = "Histogramme de f1y", ylab = "Fréquences")
t <- seq(-1, 1, 0.01)
lines(t, f1y(t), col = "mediumseagreen", lwd = 3)
legend("topright", "f1y(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)

densite_val_rejetes2<-function(y){
  gy<-exp(-(1/8)*y^2)/(2*sqrt(2*pi))
  fy<-exp(-(0.5)*y^2)*(abs(y)<= 1)/(sqrt(2*pi)*(pnorm(1)-pnorm(-1)))
  M<-4/(pnorm(1)-pnorm(-1))
  return((M*gy-fy)/(M-1))
}

densite_val_rejetes21<-function(y){
  gy<-(exp(-(1/8)*y^2)-0.5*exp(-(0.5)*y^2)*(abs(y)<= 1))/(2*sqrt(2*pi))
  M<-4/(pnorm(1)-pnorm(-1))
  return((M*gy)/(M-1))
}

hist(rej2,freq = F,breaks="FD",col="skyblue", main = "Histogramme de rejeté")
t<- seq(-max(rej2), max(rej2), 0.01)
lines(t,densite_val_rejetes21(t),col=53,lwd=2)
legend("topright", "f1x(x)", box.lty = 0, bg = "gray90", col = "skyblue", lty = 1, lwd = 3, inset = 0.05)

#construction de delta1 et delta2 à t fixé, observé comme la longueur de acc1 et acc2
par(mfrow = c(2,2))

evol_estim_delta1<-function(y){
  h<-(cumsum(y)/(2:(length(y)+1)))
  x<-rnorm(1, mean =0, sd=2)
  x_<-rejet.basique.f1y(1)
  fixe<-exp(x)+exp(x_)>5
  fixe2<-fixe*((f1x(x)*f1y(x_))/g(x,x_))/(2:(length(y)+1))
  return (fixe2+h)
}
evol_estim_delta2<-function(y){
  k<-length(y)
  h<-(cumsum(y)/(2:(k+1)))
  x<-rnorm(1, mean =0, sd=2)
  x_<-rejet.basique.f1y(1)
  fixe<-exp(x)+exp(x_)>5
  fixe2<-fixe*((f1x(x)*f1y(x_))/g(x,x_))/(2:(k+1))
  return (fixe2+h)
}

pre_delta1<-exp(acc1)+exp(acc2)>5
delta1<-evol_estim(pre_delta1)
delta1_IC<-evol_IC(pre_delta1, 0.95)
plot(1:length(delta1), delta1,col='mediumseagreen',type='l',ylim = c(0.2, 0.4) ,main = "Monte Carlo delta1", xlab="n iter", ylab ="estimation")
lines(2:length(delta1), delta1_IC$min,col='blue')
lines(2:length(delta1), delta1_IC$max,col='red')
legend("topright", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)
gbis<-function(x,y){#densité N(0,4)*N(0,4)
  return (dnorm(x, mean=0, sd=2)*dnorm(y, mean=0, sd=2))
}

M<-4/(pnorm(1)-pnorm(-1)) #on a multiplié M par 2 par rapport au rejet classique
pre_delta2.1<-exp(rej1)+exp(rej2)>5
pre_delta2.2<-pre_delta2.1*((M-1)*f1x(rej1)*f1y(rej2)/(M*gbis(rej1,rej2)-f1x(rej1)*f1y(rej2)))
delta2<-evol_estim(pre_delta2.2)
delta2_IC<-evol_IC(pre_delta2.2, 0.95)
plot(1:length(delta2), delta2,col='mediumseagreen',type='l',ylim = c(0.2, 0.4) ,main = "Monte Carlo delta2", xlab="n iter", ylab ="estimation")
lines(2:length(delta2), delta2_IC$min,col='blue')
lines(2:length(delta2), delta2_IC$max,col='red')
legend("topright", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)







####
####
#On veut poser 'a' egal à var(d2)/(var(d1)+var(d2)).
#on estime a sur les n/20 premières valeurs pour éviter d'avoir des corrélations et du biais.
p1<-n/20
delta1_av<-delta1[1:p1]
length(delta1_av)<-p1
delta2_av<-delta2[1:p1]
length(delta2_av)<-p1
esp_d1<-sum(delta1_av)/(length(delta1_av))
esp_d2<-sum(delta2_av)/(length(delta2_av))
var_d1<-sum(((delta1_av-esp_d1)^2))/(length(delta1_av)-1)
var_d2<-sum(((delta2_av-esp_d2)^2))/(length(delta2_av)-1)
a<-var_d2/(var_d1+var_d2)

var(pre_delta1)
var(pre_delta2.2)


pre_delta1_tronc<-pre_delta1[p1+1:length(delta1)]
length(pre_delta1_tronc)<-length(delta1)-p1
v1<-var(pre_delta1_tronc)
length(pre_delta1_tronc)

pre_delta2_tronc<-pre_delta2.2[p1+1:length(delta2)]
length(pre_delta2_tronc)<-length(delta2)-p1
v2<-var(pre_delta2_tronc)


estim_d1_d2_a<-a*sum(pre_delta1_tronc)/(length(pre_delta1_tronc))+(1-a)*sum(pre_delta2_tronc)/length(pre_delta2_tronc)
estim_d1_d2_a
esp_a_d1_d2<-a*esp_d1+(1-a)*esp_d2
var_a_d1_d2<-a*a*v1+(1-a)*(1-a)*v2

#les différentes variances
v1
v2
var_a_d1_d2 #on observe que la variance est réduite en recyclant



####
####
par(mfrow = c(1,2))
rejet.recycle.f3 <- function(n) { #méthode de rejet recyclage pour simuler f3
  x1<- numeric(n)
  x2<- numeric(n)
  z1<- numeric(n)
  z2<- numeric(n)
  for (i in 1:n) {
    x1[i] <- 2*rnorm(1)
    x2[i] <- rnorm(1)
    u <- runif(1)
    if (u > (1/3)*((abs(x2[i]) <= 1)+0.5)) { #on prend M1=2M pour avoir f=M1*g sur un espace de mesure nulle
      z1[i]<-x1[i]
      z2[i]<-x2[i]
      x1[i]<-0
      x2[i]<-0
    }
  }
  return(data.frame(accept3 = x1,accept4 = x2, rejet3 = z1, rejet4 = z2))
}

#on trace la densité de f3y
f3<-function(x,y){
  return(exp(-(0.5)*((x^2)/4+(y^2)))*((abs(y)<= 1)+0.5)/(4*pi*(pnorm(1)-pnorm(-1)+0.5)))
}

f3x<-function(x){
  return(dnorm(x, mean=0, sd=2))
}
f3y<-function(x){
  return(dnorm(x)*((abs(x)<= 1)+0.5)/(pnorm(1)-pnorm(-1)+0.5))
}

wf3<-rejet.recycle.f3(n)
acc3<-wf3$accept3
acc4<-wf3$accept4
rej3<-wf3$rejet3
rej4<-wf3$rejet4
acc3<-acc3[acc3!=0]
acc4<-acc4[acc4!=0]
rej3<-rej3[rej3!=0]
rej4<-rej4[rej4!=0]



hist(acc4, freq = FALSE, breaks = 50, main = "Histogramme de f3y", ylab = "Fréquences")
s <- seq(-max(acc4), max(acc4), 0.01)
lines(s, f3y(s), col = "mediumseagreen", lwd = 3)
legend("topright", "f3y(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)

densite_val_rejetes3<-function(y){#densité y de (M*g-f)/(M-1)
  g3<-exp(-(0.5)*(y^2))/sqrt(2*pi)
  fy3<-(1/3)*(exp(-(0.5)*(y^2))*((abs(y)<= 1)+0.5))/(sqrt(2*pi))
  M3<-3/(0.5+pnorm(1)-pnorm(-1))
  return((M3*(g3-fy3))/(M3-1)) #retourne après calcul (Mg-f)/(M-1) en y
}
hist(rej4,freq = F,breaks="FD",col="skyblue", main = "Histogramme de f3y rejeté ")
r<- seq(-max(rej4), max(rej4), 0.01)
lines(r,densite_val_rejetes3(r),col=53,lwd=2)
legend("topright", "f1x(x)", box.lty = 0, bg = "gray90", col = "skyblue", lty = 1, lwd = 3, inset = 0.05)

par(mfrow = c(1,2))
#calcul des estimateurs delta1 et delta2 pour f3
pre_delta1.3<-exp(acc3)+exp(acc4)>5
delta1.3<-evol_estim(pre_delta1.3)
delta1.3_IC<-evol_IC(pre_delta1.3, 0.95)
plot(1:length(delta1.3), delta1.3,col='mediumseagreen',type='l' ,ylim = c(0.2, 0.4) ,main = "Monte Carlo delta1.3", xlab="n iter", ylab ="estimation")
lines(2:length(delta1.3), delta1.3_IC$min,col='blue')
lines(2:length(delta1.3), delta1.3_IC$max,col='red')
legend("topright", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)
aide_d2<- function(x,y){
  h<-dnorm(x, mean=0, sd=2)*densite_val_rejetes3(y)
  return(f3(x,y)/h)
}
M2<-3/(0.5+pnorm(1)-pnorm(-1)) #on met 3 au lieu de 3/2 car on a pris 2M
pre_delta2.3<-exp(rej3)+exp(rej4)>5
pre_delta2.4<-pre_delta2.3*aide_d2(rej3,rej4)
delta2.3<-evol_estim(pre_delta2.4)
delta2.3_IC<-evol_IC(pre_delta2.4, 0.95)
plot(1:length(delta2.3), delta2.3,col='mediumseagreen',type='l',ylim = c(0.2, 0.4) ,main = "Monte Carlo delta2.3", xlab="n iter", ylab ="estimation")
lines(2:length(delta2.3), delta2.3_IC$min,col='blue')
lines(2:length(delta2.3), delta2.3_IC$max,col='red')
legend("topright", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)






delta1.3_av<-delta1.3[1:p]
length(delta1.3_av)<-p
delta2.3_av<-delta2.3[1:p]
length(delta2.3_av)<-p
esp_d1.3<-sum(delta1.3_av)/(length(delta1.3_av))
esp_d2.3<-sum(delta2.3_av)/(length(delta2.3_av))
var_d1.3<-sum(((delta1.3_av-esp_d1.3)^2))/(length(delta1.3_av)-1)
var_d2.3<-sum(((delta2.3_av-esp_d2.3)^2))/(length(delta2.3_av)-1)
a.3<-var_d2.3/(var_d1.3+var_d2.3)

var(pre_delta1.3)
var(pre_delta2.4)


pre_delta1.3_tronc<-pre_delta1.3[p+1:length(delta1.3)]
length(pre_delta1.3_tronc)<-length(delta1.3)-p
pre_delta2.3_tronc<-pre_delta2.4[p+1:length(delta2.3)]
length(pre_delta2.3_tronc)<-length(delta2.3)-p

estim_d1_d2_a.3<-a.3*sum(pre_delta1.3_tronc)/(length(pre_delta1.3_tronc))+(1-a.3)*sum(pre_delta2.3_tronc)/length(pre_delta2.3_tronc)
var_a_d1_d2.3<-a.3*a.3*var(pre_delta1.3_tronc)+(1-a.3)*(1-a.3)*var(pre_delta2.3_tronc)



#les différentes variances
var(pre_delta1.3_tronc)
var(pre_delta2.3_tronc)
var_a_d1_d2.3 #on observe une variance réduite en recyclant 



# PARTIE 4
#
#
#simulation chaine de markov f1
par(mfrow = c(1,2))
simu_markov_MH1<-function(n){
  x1<-numeric(n)
  x2<-numeric(n)
  temp<-rejet.basique.f1(1)
  x1[1]<-temp[1,1]
  x2[1]<-temp[2,1]
  for (i in 2:n){
    eps1<-2*rnorm(1)
    eps2<-rnorm(1)
    obj<-f1x(eps1)*f1y(eps2)*g(x1[i-1],x2[i-1])/(f1x(x1[i-1])*f1y(x2[i-1])*g(eps1,eps2))
    alpha<-min(1,obj)
    u<-rbinom(1,1,alpha)
    if (u==1){ #U est une bernoulli de paramètre alpha i.e. P(U=1)=alpha
      x1[i]<- eps1
      x2[i]<- eps2
    }
    else{
      x1[i]<- x1[i-1]
      x2[i]<- x2[i-1]
    }
  }
  v<-rbind(x1,x2)
  return(v)
}


vf_markov1<-simu_markov_MH1(n)
hist(vf_markov1[2,], freq = FALSE, breaks = 25, main = "Histogramme de f1y_markov", ylab = "Fréquences")
t <- seq(-1, 1, 0.01)
lines(t, f1y(t), col = "mediumseagreen", lwd = 3)
legend("topright", "f1y_M(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)
##############
hist(vf_markov1[1,], freq = FALSE, breaks = 30, main = "Histogramme de f1x_markov", ylab = "Fréquences")
s <- seq(-max(vf_markov1[1,]), max(vf_markov1[1,]), 0.01)
lines(s, f1x(s), col = "mediumseagreen", lwd = 3)
legend("topright", "f1x_M(x)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.05)


#
#vérification de l'approximation (1)
#
par(mfrow = c(1,3))

X11_M<-vf_markov1[1,]
X12_M<-vf_markov1[2,]
#estimateur de Monte Carlo classique pour Markov
pre_mc1_M<-exp(X11_M)+exp(X12_M)
pre_mc1_M<-pre_mc1_M>5
mc1_M<-evol_estim(pre_mc1_M)
mc1_M_IC<-evol_IC(pre_mc1_M, 0.95)
plot(1:n, mc1_M,col='mediumseagreen',type='l',ylim = c(0.2, 0.4), main = "MC Classique_Markov", xlab="n iter", ylab ="estimation")
lines(2:n, mc1_M_IC$min,col='blue')
lines(2:n, mc1_M_IC$max,col='red')
legend("topright", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur de Monte Carlo par variable antithétique pour Markov
pre_mc1_M.ant<-exp(-X11_M)+exp(-X12_M)
pre_mc1_M.ant<-pre_mc1_M.ant>5
pre_mc1_M.ant<-(pre_mc1_M.ant+pre_mc1_M)/2
mc1_M.ant<-evol_estim(pre_mc1_M.ant)
mc1_M.ant_IC<-evol_IC(pre_mc1_M.ant, 0.95)
plot(1:n, mc1_M.ant,col='mediumseagreen',type='l',ylim = c(0.2, 0.4),main = "MC Va Antithétiques_Markov", xlab="n iter", ylab ="estimation")
lines(2:n, mc1_M.ant_IC$min,col='blue')
lines(2:n, mc1_M.ant_IC$max,col='red')
legend("topright", legend=c("MC.ant", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur avec une variable de controle
#on va noter Z la fonction indicatrice(X11+X12>log(25)) 
#on note g la fonction indicatrice(exp(X11)+exp(X12)>5)
p<-n/50
Z_M<-X11_M+X12_M>log(2)
esp_Z_M<-sum(Z_M[1:p])/(p)
esp_g_M<-sum(pre_mc1_M[1:p])/(p)
cov_gZ_M<-sum((pre_mc1_M[1:p]-esp_g_M)*(Z_M[1:p]-esp_Z_M))
var_Z_M<-sum((Z_M[1:p]-esp_Z_M)*(Z_M[1:p]-esp_Z_M))
b_M<-cov_gZ_M/var_Z_M
Z_Mtronc<-Z_M[p+1:n]
length(Z_Mtronc)<-n-p #changement manuel de la taille de Ztronc
pre_mc1_M.tronc<-pre_mc1_M[p+1:n]
length(pre_mc1_M.tronc)<-n-p #changement manuel de la taille de pre_mc.tronc
pre_mc1_M.controle<-pre_mc1_M.tronc-b_M*(Z_Mtronc-esp_Z_M)
mc1_M.controle<-evol_estim(pre_mc1_M.controle)
mc1_M.controle_IC<-evol_IC(pre_mc1_M.controle, 0.95)
m_M<-p+1
plot(m_M:n, mc1_M.controle,col='mediumseagreen',type='l',ylim = c(0.2, 0.4) ,main = "MC avec controle_Markov", xlab="n iter", ylab ="estimation")
lines((m_M+1):n, mc1_M.controle_IC$min,col='blue')
lines((m_M+1):n, mc1_M.controle_IC$max,col='red')
legend("topright", legend=c("MC.controle", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#les différentes variances
var(pre_mc1_M)
var(pre_mc1_M.ant)
var(pre_mc1_M.controle)


#simulation chaine de markov f2
par(mfrow = c(1,2))
simu_markov_MH2<-function(n){
  y1<-numeric(n)
  y2<-numeric(n)
  temp<-rejet.basique.f2(1)
  y1[1]<-temp[1,1]
  y2[1]<-temp[2,1]
  for (i in 2:n){
    eps1<-2*rnorm(1)
    eps2<-rnorm(1)
    obj<-f2(eps1,eps2)*g(y1[i-1],y2[i-1])/(f2(y1[i-1],y2[i-1])*g(eps1,eps2)) #alpha est défini avec f2, mais comme f2=c*f2tild
    #et que l'on divise dans alpha f2(eps) par f2(x[i-1]), cela reviens à diviser f2tild(eps) par f2tild(x[i-1])
    alpha<-min(1,obj)
    u<-rbinom(1,1,alpha)
    if (u==1){ #U est une bernouilli de parametre alpha i.e. P(U=1)=alpha
      y1[i]<- eps1
      y2[i]<- eps2
    }
    else{
      y1[i]<- y1[i-1]
      y2[i]<- y2[i-1]
    }
  }
  w<-rbind(y1,y2)
  return(w)
}


vf_markov2<-simu_markov_MH2(n)
plot(vf2[1,], vf2[2,]) #pour comparer avec le tirage en rejet normal de f2
plot(vf_markov2[1,], vf_markov2[2,]) #on remarque que les deux se ressemblent bien 

par(mfrow = c(1,2))

#estimateur de Monte Carlo classique  avec Markov
X21_M<-vf_markov2[1,]
X22_M<-vf_markov2[2,]
pre_mc2_M<-cos(X21_M*X22_M)*sin(X21_M)*exp(sin(X21_M+X22_M))
mc2_M<-evol_estim(pre_mc2_M)
mc2_M_IC<-evol_IC(pre_mc2_M, 0.95)
plot(1:n, mc2_M,col='mediumseagreen',type='l',ylim = c(0, 0.2), main = "MC Classique Markov", xlab="n iter", ylab ="estimation")
lines(2:n, mc2_M_IC$min,col='blue')
lines(2:n, mc2_M_IC$max,col='red')
legend("topright", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)

#estimateur de Monte Carlo par variable antithétique
#Soit (X,Y) qui suit la loi de densité f2
#on remarque en étudiant E(h(-X,-y)) pour tout fonction h mesurable bornée que E(h(-X,-Y))=E(h(X,y)).
#Ceci est du à la parité de cos et l'imparité de sin
Xbis21_M<-(-vf_markov2[1,])
Xbis22_M<-(-vf_markov2[2,])
pre_mc2_M.ant<-cos(Xbis21_M*Xbis22_M)*sin(Xbis21_M)*exp(sin(Xbis21_M+Xbis22_M))
pre_mc2_M.ant<-(pre_mc2_M+pre_mc2_M.ant)/2
mc2_M.ant<-evol_estim(pre_mc2_M.ant)
mc2_M.ant_IC<-evol_IC(pre_mc2_M.ant,0.95)
plot(1:n, mc2_M.ant,col='mediumseagreen',type='l',ylim = c(0, 0.2), main = "MC Va Antithétiques Markov", xlab="n iter", ylab ="estimation")
lines(2:n, mc2_M.ant_IC$min,col='blue')
lines(2:n, mc2_M.ant_IC$max,col='red')
legend("topright", legend=c("MC", "IC_min", "IC_max"), box.lty = 0, 
       col = c("mediumseagreen", "blue", "red"),lty = 1, lwd = 3)


#les différentes variances
var(pre_mc2_M)
var(pre_mc2_M.ant)