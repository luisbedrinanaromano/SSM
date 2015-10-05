#===================Librerias requeridas=======================
#==============================================================
library(rjags)
library(R2jags)
library(CircStats)
library(runjags)
library(coda)
library(lattice)
library(adehabitatLT)
library(sp)


#Leer el archivo con los datos
B84481.dat<-read.csv2("84481-Locs2.csv")
#Elimino observaciones con fechas identicas
B84481.dat<-B84481.dat[-c(47,159),]
head(B84481.dat)



#Ingresar el tiempo en la clase POSIXct
da <- as.POSIXct(B84481.dat$Date)
da

#Ahora transformo los datos a un objeto clase ltraj, lo cual
#inmediatamente calcula estadisticas descritivas necesarias
#como el angulo en radianes y la diferencia de tiempo "dt"
#entre relocaciones.
B84481.data <- as.ltraj(xy = B84481.dat[,c("X","Y")], date = da, id = B84481.dat$DeployID)

#Ahora transformo los datos a un data.frame
#normal, con el fin de filtrarlos.
B84481.data<-ld(B84481.data)
B84481.data<-na.omit(B84481.data)


#Ahora genero una nueva variable "mr", dividiendo la 
#distancia entre relocaciones y el tiempo entre ellas. 
#Expresado en m/s
B84481.data$mr<-B84481.data$dist/B84481.data$dt


#========================B84481=============================
#============================================================
require(runjags)
ones <- rep(1,500)
B84481<- dump.format(list(N=length(B84481.data$abs.angle), ta=B84481.data$abs.angle,
                          mr=B84481.data$mr, ones=ones))

B84481.data.l<-list(N=length(B84481.data$abs.angle), 
                    ta=B84481.data$abs.angle,mr=B84481.data$mr, ones=ones)

#Modelo doble
doble.jg<-function(){
  ## likelihood 
  for(t in 1:N){
    idx[t] ~ dcat(nu[])
    mr[t] ~ dgamma(b[idx[t]], a[idx[t]])
    ones[t] ~ dbern(wc[t])
    wc[t] <- (1/(2*Pi) * (1-rho[idx[t]]*rho[idx[t]]) / (1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(ta[t]-mu[idx[t]])))/500 
  }
  mn[1] ~ dunif(0,14)
  eps ~ dunif(0,14)
  mn[2] <- mn[1] + eps
  b[1] ~ dgamma(0.001, 0.001)
  b[2] ~ dgamma(0.001, 0.001)
  a[1] <- b[1]/mn[1]
  a[2] <- b[2]/mn[2]
  mu[1] ~ dunif(-32.98672,32.98672)
  mu[2] ~ dunif(-32.98672,32.98672)
  eps3 ~ dunif(0,1)
  eps4 ~ dunif(0,1)
  rho[1] <- min(eps3,eps4)
  rho[2] <- max(eps3,eps4)
  p ~ dunif(0,1)
  nu[1] <- p
  nu[2] <- 1 - p
  Pi <- 3.14159265359
}

#=================Especificaciones del MCMC===============
#=========================================================

nc<-3 #numero de cadenas
ni<-60000 #numero de iteraciones
nb<-10000 #number of burnin, numero de iteraciones iniciales descartadas para la
#estimacion de los parametros
nt<-25 #thinning rate, el valor de nt detemina cada cuantos valores el programa
# guardara valores para estimar los parametros. 1 implica guardar todos los valores
#en este caso guardar? un valor cada 25 datos

#Parametros para el modelo doble
d.params<-c("a[1]","a[2]","b[1]","b[2]","mu[1]","mu[2]","rho[1]","rho[2]","idx")

B84481.out<-jags(data = B84481.data.l, model.file = doble.jg, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parameters.to.save = d.params)
B84481.out