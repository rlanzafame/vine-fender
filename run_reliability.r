library(base)
library(copula)
library(datasets)
library(DISTRIB)
library(distributions3)
library(fitdistrplus)
library(graphics)
library(grDevices)
library(invgamma)
library(lattice)
library(lestat)
library(lme4)
library(lmomco)
library(Lmoments)
library(lsei)
library(MASS)
library(Matrix)
library(methods)
library(mvtnorm)
library(npsurv)
library(pracma)
library(rvinecopulib)
library(stats)
library(stats4)
library(survival)
library(utils)
library(VineCopula)
library(xlsx)

# Coba form manual
mass <- read.csv("mass1.csv", header=FALSE)
mass <- as.numeric(unlist(mass))
velocity <- read.csv("vel1.csv", header=FALSE)
velocity <- as.numeric(unlist(velocity))
angle <- read.csv("angle1.csv", header=FALSE)
angle <- as.numeric(unlist(angle))
matrix <- cbind(mass,angle,velocity)

RVM <- RVineStructureSelect(matrix, familyset=1:5, selectioncrit = "AIC",progress = TRUE)

#Change the value here
Ecv <- 1950 #Fender rated capacity in kNm
deflection <- 1600 #Fender rated deflection in mm

#Define weibull parameters for the velocity
a=1.7
b=4.38


#first trial
y1 <- 1
y2 <- 1.6
y3 <- 0
y4 <- 0
y5 <- 0

y <- c(y1,y2,y3,y4,y5)
e1 <- 1
g<-10

step <- 0.1

while (g>0.09)
{
  #Find CDF
  
  U1 <- pnorm(y[1],0,1)
  U21 <- pnorm(y[2],0,1)
  U2  <- BiCopHinv1(U1,U21,family=34,-1.69)
  U312 <- pnorm(y[3],0,1)
  U12 <- BiCopHfunc2(U1,U2,family=34,-1.69)
  U32 <- BiCopHinv1(U12,U312,family=5,1.35)
  U3 <- BiCopHinv1(U2,U32,family=33,-1.12)
  U4 <- pnorm(y[4],0,1)
  U5 <- pnorm(y[5],0,1)
  
  #Find X
  
  x1 <- (b*(-log(1-U1))^(1/a))
  x2 <- invcdf(uniformdistribution(8000,260000),U2)
  x3 <- invcdf(gammadistribution(1.19,1/0.27),U3)
  x4 <- qnorm(U4,1,0.033)
  x5 <- qnorm(U5,15,5)
  
  #iterasi mulai dari sini
  v     <- x1
  M     <- x2
  angles<- x3
  MF    <- x4
  temp  <- x5
  
  # Construct jacobian matrix
  Jyx11 <- dweibull(x1,a,b)/dnorm(y1,0,1)
  Jyx22 <- BiCopPDF(U1,U2,family=34,par=-1.69)*dunif(x2,8000,260000)/dnorm(y2,0,1)
  Jyx33 <- BiCopPDF(U12,U32,family=5,1.35)*BiCopPDF(U2,U3,family=33,-1.12)*dgamma(x3,1.19,0.27)/dnorm(y3,0,1)
  Jyx44 <- dnorm(x4,1,0.033)/dnorm(y4,0,1)
  Jyx55 <- dnorm(x5,15,5)/dnorm(y5,0,1)
  
  # cari Jyx12 pake numerical differentiation (turunan terhadap x1)
  h       <-  0.01
  U1h     <-  pweibull(x1+h,a,b)
  U21h    <-  BiCopHfunc1(U1h,U2,family=34,-1.69)
  dF21dx1 <-  (U21h-U21)/h
  Jyx12   <-  dF21dx1/dnorm(y2,0,1)
  
  # cari Jyx13 pake numerical differentiation (turunan terhadap x1)
  U12h     <- BiCopHfunc2(U1h,U2,family=34,-1.69)
  U312h    <- BiCopHfunc1(U12h,U32,family=5,1.35)
  df312dx1 <- (U312h-U312)/h
  Jyx13    <- df312dx1/dnorm(y3,0,1)
  
  # cari Jyx23 numerical differentiation (turunan terhadap x2)
  h        <- 10
  U2h      <- punif(x2+h,8000,260000)
  U12h     <- BiCopHfunc2(U1,U2h,family=34,-1.69)
  U32h     <- BiCopHfunc1(U2h,U3,family=33,-1.12)
  U312h    <- BiCopHfunc1(U12h,U32h,family=5,1.35)
  df312dx2 <- (U312h-U312)/h
  Jyx23  <- df312dx2/dnorm(y3,0,1)
  
  Jyx <- rbind(c(Jyx11,0,0,0,0),c(Jyx12,Jyx22,0,0,0),c(Jyx13,Jyx23,Jyx33,0,0),c(0,0,0,Jyx44,0),c(0,0,0,0,Jyx55))
  Jxy <- inv(Jyx)
  
  Lpp   <- ((M)^0.3858)*3.279
  D     <- ((M)^0.2321)*0.9519
  B     <- ((M)^0.33339)*0.8528
  Cb    <- M/(Lpp*B*D*1.02)
  R     <- sqrt((Lpp/4)^2+(B/2)^2)
  K     <- (0.19*Cb+0.11)*Lpp
  gamma <- 90-angles-asin(B/(2*R))
  Ce    <- (K^2+R^2*(cospi(gamma/180))^2)/(K^2+R^2)
  Cm    <- 1+pi*D/(2*Cb*B)
  VF    <- -0.056*log(0.72*deflection/(0.74*v*10))+1.2011
  TF    <- -3*10^-6*(temp^3)+0.0002*(temp^2)-0.0058*temp+1.0533
  
  g     <- Ecv*VF*MF*TF-0.5*M*(v/100)^2*Cm*Ce
  
  #gradg1
  hv     <- 0.1
  vh     <- v+hv
  VFh    <- -0.056*log(0.72*deflection/(0.74*vh*10))+1.2011
  gv1    <- Ecv*VFh*MF*TF-0.5*M*(vh/100)^2*Cm*Ce
  gradg1 <- (gv1-g)/hv
  
  #gradg2
  hM <- 10
  Mh <- M+hM
  Lpph  <- ((Mh)^0.3858)*3.279
  Dh     <- ((Mh)^0.2321)*0.9519
  Bh     <- ((Mh)^0.33339)*0.8528
  Cbh    <- Mh/(Lpph*Bh*Dh*1.02)
  Rh     <- sqrt((Lpph/4)^2+(Bh/2)^2)
  Kh     <- (0.19*Cbh+0.11)*Lpph
  gammah <- 90-angles-asin(Bh/(2*Rh))
  Ceh    <- (Kh^2+Rh^2*(cospi(gamma/180))^2)/(Kh^2+Rh^2)
  Cmh    <- 1+pi*Dh/(2*Cbh*Bh)
  gM2     <- Ecv*VF*MF*TF-0.5*Mh*(v/100)^2*Cmh*Ceh
  gradg2 <- (gM2-g)/hM
  
  #gradg3
  hang      <- 0.01
  anglesh   <- angles+hang
  gammahang <- 90-anglesh-asin(B/(2*R))
  Cehang    <- (K^2+R^2*(cospi(gammahang/180))^2)/(K^2+R^2)
  ghang     <- Ecv*VF*MF*TF-0.5*M*(v/100)^2*Cm*Cehang
  gradg3    <- (ghang-g)/hang
  
  #gradg4 terhadap MF
  hMF        <- 0.1
  MFgrad4    <- MF+hMF
  ggrad4     <- Ecv*VF*TF*MFgrad4-0.5*M*(v/100)^2*Cm*Ce
  gradg4     <- (ggrad4-g)/hMF
  
  #gradg5 terhadap temperature
  htemp      <- 0.1
  tempgrad5  <- temp+htemp
  TFgrad5    <- -3*10^-6*(tempgrad5^3)+0.0002*(tempgrad5^2)-0.0058*tempgrad5+1.0533
  ggrad5     <- Ecv*VF*TFgrad5 *MF-0.5*M*(v/100)^2*Cm*Ce
  gradg5     <- (ggrad5-g)/htemp
  
  gradgx <- c(gradg1,gradg2,gradg3,gradg4,gradg5)
  
  Grad_G <- gradgx%*%Jxy
  alpha  <- -Grad_G/(sqrt(sum(Grad_G^2)))
  
  d1     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[1])-y[1]
  d2     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[2])-y[2]
  d3     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[3])-y[3]
  d4     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[4])-y[4]
  d5     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[5])-y[5]
  d      <- c(d1,d2,d3,d4,d5)
  
  
  ynew   <-  y+step*d
  
  # Find new point xnew
  U1new   <- pnorm(ynew[1],0,1)
  U21new  <- pnorm(ynew[2],0,1)
  U2new   <- BiCopHinv1(U1new,U21new,family=34,-1.69)
  U312new <- pnorm(ynew[3],0,1)
  U12new  <- BiCopHfunc2(U1new,U2new,family=34,-1.69)
  U32new  <- BiCopHinv1(U12new,U312new,family=5,1.35)
  U3new   <- BiCopHinv1(U2new,U32new,family=33,-1.12)
  U4new   <- pnorm(ynew[4],0,1)
  U5new   <- pnorm(ynew[5],0,1)
  
  #Find X
  
  x1new <- (b*(-log(1-U1new))^(1/a))
  x2new <- invcdf(uniformdistribution(8000,260000),U2new)
  x3new <- invcdf(gammadistribution(1.19,1/0.27),U3new)
  x4new <- qnorm(U4new,1,0.033)
  x5new <- qnorm(U5new,15,5)
  
  vnew      <- x1new
  Mnew      <- x2new
  anglesnew <- x3new
  MFnew     <- x4new
  tempnew   <- x5new
  
  Lppnew   <- ((Mnew)^0.3858)*3.279
  Dnew     <- ((Mnew)^0.2321)*0.9519
  Bnew     <- ((Mnew)^0.33339)*0.8528
  Cbnew    <- Mnew/(Lppnew*Bnew*Dnew*1.02)
  Rnew     <- sqrt((Lppnew/4)^2+(Bnew/2)^2)
  Knew     <- (0.19*Cbnew+0.11)*Lppnew
  gammanew <- 90-anglesnew-asin(Bnew/(2*Rnew))
  Cenew    <- (Knew^2+Rnew^2*(cospi(gammanew/180))^2)/(Knew^2+Rnew^2)
  Cmnew    <- 1+pi*Dnew/(2*Cbnew*Bnew)
  VFnew    <- -0.056*log(0.72*deflection/(0.74*vnew*10))+1.2011
  TFnew    <- -3*10^-6*(tempnew^3)+0.0002*(tempnew^2)-0.0058*tempnew+1.0533
  gnew     <- Ecv*VFnew*MFnew*TFnew-0.5*Mnew*(vnew/100)^2*Cmnew*Cenew
  
  ci    <- sqrt(sum(y^2))/sqrt(sum(Grad_G ^2))+10
  
  mold <- 0.5*(sqrt(sum(y^2)))^2+ci*abs(g)
  mnew <- 0.5*(sqrt(sum(ynew^2)))^2+ci*abs(gnew)
  
  while(mold<mnew){stop("stop")}
    
  print(g)
  
  y  <- ynew
  }


beta=(sqrt(sum(ynew^2)))

#hitung alpha asli
alpha1<-qnorm(U1,0,1)/beta
alpha2<-qnorm(U2,0,1)/beta
alpha3<-(qnorm(U3,0,1))/beta

alphacorrected<-c(alpha[1],alpha2,alpha3,alpha[4],alpha[5])
designpointsx<-c(x1new,x2new,x3new,x4new,x5new)
designpointsu<-c(U1new,U2new,U3new,U4new,U5new)

sigmahat<- Jxy%*%t(Jxy)
diagonal<-diag(sigmahat)
dakar<-sqrt(diagonal)
diagonal<- rbind(c(dakar[1],0,0,0,0),c(0,dakar[2],0,0,0),c(0,0,dakar[3],0,0),c(0,0,0,dakar[4],0),c(0,0,0,0,dakar[5]))

atas<- alpha%*%Jyx%*%diagonal

gamma1<-atas/(sqrt(sum(atas^2)))


#Printresult
print(beta)
print(alphacorrected)
print(designpointsu)
print(designpointsx)
