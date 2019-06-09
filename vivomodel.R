library(RxODE)
library(nlmixr)
library(xpose)
library(xpose.nlmixr)
library(ggplot2)

setwd("~/StaffFellow/Projects/MIDD/Simulations/2")

rm(list=ls(all=TRUE)) 
graphics.off()

dat <- read.table("EBOV.txt", head=TRUE)
str(dat)
ggplot(dat, aes(TIME, DV)) + geom_line(aes(group=id), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 Viremia (cp/mL)") + labs(title="Viremia", subtitle="Viremia vs. time by individual")

vivo <- function() {
  ini({

    ld<-log(0.224)
    lp<-log(4.15E4)
    lV0<-log(7.9E-5)
    lb<-log(7.9E-11)
    lq<-log(0.0074)
    lTt<-log(1.73)
    lfi<-log(2.67)
    lE0<-log(36900)
    lzeta<-log(0.455)
    lr<-log(0.338)
    ls<-log(3050)
    lTe<-log(6.5E-4)
    lka<-log(2.08E-5)
    lql<-log(0.0097)
    lqn<-log(0.0046)

    eta.d ~ 0.028561
    eta.p ~ 3.0976
    eta.V0 ~ 2.25
    eta.Tt ~ 1.1025
    eta.E0 ~ 0.600625
    eta.ql ~ 0.5776
    eta.qn ~ 1.0609
    
    add.err <- 1
  })
  model({
    
    ## Parameters------------------------
    d <- exp(ld + eta.d)
    p <- exp(lp + eta.p)
    V0 <- exp(lV0 + eta.V0)    
    b<-exp(lb)
    q<-exp(lq)
    Tt <- exp(lTt + eta.Tt)
    fi<-exp(lfi)
    E0 <- exp(lE0 + eta.E0)
    zeta<-exp(lzeta)
    r<-exp(lr)
    s<-exp(ls)
    Te<-exp(lTe)
    ka<-exp(lka)
    ql <- exp(lql + eta.ql)
    qn <- exp(lqn + eta.qn)

    ## Fixed parameters------------------------ 
    c = 20;      
    T0 = 1E8;    
    P0 = 0.001; 
    k = 4; 
    de = 0.001; 
    dfe = 0.4; 
    dl = 0.4; 
    dn = 0.4;
    sca=log(10.0);

    ## Initial Values------------------------ 
    T(0) <- T0
    V(0) <- V0
    E1(0) <- E0
    E2(0) <- P0*E0
    
    ##-----------------------------------------       
    d/dt(T) =-b*T*V - fi*T*F/(F+Tt);
    d/dt(I1) = b*T*V - k*I1;
    d/dt(I2) = k*I1 - d*I2 - ka*I2*E2;
    d/dt(R) = fi*T*F/(F+Tt);
    d/dt(V) = p*I2 - c*V;
    d/dt(F) = q*I2 - dfe*F;
    d/dt(L) = ql*I2 - dl*L;
    d/dt(N) = qn*I2 - dn*N;
    d/dt(E1) = s - zeta*F*E1/(F+Te) - de*E1;
    d/dt(E2) = r*E2*(1-E2/E0)-de*E2;
    
    Vir=log(V)/sca
    Vir ~ add(add.err)
  })
}

# fit <- nlmixr(vivo, dat, est="saem")
fit <- nlmixr(vivo, dat, est="saem", saemControl(n.burn=90, n.em=50))

plot(fit)
traceplot(fit)

## get xpdb
xpdb <- xpose_data_nlmixr(fit)

# Plot VPC
vpc.ui(fit,n=500)
