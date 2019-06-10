library(RxODE)
library(nlmixr)
library(xpose)
library(xpose.nlmixr)
library(ggplot2)

setwd("~/StaffFellow/Projects/MIDD/Simulations/2")

rm(list=ls(all=TRUE)) 
graphics.off()

dat <- read.table("EBOV.txt", head=TRUE)
# dat$DV = (10^dat$DV)
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
    #lepsi<-log(0)
    
    eta.d ~ 0.028561
    eta.p ~ 3.0976
    eta.V0 ~ 2.25
    # eta.b ~ 0.01 ##!
    # eta.q ~ 0.01 ##!
    eta.Tt ~ 1.1025
    # eta.fi ~ 0.01 ##!
    eta.E0 ~ 0.600625
    # eta.zeta ~ 0.01 ##!
    # eta.r ~ 0.01 ##!
    # eta.s ~ 0.01 ##!
    # eta.Te ~ 0.01 ##!
    # eta.ka ~ 0.01 ##!
    eta.ql ~ 0.5776
    eta.qn ~ 1.0609
    
    add.err <- 1
  })
  model({
    
    d <- exp(ld + eta.d)
    p <- exp(lp + eta.p)
    viremia0 <- exp(lV0 + eta.V0)    
    b<-exp(lb)# + eta.b) ##!
    q<-exp(lq)# + eta.q) ##!
    Tt <- exp(lTt + eta.Tt)
    fi<-exp(lfi)# + eta.fi) ##!
    CD8T <- exp(lE0 + eta.E0)
    zeta<-exp(lzeta)# + eta.zeta) ##!
    r<-exp(lr)# + eta.r) ##!
    s<-exp(ls)# + eta.s) ##!
    Te<-exp(lTe)# + eta.Te) ##!
    ka<-exp(lka)# + eta.ka) ##!
    ql <- exp(lql + eta.ql)
    qn <- exp(lqn + eta.qn)
    #epsi<-exp(lepsi)
    
    ## Fixed parameters
    c = 20;      # the free virion elimination rate, c, was set to 20 per day
    target0 = 1E8;    # the initial concentration of target cells, T0
    PspectT = 0.001; # initial proportion of specific EBOV CD8 T cells
    k = 4; # the eclipse phase duration, noted 1/k, ranges between 2 and 15h
    de = 0.001; #CD8 T cell perforin+ elimination rate
    dfe = 0.4; # TNF clearance
    dl = 0.4; # IL clearance
    dn = 0.4; # IFN clearance
    sca=2.30258509299;#log(10.0);
    
    ##----------------------------------------- 
    target(0) <- target0
    # infect1(0) <- 0.0
    # infect2(0) <- 0.0
    # refractory(0) <- 0.0
    viremia(0) <- viremia0
    # ifn(0) <- 0.0
    # il6(0) <- 0.0
    # tnf(0) <- 0.0
    nonspecT(0) <- CD8T
    specT(0) <- PspectT*CD8T
    
    ##-----------------------------------------       
    d/dt(target) =-b*target*viremia - fi*target*ifn/(ifn+Tt);
    d/dt(infect1) = b*target*viremia - k*infect1;
    d/dt(infect2) = k*infect1 - d*infect2 - ka*infect2*specT;
    d/dt(refractory) = fi*target*ifn/(ifn+Tt);
    d/dt(viremia) = p*infect2 - c*viremia; ##p*(1-epsi)*infect2 - c*viremia;
    d/dt(ifn) = q*infect2 - dfe*ifn;
    d/dt(il6) = ql*infect2 - dl*il6;
    d/dt(tnf) = qn*infect2 - dn*tnf;
    d/dt(nonspecT) = s - zeta*ifn*nonspecT/(ifn+Te) - de*nonspecT;
    d/dt(specT) = r*specT*(1.0-specT/CD8T)-de*specT;
    
    endpoint=log(viremia)/sca
    endpoint ~ add(add.err)
  })
}

# fit <- nlmixr(vivo, dat, est="nlme")
# fit <- nlmixr(vivo, dat, est="saem")
fit <- nlmixr(vivo, dat, est="saem", saemControl(n.burn=90, n.em=50))

plot(fit)
traceplot(fit)
##
xpdb <- xpose_data_nlmixr(fit)
# 
vpc.ui(fit,n=500)