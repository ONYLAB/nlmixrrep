library(RxODE)
library(nlmixr)
library(xpose)
library(xpose.nlmixr)
library(ggplot2)
library(pdftools)

rm(list=ls(all=TRUE)) 
graphics.off()

dat <- read.csv("EBOV.csv", head=TRUE)
str(dat)
ggplot(dat %>% filter(dat$DVID=='V0'), aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 Viremia (cp/mL)") + labs(title="Viremia", subtitle="Viremia vs. time by individual")
ggplot(dat %>% filter(dat$DVID=='CD8perforin'), aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 CD8 Perforin+ (cell/mL/1000)") + labs(title="CD8", subtitle="CD8+ vs. time by individual")
ggplot(dat %>% filter(dat$DVID=='IFNa'), aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 IFNa (pg/mL)") + labs(title="IFNa", subtitle="IFNa vs. time by individual")
ggplot(dat %>% filter(dat$DVID=='IL6'), aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 IL6 (pg/mL)") + labs(title="IL6", subtitle="IL6 vs. time by individual")
ggplot(dat %>% filter(dat$DVID=='TNFa'), aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (d)") + scale_y_continuous("log10 TNFa (pg/mL)") + labs(title="TNFa", subtitle="TNFa vs. time by individual")

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

    eta.d ~ .1 #log(0.0004064256) #(0.224*(9/100))^2 #0.028561
    eta.p ~ 15 #(4.15E4*(71/100))^2  #3.0976
    eta.V0 ~ 1E-7 # (7.9E-5*(11/100))^2  #2.25
    # eta.b ~ -23.26157326346153 ##!
    # eta.q ~ -4.9062752787720125 ##!
    eta.Tt ~ 0.2 #(1.73*(35/100))^2 #1.1025
    # eta.fi ~ 0.9820784724121581 ##!
    eta.E0 ~ 15 # (36900*(30/100))^2 #0.600625
    # eta.zeta ~ -0.7874578600311866 ##!
    # eta.r ~ -1.0847093834991182 ##!
    # eta.s ~ 8.022896869601457 ##!
    # eta.Te ~ -7.338538195074591 ##!
    # eta.ka ~ -10.780557571257003 ##!
    eta.ql ~ .01 # (0.0097*(65/100))^2 #0.5776
    eta.qn ~ .01 #(0.0046*(69/100))^2 #1.0609

    add.err <- 1 #
  })
  model({
    
    d <- exp(ld + eta.d)
    p <- exp(lp + eta.p)
    viremia0 <- exp(lV0 + eta.V0)    
    b<-exp(lb)
    q<-exp(lq)
    Tt <- exp(lTt + eta.Tt)
    fi<-exp(lfi)
    CD8T <- exp(lE0 + eta.E0)
    zeta<-exp(lzeta)
    r<-exp(lr)
    s<-exp(ls)
    Te<-exp(lTe)
    ka<-exp(lka)
    ql <- exp(lql + eta.ql)
    qn <- exp(lqn + eta.qn)
    #epsi<-exp(lepsi)

    ## Fixed parameters
    c = 20;           # the free virion elimination rate, c, was set to 20 per day
    target0 = 1E8;    # the initial concentration of target cells, T0
    PspectT = 0.001;  # initial proportion of specific EBOV CD8 T cells
    k = 4;            # the eclipse phase duration, noted 1/k, ranges between 2 and 15h
    de = 0.001;       # CD8 T cell perforin+ elimination rate
    dfe = 0.4;        # TNF clearance
    dl = 0.4;         # IL clearance
    dn = 0.4;         # IFN clearance
    sca= 2.30258509299;    # log(10.0);
    epsi=0;

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
    d/dt(viremia) = p*(1-epsi)*infect2 - c*viremia;
    d/dt(ifn) = q*infect2 - dfe*ifn;
    d/dt(il6) = ql*infect2 - dl*il6;
    d/dt(tnf) = qn*infect2 - dn*tnf;
    d/dt(nonspecT) = s - zeta*ifn*nonspecT/(ifn+Te) - de*nonspecT;
    d/dt(specT) = r*specT*(1.0-specT/CD8T)-de*specT;
    
    Vir0=log(viremia)/sca
    Vir0 ~ add(add.err) | viremia

    CD8perforin=log((nonspecT+specT)/1000)/sca
    CD8perforin ~ add(add.err) | (nonspecT+specT)

    IFNa=log(ifn)/sca
    IFNa ~ add(add.err) | ifn

    IL6=log(il6)/sca
    IL6 ~ add(add.err) | il6

    TNFa=log(tnf)/sca
    TNFa ~ add(add.err) | tnf
    
  })
}

# fit <- nlmixr(vivo, dat, est="nlme")
# fit <- nlmixr(vivo, dat, est="foce")
fit <- nlmixr(vivo, dat, est="saem")
# fit <- nlmixr(vivo, dat, est="saem", saemControl(n.burn=90, n.em=50))

##
xpdb <- xpose_data_nlmixr(fit)

##Goodness-of-fit plots / Analysis
pdf(file="GOF plots case EBOV NHP.pdf", height=9, width=9, paper="letter")
plot(augPred(fit))
tp <- traceplot(fit)
tp + ylab("Log(parameter)") + xlab("Iteration")
dv_vs_pred(xpdb) + ylab("Observed Viremia (log10 cp/mL)") + xlab("Population Predicted Viremia (log10 cp/mL)")
dv_vs_ipred(xpdb) + ylab("Observed Viremia (log10cp/mL)") + xlab("Individual Predicted Viremia (log10 cp/mL)")
res_vs_pred(xpdb) + ylab("Conditional Weighted Residuals") + xlab("Population Predicted Viremia (log10 cp/mL)")
res_vs_idv(xpdb) + ylab("Conditional Weighted Residuals") + xlab("Time (d)")
prm_vs_iteration(xpdb)
absval_res_vs_idv(xpdb, res = 'IWRES') + ylab("Individual Weighted Residuals") + xlab("Time (d)")
absval_res_vs_pred(xpdb, res = 'IWRES') + ylab("Individual Weighted Residuals") + xlab("Population Predicted Viremia (log10 cp/mL)")
ind_plots(xpdb, nrow=3, ncol=4) + ylab("Predicted and Observed Viremia (log10 cp/mL)") + xlab("Time (d)")
res_distrib(xpdb) + ylab("Density") + xlab("Conditional Weighted Residuals")
##Visual Predictive Checks
# vpc.ui(fit,n=500)
# vpc.ui(fit,n=500,stratify=NULL, show=list(obs_dv=T), bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50), ylab = "Viremia (cp/mL)", xlab = "Time (days)")
# vpc.ui(fit,n=500, show=list(obs_dv=T), bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50), ylab = "Viremia (cp/mL)", xlab = "Time (days)")
dev.off()

 
 
# # 
# # ev <- eventTable(amount.units='ug', time.units='Days')
# # ev$add.sampling(seq(0,30,length.out=100)) #ev$add.sampling(0:30)
# # 
# # x <- solve(mod1,theta, ev, inits);
# # 
# # library(ggplot2)
# # x <- as.data.frame(x)
# # ggplot(x,aes(time,log10(V))) + geom_line() + ylab(expression("Viral Load (log10 cp mL"*{}^{-1}*")")) + xlab("Time") + ylim(-2,10) #Viral load
# # ggplot(x,aes(time,1E-3*(E1+E2))) + geom_line() + ylab(expression(paste("CD8 perforin (cells ", mu, "L"*{}^{-1}*")"))) + xlab("Time") + ylim(0,150) #CD8+Perforin+ concentration
# # 
# # #Cytokines
# # ggplot(x,aes(time,log10(N))) + geom_line() + ylab(expression(paste("IFN", alpha, " (log10 pg mL"*{}^{-1}*")"))) + xlab("Time") + ylim(-2,4)
# # ggplot(x,aes(time,log10(F))) + geom_line() + ylab(expression(paste("TNF", alpha, " (log10 pg mL"*{}^{-1}*")"))) + xlab("Time") + ylim(0,4)
# # ggplot(x,aes(time,log10(L))) + geom_line() + ylab(expression(paste("IL6", " (log10 pg mL"*{}^{-1}*")"))) + xlab("Time") + ylim(0,4) 
# # 
