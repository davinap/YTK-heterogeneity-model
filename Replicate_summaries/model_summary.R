setwd("C:/Users/davin/Documents/PhD/Writing/Heterogeneity/Supp_Data/Tool_promoter_flow/r1-r5_means/KS_all_wneg/")

#load packages

library(gtools)
library(reshape2)
library(gridExtra)
library(openxlsx)
library(readxl)
library(tictoc)
library(scales)
library(directlabels)
library(ggpubr)
library(plotrix)
library(extrafont)
library(grafify)
library(fitdistrplus)
library(flexmix)
library(tidyverse)
library(conflicted)

#Dates of replicate runs:
#r1 = 140623_rep1
#r2 = 140623_rep2
#r3 = 210623_rep1
#r4 = 210623_rep1
#r5 = 300623


####FUNCTIONS AND TEST RUN WITH Y1####
#reads in excel files and turns them into a dataframe
#col_types makes the columns numeric as sometimes they are 'character's instead
#'guess' will guess the type of column - should be character for column 1
r1_y1<- as.data.frame(read_excel("r1_y1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric")))

r2_y1<- as.data.frame(read_excel("r2_y1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_y1<- as.data.frame(read_excel("r3_y1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                              "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_y1<- as.data.frame(read_excel("r4_y1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric")))

r5_y1<- as.data.frame(read_excel("r5_y1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric")))




##View(r1_y1) #good, ignore errors



####ORDER.MUS FUNCTION####
#order and sort each mu into highest to lowest
#cut = sample names and mu1-3 columns only
#cannot be a function because need index to work with the sd and lambda parts

order.mus <- function(x){
  cutm <- select(x, 1:4)
  mltm<- melt(cutm, id.vars="Promoter")
  inx<-order(-mltm$value)
  ordm<-mltm[inx,] #ordered mus
  cuts <- select(x, 1,5:7)
  mlts<- melt(cuts, id.vars="Promoter") 
  ords<- mlts[inx,] #ordered sds
  cutl<-select(x, 1, 8:10)
  mltl <- melt(cutl, id.vars = "Promoter")
  ordl <- mltl[inx,] #ordered lambdas
  
  sam<-rep(c("pFOX2","pHCS1","pHHT1","pTDH2","pHHT2","pCDC19","pRPL28","BP"),1, each=3)
  ordm$Promoter <-factor(ordm$Promoter, levels=unique(sam))
  ordm_a <-arrange(ordm,Promoter) #arrange by column(dataframe of ordered mus to arrange,column)
  ords$Promoter <-factor(ords$Promoter, levels=unique(sam)) #factorise the Promoter column so it can be arranged in order of Promoter name
  ords_a <-arrange(ords,Promoter) #
  ordl$Promoter <-factor(ordl$Promoter, levels=unique(sam))
  ordl_a <-arrange(ordl,Promoter) #
  
  mus<-rep(c("mu1","mu2","mu3"),8) #new dataframe for ordered mus
  sds<-rep(c("sd1","sd2","sd3"),8) #new dataframe for ordered sds
  lams<-rep(c("lam1","lam2","lam3"),8) #new dataframe for ordered lams
  
  dfmus<-data.frame(Promoter=sam, var=mus, value=ordm_a[,3]) #new dataframe to write mus in to
  dfsds<-data.frame(Promoter=sam, var=sds, value=ords_a[,3]) #new dataframe to write sds in to
  dflams<-data.frame(Promoter=sam,var=lams,value=ordl_a[,3]) #new dataframe to write lams in to
  
  bind<- rbind(dfmus, dfsds, dflams) #put all the mus, sds and lambdas together after ordering in desc order
  wide <- reshape(bind, idvar="Promoter", timevar="var", direction="wide") #put mus. sds and lams in their own columns again
  names <- wide %>% rename(mu1 = value.mu1,
                           mu2 = value.mu2,
                           mu3 = value.mu3,
                           sd1 = value.sd1,
                           sd2 = value.sd2,
                           sd3 = value.sd3,
                           lam1 = value.lam1,
                           lam2 = value.lam2,
                           lam3 = value.lam3)
  n <- mutate(names, n = x[,11])
  return(n)
}


y1.1 <- order.mus(r1_y1)
y1.2 <- order.mus(r2_y1)
y1.3 <- order.mus(r3_y1)
y1.4 <- order.mus(r4_y1)
y1.5 <- order.mus(r5_y1)

View(r1_y1)


together_y1<-bind_rows(y1.1,y1.2,y1.3,y1.4,y1.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)


sam2<-c("pFOX2","pHCS1","pHHT1","pTDH2","pHHT2","pCDC19","pRPL28","BP")
together_y1$Promoter<-factor(together_y1$Promoter, levels=sam2)
y1all<- together_y1 %>% arrange(factor(Promoter, levels = sam2))
View(y1all) #good
write.xlsx(y1all, "y1all_ks.xlsx")



####START OF ALL RUNS####
####YPD STAT - Y2####
r1_y2<- as.data.frame(read_excel("r1_y2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r2_y2<- as.data.frame(read_excel("r2_y2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_y2<- as.data.frame(read_excel("r3_y2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_y2<- as.data.frame(read_excel("r4_y2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r5_y2<- as.data.frame(read_excel("r5_y2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

y2.1 <- order.mus(r1_y2)
y2.2 <- order.mus(r2_y2)
y2.3 <- order.mus(r3_y2)
y2.4 <- order.mus(r4_y2)
y2.5 <- order.mus(r5_y2)

together_y2<-bind_rows(y2.1,y2.2,y2.3,y2.4,y2.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)
together_y2$Promoter<-factor(together_y2$Promoter, levels=sam2)
y2all<- together_y2 %>% arrange(factor(Promoter, levels = sam2))
write.xlsx(y2all, "y2all_ks.xlsx")



####YPD10 STAT - x1####
r1_x1<- as.data.frame(read_excel("r1_x1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r2_x1<- as.data.frame(read_excel("r2_x1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_x1<- as.data.frame(read_excel("r3_x1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_x1<- as.data.frame(read_excel("r4_x1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r5_x1<- as.data.frame(read_excel("r5_x1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

x1.1 <- order.mus(r1_x1)
x1.2 <- order.mus(r2_x1)
x1.3 <- order.mus(r3_x1)
x1.4 <- order.mus(r4_x1)
x1.5 <- order.mus(r5_x1)

together_x1<-bind_rows(x1.1,x1.2,x1.3,x1.4,x1.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)
together_x1$Promoter<-factor(together_x1$Promoter, levels=sam2)
x1all<- together_x1 %>% arrange(factor(Promoter, levels = sam2))
write.xlsx(x1all, "x1all_ks.xlsx")



####YPD10 STAT - x2####
r1_x2<- as.data.frame(read_excel("r1_x2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r2_x2<- as.data.frame(read_excel("r2_x2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_x2<- as.data.frame(read_excel("r3_x2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_x2<- as.data.frame(read_excel("r4_x2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r5_x2<- as.data.frame(read_excel("r5_x2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

x2.1 <- order.mus(r1_x2)
x2.2 <- order.mus(r2_x2)
x2.3 <- order.mus(r3_x2)
x2.4 <- order.mus(r4_x2)
x2.5 <- order.mus(r5_x2)

together_x2<-bind_rows(x2.1,x2.2,x2.3,x2.4,x2.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)
together_x2$Promoter<-factor(together_x2$Promoter, levels=sam2)
x2all<- together_x2 %>% arrange(factor(Promoter, levels = sam2))
write.xlsx(x2all, "x2all_ks.xlsx")



####YNB EXP - m1####
r1_m1<- as.data.frame(read_excel("r1_m1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r2_m1<- as.data.frame(read_excel("r2_m1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_m1<- as.data.frame(read_excel("r3_m1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_m1<- as.data.frame(read_excel("r4_m1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r5_m1<- as.data.frame(read_excel("r5_m1.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

m1.1 <- order.mus(r1_m1)
m1.2 <- order.mus(r2_m1)
m1.3 <- order.mus(r3_m1)
m1.4 <- order.mus(r4_m1)
m1.5 <- order.mus(r5_m1)

together_m1<-bind_rows(m1.1,m1.2,m1.3,m1.4,m1.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)
together_m1$Promoter<-factor(together_m1$Promoter, levels=sam2)
m1all<- together_m1 %>% arrange(factor(Promoter, levels = sam2))
write.xlsx(m1all, "m1all_ks.xlsx")



####YNB STAT - m2####
r1_m2<- as.data.frame(read_excel("r1_m2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r2_m2<- as.data.frame(read_excel("r2_m2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r3_m2<- as.data.frame(read_excel("r3_m2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r4_m2<- as.data.frame(read_excel("r4_m2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

r5_m2<- as.data.frame(read_excel("r5_m2.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric")))

m2.1 <- order.mus(r1_m2)
m2.2 <- order.mus(r2_m2)
m2.3 <- order.mus(r3_m2)
m2.4 <- order.mus(r4_m2)
m2.5 <- order.mus(r5_m2)

together_m2<-bind_rows(m2.1,m2.2,m2.3,m2.4,m2.5) %>% group_by(Promoter) %>% summarise_all(mean, na.rm=TRUE)
together_m2$Promoter<-factor(together_m2$Promoter, levels=sam2)
m2all<- together_m2 %>% arrange(factor(Promoter, levels = sam2))
write.xlsx(m2all, "m2all_ks.xlsx")






####long tabling####
#making long versions of table for verthis script
#want sample/var/value
melt(y1all) %>% filter(str_detect(variable, "mu")) #good


####JEFLABS example####
mo1 <- FLXMRglm(family = "gaussian")
mo2 <- FLXMRglm(family = "gaussian")
mo3 <- FLXMRglm(family = "gaussian")

#function for making plots with area relative to population size
plot_mix_comps <- function(x, mu, sigma, lam) {lam * dnorm(x, mu, sigma)}

####tgg####
tgg <-theme(axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(family="Calibri Light", size=15),
            axis.ticks.y = element_line(colour="gray", linewidth =0.5),
            axis.ticks.x = element_line(colour="gray", linewidth =0.5),
            axis.line.y = element_line(colour="gray", linewidth =0.5),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank(),
            plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))




#####colour test####
make.3dist.ks <- function(m1,m2,m3,sd1,sd2,sd3,n1,n2,n3,f,neg){
  #read in raw data for making raw flow data plot
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  
  ctrl <- read.csv(neg, header = TRUE, stringsAsFactors = FALSE) #to plot negative control (BP strain)
  ctrlpos <- filter(ctrl, ctrl[1]>0) #remove rows containing negative values which would break log transform
  logctrl <- log10(ctrlpos)
  
  m1 <- m1
  m2 <- m2
  m3 <- m3
  sd1 <- sd1
  sd2 <- sd2 
  sd3 <- sd3
  p1 <- (n1/(n1+n2+n3))*2000
  p2 <- (n2/(n1+n2+n3))*2000
  p3 <- (n3/(n1+n2+n3))*2000
  
  set.seed(1)
  a <- rnorm(n=p1, mean=m1, sd=sd1)
  b <- rnorm(n=p2, mean=m2, sd=sd2)
  c <- rnorm(n=p3, mean=m3, sd=sd3)
  
  x <- c(a,b,c)
  class <- c(rep('a', p1), rep('b', p2),rep('c', p3))
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  ff<-flexmix(data$x~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  
  #plot
  p <- ggplot(data, aes(x=x)) + 
    geom_density(data = logged, aes(x=`Comp.BL1.H....BL1.H`), color=alpha("black",0.8)) +
    stat_function(fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#7D4600", fill="#7D4600",alpha=0.5, geom="polygon") + #brown
    stat_function(fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#565264",fill="#565264",alpha=0.5,geom="polygon") + #blue-grey
    stat_function(fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#69995D",fill="#69995D",alpha=0.5,geom="polygon") + #green
    ##  geom_vline(xintercept = (pr), linetype = 2, colour = "#00d400ff", size=1.25) +
    ##geom_density(data=ctrl, aes(x=neg), color="#dddddd",size=2,fill="#dddddd", alpha=0.3) + #plot negative control
    geom_density(data = logctrl, aes(x=`Comp.BL1.H....BL1.H`), color="#A0A0A0",linetype = "dashed")+

  
    geom_hline(yintercept = 0, colour = "gray") +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0,3.75), breaks = c(1,2,3), expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    tgg +
    xlab("")
  return(p)
}

make.1dist.ks <- function(m1, sd1, n1, f,neg){
  
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  
  ctrl <- read.csv(neg, header = TRUE, stringsAsFactors = FALSE) #to plot negative control (BP strain)
  ctrlpos <- filter(ctrl, ctrl[1]>0) #remove rows containing negative values which would break log transform
  logctrl <- log10(ctrlpos)
  
  a <- rnorm(n=n1, mean=m1, sd=sd1)
  x <- a
  class <- c(rep('a', n1))
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  p <- ggplot(data, aes(x=x)) + 
    geom_density(data = logged, aes(x=`Comp.BL1.H....BL1.H`), color=alpha("black",0.8)) +
    stat_function(fun=plot_mix_comps,
                  args=list(mu=m1, 
                            sigma=sd1,
                            lam=n1/n1),
                  geom = "area",
                  color="#d6cfcb",
                  fill="#d6cfcb", alpha=0.65) +
   ## geom_vline(xintercept = pr, linetype = 2, colour = "#00d400ff", size=1.25) +
    geom_density(data = logctrl, aes(x=`Comp.BL1.H....BL1.H`), color="#A0A0A0",linetype = "dashed")+
    
    geom_hline(yintercept = 0, colour = "gray") +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0,3.75), breaks = c(1,2,3), expand = expansion(mult = c(0, 0.05))) +
    xlab("")+
    ylab("")+
    tgg
  return(p)
}




####making all x1 graphs####

#for ypd10 rep5
x1df <- as.data.frame(x1all)



p1<- make.1dist.ks(x1df[1,2], x1df[1,5], x1df[1,8], "export_A1_live3.csv","export_A8_live3.csv") #pFOX2
p2<- make.1dist.ks(x1df[2,2], x1df[2,5], x1df[2,8], "export_A2_live3.csv","export_A8_live3.csv") #pHCS1
p3<- make.3dist.ks(x1df[3,2], x1df[3,3], x1df[3,4],
                   x1df[3,5], x1df[3,6],x1df[3,7],
                   x1df[3,8], x1df[3,9], x1df[3,10], "export_A3_live3.csv","export_A8_live3.csv") #pHHT1
p4<- make.1dist.ks(x1df[4,2], x1df[4,5], x1df[4,8], "export_A4_live3.csv","export_A8_live3.csv") #pTDH2
p5<- make.3dist.ks(x1df[5,2], x1df[5,3], x1df[5,4],
                   x1df[5,5], x1df[5,6],x1df[5,7],
                   x1df[5,8], x1df[5,9], x1df[5,10], "export_A5_live3.csv","export_A8_live3.csv") #pHHT2
p6<- make.1dist.ks(x1df[6,2], x1df[6,5], x1df[6,8], "export_A6_live3.csv","export_A8_live3.csv") #pCDC19
p7<- make.1dist.ks(x1df[7,2], x1df[7,5], x1df[7,8], "export_A7_live3.csv","export_A8_live3.csv") #pRPL28
  ##theme(axis.text.x = element_text(family="Calibri Light", size=7)) 

p3
p5
p1
p2


test1 <- ggarrange(p2,p3,p4,p5,p6,p7, 
                   ncol = 1, nrow=6)

test1

ypdxexp_plots<-ggarrange(p2,p3,p4,p5,p6,p7, 
                         ncol = 1, nrow=6)
ypdxexp_plots



####NEXT - 290823 - x2####
#for ypd10 stat rep5
x2df <- as.data.frame(x2all)
##View(x2df)

q1<- make.1dist.ks(x2df[1,2], x2df[1,5], x2df[1,8], "export_B1_live4.csv","export_B8_live4.csv") 
q2<- make.1dist.ks(x2df[2,2], x2df[2,5], x2df[2,8], "export_B2_live4.csv","export_B8_live4.csv") 
q3<- make.1dist.ks(x2df[3,2], x2df[3,5], x2df[3,8], "export_B3_live4.csv","export_B8_live4.csv")
q4<- make.1dist.ks(x2df[4,2], x2df[4,5], x2df[4,8], "export_B4_live4.csv","export_B8_live4.csv")
q5<- make.3dist.ks(x2df[5,2], x2df[5,3], x2df[5,4],
                   x2df[5,5], x2df[5,6], x2df[5,7],
                   x2df[5,8], x2df[5,9], x2df[5,10],"export_B5_live4.csv","export_B8_live4.csv")
q6<- make.1dist.ks(x2df[6,2], x2df[6,5], x2df[6,8], "export_B6_live4.csv","export_B8_live4.csv")
q7<- make.1dist.ks(x2df[7,2], x2df[7,5], x2df[7,8], "export_B7_live4.csv","export_B8_live4.csv") ##+ 
 ##theme(axis.text.x = element_text(family="Calibri Light", size=7))




test2 <- ggarrange(q2,q3,q4,q5,q6,q7,
          ncol = 1, nrow = 6)
##ggsave('test.svg', width=490, height=375, units='px')

ypdx_staplots <-ggarrange(q2,q3,q4,q5,q6,q7,
                          ncol = 1, nrow = 6)

ypdxexp_plots
ypdx_staplots 

#for ypd exp 
#raw data = rep1
y1df <- as.data.frame(y1all)

#renamed plots to s1-7 instead of r1-7 as 'r1' has been used in script previously
s1<- make.1dist.ks(y1df[1,2], y1df[1,5], y1df[1,8], "export_A1_live1.csv","export_A8_live1.csv") 
s2<- make.1dist.ks(y1df[2,2], y1df[2,5], y1df[2,8], "export_A2_live1.csv","export_A8_live1.csv") 
s3<- make.3dist.ks(y1df[3,2], y1df[3,3], y1df[3,4],
                   y1df[3,5], y1df[3,6],y1df[3,7],
                   y1df[3,8], y1df[3,9], y1df[3,10], "export_A3_live1.csv","export_A8_live1.csv")
s4<- make.1dist.ks(y1df[4,2], y1df[4,5], y1df[4,8], "export_A4_live1.csv","export_A8_live1.csv")
s5<- make.3dist.ks(y1df[5,2], y1df[5,3], y1df[5,4],
                   y1df[5,5], y1df[5,6],y1df[5,7],
                   y1df[5,8], y1df[5,9], y1df[5,10], "export_A5_live1.csv","export_A8_live1.csv")
s6<- make.1dist.ks(y1df[6,2], y1df[6,5], y1df[6,8], "export_A6_live1.csv","export_A8_live1.csv")
s7<- make.1dist.ks(y1df[7,2], y1df[7,5], y1df[7,8], "export_A7_live1.csv","export_A8_live1.csv") ##+ 
 ## theme(axis.text.x = element_text(family="Calibri Light", size=7))


ypd_expplots <- ggarrange(s2,s3,s4,s5,s6,s7,
          ncol = 1, nrow = 6)

ypd_expplots


#for ypd stationary
#raw data = rep1 140623
y2df <- as.data.frame(y2all)

t1<- make.1dist.ks(y2df[1,2], y2df[1,5], y2df[1,8], "export_A1_live2.csv","export_A8_live2.csv") 
t2<- make.1dist.ks(y2df[2,2], y2df[2,5], y2df[2,8], "export_A2_live2.csv","export_A8_live2.csv") 
t3<- make.1dist.ks(y2df[3,2], y2df[3,5], y2df[3,8], "export_A3_live2.csv","export_A8_live2.csv")
t4<- make.1dist.ks(y2df[4,2], y2df[4,5], y2df[4,8], "export_A4_live2.csv","export_A8_live2.csv")
t5<- make.3dist.ks(y2df[5,2], y2df[5,3], y2df[5,4],
                   y2df[5,5], y2df[5,6], y2df[5,7],
                   y2df[5,8], y2df[5,9], y2df[5,10],"export_A5_live2.csv","export_A8_live2.csv")
t6<- make.1dist.ks(y2df[6,2], y2df[6,5], y2df[6,8], "export_A6_live2.csv","export_A8_live2.csv")
t7<- make.1dist.ks(y2df[7,2], y2df[7,5], y2df[7,8], "export_A7_live2.csv","export_A8_live2.csv") ##+ 
 ## theme(axis.text.x = element_text(family="Calibri Light", size=7))


ypd_staplots <- ggarrange(t2,t3,t4,t5,t6,t7,
                   ncol = 1, nrow = 6)

ypd_staplots


#for ynb exp 
#raw data = rep5 300623
m1df <- as.data.frame(m1all)

#renamed plots to s1-7 instead of r1-7 as 'r1' has been used in script previously
u1<- make.1dist.ks(m1df[1,2], m1df[1,5], m1df[1,8], "export_B1_live5.csv","export_B8_live5.csv") 
u2<- make.1dist.ks(m1df[2,2], m1df[2,5], m1df[2,8], "export_B2_live5.csv","export_B8_live5.csv") 
u3<- make.3dist.ks(m1df[3,2], m1df[3,3], m1df[3,4],
                   m1df[3,5], m1df[3,6],m1df[3,7],
                   m1df[3,8], m1df[3,9], m1df[3,10], "export_B3_live5.csv","export_B8_live5.csv")
u4<- make.1dist.ks(m1df[4,2], m1df[4,5], m1df[4,8], "export_B4_live5.csv","export_B8_live5.csv")
u5<- make.3dist.ks(m1df[5,2], m1df[5,3], m1df[5,4],
                   m1df[5,5], m1df[5,6],m1df[5,7],
                   m1df[5,8], m1df[5,9], m1df[5,10], "export_B5_live5.csv","export_B8_live5.csv")
u6<- make.1dist.ks(m1df[6,2], m1df[6,5], m1df[6,8], "export_B6_live5.csv","export_B8_live5.csv")
u7<- make.1dist.ks(m1df[7,2], m1df[7,5], m1df[7,8], "export_B7_live5.csv","export_B8_live5.csv") ##+ 
  ##theme(axis.text.x = element_text(family="Calibri Light", size=7))


ynb_expplot <- ggarrange(u2,u3,u4,u5,u6,u7,
                   ncol = 1, nrow = 6)
ynb_expplot


#for ynb exp 
#raw data = rep5 300623
m2df <- as.data.frame(m2all)

#renamed plots to s1-7 instead of r1-7 as 'r1' has been used in script previously
v1<- make.1dist.ks(m2df[1,2], m2df[1,5], m2df[1,8], "export_C1_live6.csv","export_C8_live6.csv") 
v2<- make.1dist.ks(m2df[2,2], m2df[2,5], m2df[2,8], "export_C2_live6.csv","export_C8_live6.csv") 
v3<- make.3dist.ks(m2df[3,2], m2df[3,3], m2df[3,4],
                   m2df[3,5], m2df[3,6],m2df[3,7],
                   m2df[3,8], m2df[3,9], m2df[3,10], "export_C3_live6.csv","export_C8_live6.csv")
v4<- make.1dist.ks(m2df[4,2], m2df[4,5], m2df[4,8], "export_C4_live6.csv","export_C8_live6.csv")
v5<- make.3dist.ks(m2df[5,2], m2df[5,3], m2df[5,4],
                   m2df[5,5], m2df[5,6],m2df[5,7],
                   m2df[5,8], m2df[5,9], m2df[5,10], "export_C5_live6.csv","export_C8_live6.csv")
v6<- make.1dist.ks(m2df[6,2], m2df[6,5], m2df[6,8], "export_C6_live6.csv","export_C8_live6.csv")
v7<- make.1dist.ks(m2df[7,2], m2df[7,5], m2df[7,8], "export_C7_live6.csv","export_C8_live6.csv") ##+ 
  ##theme(axis.text.x = element_text(family="Calibri Light", size=7))


ynb_staplots <- ggarrange(v2,v3,v4,v5,v6,v7,
                   ncol = 1, nrow = 6)
ynb_staplots


#put all plots together:
ggarrange(ypdxexp_plots,ypdx_staplots,ypd_expplots,ypd_staplots,ynb_expplot,ynb_staplots,
          ncol=2,
          nrow=3)

ypdxexp_plots
ypdx_staplots
ypd_expplots
ypd_staplots
ynb_expplot
ynb_staplots




#plot hht1 separately with make clear function with a smaller y axis for magnified plot
v3
u3


tgg2 <-theme(axis.title.y = element_blank(),
             axis.title.x = element_blank(),
             axis.text.y = element_text(family="Calibri Light", size=9),
             axis.text.x = element_text(family="Calibri Light", size=10),
             axis.ticks.y = element_line(colour="gray", linewidth =0.5),
             axis.ticks.x = element_line(colour="gray", linewidth =0.5),
             axis.line.y = element_line(colour="gray", linewidth =0.5),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             plot.margin = unit(c(2, 5.5, 2, 5.5), "pt"))


make.3dist.clear <- function(m1,m2,m3,sd1,sd2,sd3,n1,n2,n3,f){
  #read in raw data for making raw flow data plot
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  
  m1 <- m1
  m2 <- m2
  m3 <- m3
  sd1 <- sd1
  sd2 <- sd2 
  sd3 <- sd3
  p1 <- (n1/(n1+n2+n3))*2000
  p2 <- (n2/(n1+n2+n3))*2000
  p3 <- (n3/(n1+n2+n3))*2000
  
  set.seed(1)
  a <- rnorm(n=p1, mean=m1, sd=sd1)
  b <- rnorm(n=p2, mean=m2, sd=sd2)
  c <- rnorm(n=p3, mean=m3, sd=sd3)
  
  x <- c(a,b,c) #good
  class <- c(rep('a', p1), rep('b', p2),rep('c', p3)) #good, 2000 entries with a/b/c proportions 
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  ff<-flexmix(data$x~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  
  #plot
  p <- ggplot(data, aes(x=x)) + 
    geom_density(data = logged, aes(x=`Comp.BL1.H....BL1.H`), color="#A0A0A0") +
    stat_function(fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#7D4600", fill="#7D4600",alpha=0.5, geom="polygon") + #brown
    stat_function(fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#565264",fill="#565264",alpha=0.5,geom="polygon") + #grey
    stat_function(fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#69995D",fill="#69995D",alpha=0.5,geom="polygon") + #green
    ##  geom_vline(xintercept = (pr), linetype = 2, colour = "#00d400ff", size=1.25) +
    geom_hline(yintercept = 0, colour = "gray") +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0,1.8), breaks = c(0.5,1,1.5), expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    tgg2 +
    xlab("")
  ##return(p)
  return(p)
}


#clear plots
clear_hht1<- make.3dist.clear(m2df[3,2], m2df[3,3], m2df[3,4],
                   m2df[3,5], m2df[3,6],m2df[3,7],
                   m2df[3,8], m2df[3,9], m2df[3,10], "export_C3_live6.csv")
##View(clear_hht1) 

clear_hht2<- make.3dist.clear(m2df[5,2], m2df[5,3], m2df[5,4],
                   m2df[5,5], m2df[5,6],m2df[5,7],
                   m2df[5,8], m2df[5,9], m2df[5,10], "export_C5_live6.csv")

ggarrange(clear_hht1, clear_hht2,
          ncol=1,
          nrow=2)




