## Respiratory physiology of diving caimans ##
## by Bautista and Damsgaard ##

# Set working directory
setwd("/Users/au231308/Dropbox/Projects/Crocodile RBC/")


# Load packages
Packages <- c("RColorBrewer","ggplot2", "gridExtra","readxl","nlme","lme4","multcomp","multcompView","cowplot")
lapply(Packages, library, character.only = TRUE) # load all required packages

# Colours
cols                     <- brewer.pal(9,"Set1")
cols<-cols[c(1,2,3,4,5,7,8,9)] # remove yellow

## Custom functions
stderr                   <- function(x) sd(x)/sqrt(length(x))                                      # Standard error of mean
lower                    <- function(x) mean(x)-stderr(x)                                          # Mean - standard error of mean
upper                    <- function(x) mean(x)+stderr(x)                                          # Mean + standard error of mean


# Constants
aco2 <-0.03762409 # CO2 solubility at 28 degrees in mmol mmHg-1
ao2 <-0.001591684 # O2 solubility at 28 degrees in mmol mmHg-1 


#Calibration constants for [O2] determination
CV <-1.3961 # chamber volume in ml
SV <-0.01954 # syringe volume in ml
DF <-(CV-SV)/CV # dilution factor
CC <-0.090929256 # calibration constant mmol l-1 mmHg-1


# Load data
data<-read_excel(path = "./Data/raw data.xlsx", sheet = "raw data") # Import data
data<-as.data.frame(data) # save as data frame
data <- data[data$treatment!="R",] # Dont include the two recovery samples
data$pHi
# Calculated oxygen and hemoglobin values
data$po2.28 <- 0.2095*(data$pbar-28.3) # atmospheric PO2 at 100% humidity and 28 oC
data$po2.40 <- 0.2095*(data$pbar-55.3) # atmospheric PO2 at 100% humidity and 40 oC (ie PO2 in the Tucker chamber)
data$po2 <- data$po2*data$po2.28/155 # # Correct arterial PO2 for correct barometric pressure
data$hb   <-data$hct/100*25 # whole blood hemoglobin concentration based on Hct and MCHC 
data$po2.1.1<-data$po2.1.1*data$po2.40/155 # Correct Tucker PO2 for correct barometric pressure
data$po2.1.2<-data$po2.1.2*data$po2.40/155 # Correct Tucker PO2 for correct barometric pressure
data$po2.2.1<-data$po2.2.1*data$po2.40/155 # Correct Tucker PO2 for correct barometric pressure
data$po2.2.2<-data$po2.2.2*data$po2.40/155 # Correct Tucker PO2 for correct barometric pressure
data$dpo2<-rowMeans(cbind(data$po2.1.2-data$po2.1.1*DF,data$po2.2.2-data$po2.2.1*DF), na.rm = T)
data$hbo2   <-data$dpo2*CC-ao2*data$po2
data$sat   <-data$hbo2/data$hb*100
data$lac<-data$`Lactate (mM)`


# Calculated acid/base parameters
#data$phi<- 0.609485885*data$pha+2.707293323 # intracellular RBC pH based on Donnan distribution ratio from Jensen et al 1998
data$pk   <- -0.0817*data$pha+6.7818 # CO2 hydration constant from Jensen et al 1998
data$pco2 <- data$co2p/(aco2*(1+10^(data$pha-data$pk))) # PCO2 calculated from Henderson Hasselbalch equation. 
data$hco3p <-data$co2p-aco2*data$pco2 # plasma HCO3 

  # calculate pH and sat-corrected HCO3i.app
  data$r<-rep(NA,dim(data)[1])
  for (i in 1:dim(data)[1]){
    sat<-data$sat[i]/100
    ph<-data$pha[i]
    desat<-1-sat
    int<-sat*13.94622439+desat*5.596302496;int # calculated the Hb-O2 sat-weighted intercept of r vs pHa from Jensen et al 1998
    slo<-sat*(-1.680718031)+desat*(-0.507047282);slo # calculated the Hb-O2 sat-weighted intercept of r vs pHa from Jensen et al 1998
    data$r[i]<-int+slo*ph;data$r[i]
    data$hco3i.app[i]<-data$r[i]*data$hco3p[i] # Calculate apparent RBC HCO3 based on R and plasma HCO3. 
  }

data$hco3i <- data$cl.rbc/data$cl.p*data$hco3p # Intracellular HCO3 based on the Donnan distribution ratio of Cl-
data$hco3hb<-data$hco3i.app-data$hco3i



## PLOT DATA ##
# parameters to be plotted
para<-c("pco2", "po2","hco3p","pha","pHi","hct","osm","hco3i","sat","lac","cl.rbc","cl.p","hco3i.app","hco3hb")


# Corresponding axis labels
label                    <- c(expression(bold("P"[a]*"CO"["2"]*" (mmHg)")),
                              expression(bold("P"[a]*"O"["2"]*" (mmHg)")),
                              expression(bold("[HCO"[3]^-{}*"]"[p]*" (mmol l"^"-1"*")")),
                              expression(bold("pH"[a])),
                              expression(bold("pH"[i])),
                              expression(bold("[Hct] (%)")),
                              expression(bold("Osmolality (mOsm kg"^"-1"*")")),
                              expression(bold("[HCO"[3]^-{}*"]"["i, free"]*" (mmol l"^"-1"*")")),
                              expression(bold("Hb-O"[2]*" sat (%)")),
                              expression(bold("[lactate] (mmol l"^"-1"*")")),
                              expression(bold("[Cl"^"-"*"]"[i]*" (mmol l"^"-1"*")")),
                              expression(bold("[Cl"^"-"*"]"[p]*" (mmol l"^"-1"*")")),
                              expression(bold("Total RBC [HCO"[3]^-{}*"]"[i]*" (mmol l"^"-1"*")")),
                              expression(bold("[Hb-HCO"[3]^-{}*"] (mmol l"^"-1"*")")))

# Generate a function that runs stats and plots each parameter as a function of salinity

plot<-function (i,xtext) {
  xtext                      <- ifelse(missing(xtext),  F,T)
  
  # Calculate means and standard error
  df<-data
  
  a = cbind(aggregate(as.formula(paste(as.name(para[i]),"~treatment")),df,mean),                         # Data frame with mean and standard error over salinities
            aggregate(as.formula(paste(as.name(para[i]),"~treatment")),df,stderr)[,2],
            aggregate(as.formula(paste(as.name(para[i]),"~treatment")),df,length)[,2])
  colnames(a)<-c("Treatment","Mean","stderr","n")
  a$Treatment<-as.factor(a$Treatment)
  
  
  # Run mixed model anova and do pair wise comparisons
  model<-lmer(paste(as.name(para[i]),"~treatment+(1|animalID)"),df);model
  fit<-summary(glht(model, linfct = mcp(treatment = "Tukey")), test = adjusted("holm"));fit
  p<-fit$test$pvalues;p
  names(p)<-sub(" - ", "-", names(p));p
  p<-multcompLetters(p);p
  p<-p$Letters;p
  p<-c(p[3],p[1:2]);p
  data[data$`animal number`==i,][para[2]]
  data.frame(
    data[data$`animal number`==i,][para[2]],
    data[data$`animal number`==i,]["treatment"])
  
  write.table(a,file = paste("./Data/Summarized data/",para[i],"_means.txt",sep = ""))
  write.table(fit$test$pvalues,file = paste("./Data/Summarized data/",para[i],"_MMAOV_posthoc_p-value.txt",sep = ""))
  write.table(fit$test$tstat,file = paste("./Data/Summarized data/",para[i],"_MMAOV_posthoc_tstatistic.txt",sep = ""))
  # Plot the data
  
  ggplot() +
    
    # Add individual data points and connecting lines
    lapply(
      X = 1:8, 
      FUN = function(k){
        temp.df<-data[data$`animal number`==k,]
        temp.df<-
          data.frame(
            data[data$`animal number`==k,][para[i]],
            data[data$`animal number`==k,]["treatment"])
        colnames(temp.df)<-c("value","treatment")
        geom_point(data = temp.df, mapping = aes(x = treatment, y = value), color =cols[k], alpha = 0.5)
      })+ 
    
    lapply(
      X = 1:8, 
      FUN = function(k){
        temp.df<-data[data$`animal number`==k,]
        temp.df<-
          data.frame(
            data[data$`animal number`==k,][para[i]],
            data[data$`animal number`==k,]["treatment"])
        colnames(temp.df)<-c("value","treatment")
        geom_line(data = temp.df, mapping = aes(x = treatment, y = value,group = 1), color =cols[k], alpha = 0.5)
      }) +
    
    
    # Add means standard error and connecting lines
    geom_point(data=a, aes(x=Treatment, y=Mean))+
    geom_path(data=a, aes(x=Treatment, y=Mean, group = 1))+
    geom_errorbar(data = a, aes(x = Treatment, y = Mean, ymin=Mean-stderr, ymax=Mean+stderr),width=0.5,lwd=0.5) +
    
    # Add pairwise comparison labels
    geom_text(
      aes(label=p),
      vjust=-0.5, color="black",
      x = a$Treatment,
      y=as.numeric(a$Mean+a$stderr), 
      size=2.8
    )+
    
    # Change x-axis text
    scale_x_discrete(breaks = c("C1","D1","D2"), name = "",
                     labels = c("Pre-dive","18 min dive","32 min dive"))+
    
    # Set up plotting theme
    theme_classic()+
    theme(
      legend.position = "none",
      panel.grid.major   = element_line(colour = "grey90",size = 0.2),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=8),
      axis.title.y = element_text(size=8, margin = margin(t=0,r = 0,b = 0,l = 10)),
      strip.background = element_blank(),
      plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))+
    ylab(label[i])
}
plot(1)
i=1
repeat {
  plot(i)
  ggsave(filename = paste("./Figures/a_",para[i],".pdf",sep = ""),width = 3.5,height = 3.5)
  #ggsave(filename = paste("./Figures/",para[i],".jpeg",sep = ""),width = 3.5,height = 3.5)
  i = i+1
  if (i == length(para)+1){
    break
  }
}

para<-c("pco2", "po2","hco3p","pha","pHi","hct","osm","hco3i","sat","lac","cl.rbc","cl.p","hco3i.app","hco3hb")


pdf("./Figures/Figure 1.pdf",width = 180/25.4,height = 4/3*140/25.4)
plot_grid(
  plot(9),plot(2),plot(6), # sat, po2, hct
  plot(4),plot(5),plot(1), #pha, phi, pco2
  plot(3),plot(8),plot(10), # hco3p, hco3, lac
  
  plot(12),plot(11),plot(7), # clp, cli, osm
          nrow = 4,ncol =3,align = 'hv',labels = "AUTO",label_x = .27,label_size = 12)
dev.off()


# Effect of Hb-O2 sat on Hb-HCO3-
  # Is there an effect of saturation on hemoglobin bound hco3-?
  model1<-lmer(hco3hb~sat+(1|animalID),data) # model that includes sat
  model2<-lmer(hco3hb~(1|animalID),data) # model that does not include sat

  anova(model1,model2)
  # model 1 has lowest AIC value, so there is a significant effect of sat on Hb-HCO3-.
  
  model1.sum<-summary(model1)
  model1.sum

 #plot sat vs hb-HCO3
ggplot()+
  #geom_point(mapping = aes(x = sat, y = hco3hb),data = data)+
  
  # Add individual data points and connecting lines
  lapply(
    X = 1:8, 
    FUN = function(k){
      temp.df<-data[data$`animal number`==k,]
      temp.df<-
        data.frame(
          data[data$`animal number`==k,]["sat"],
          data[data$`animal number`==k,]["hco3hb"])
      colnames(temp.df)<-c("sat","hco3hb")
      geom_point(data = temp.df, mapping = aes(x = sat, y = hco3hb), color =cols[k], alpha = 0.5)
    })+ 
  
  lapply(
    X = 1:8, 
    FUN = function(k){
      temp.df<-data[data$`animal number`==k,]
      temp.df<-
        data.frame(
          data[data$`animal number`==k,]["sat"],
          data[data$`animal number`==k,]["hco3hb"])
      colnames(temp.df)<-c("sat","hco3hb")
      geom_line(data = temp.df, mapping = aes(x = sat, y = hco3hb,group = 1), color =cols[k], alpha = 0.5)
    }) +
  
  
  geom_abline(slope = model1.sum$coefficients[2,1],intercept = model1.sum$coefficients[1,1])+
  theme_classic()+
  theme(
    legend.position = "none",
    panel.grid.major   = element_line(colour = "grey90",size = 0.2),
    axis.text = element_text(size=8),
    axis.title = element_text(size=8, margin = margin(t=0,r = 0,b = 0,l = 10)),
    strip.background = element_blank(),
    plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))+
  ylab(expression(bold("[Hb-HCO"[3]^-{}*"] (mmol l"^"-1"*")")))+
  xlab(expression(bold("Hb-O"[2]*" sat (%)")))-> fig3b
fig3b

pdf("./Figures/Figure 3.pdf",width = 90/25.4,height = 120/25.4)
plot_grid(plot(14),fig3b,
          nrow = 2,ncol =1,align = 'hv',labels = "AUTO",label_x = .15,label_size = 12)
dev.off()





## DAVENPORT DIAGRAM ##
df<-
  data.frame(
    HCO3 = c(data$hco3p,data$hco3i),
    pH = c(data$pha,data$pHi),
    Compartment = c(rep("a",length(data$hco3p)),rep("i",length(data$hco3p))),
    Treatment = c(data$treatment,data$treatment)
  )

a<-cbind(aggregate(pH~Treatment+Compartment, df, mean),
         aggregate(pH~Treatment+Compartment, df, stderr)[,3],
         aggregate(HCO3~Treatment+Compartment, df, mean)[,3],
         aggregate(HCO3~Treatment+Compartment, df, stderr)[,3])

colnames(a)<-c("Treatment","Compartment","pH","pH.Err","HCO3","HCO3.Err")
a$Treatment<-c("C1","D1","D2","C1","D1","D2")
a$Compartment<-as.factor(a$Compartment)

# data frame for extracellular non bicarbonate buffering effect from Jensen et al 1998
nbb<-data.frame(
  ph = c(a[1,3],7.4,7.4),
  hco3 = c(a[1,5] ,a[1,5] -22.7*(7.4-a[1,3] ),a[1,5] -20.9*(7.4-a[1,3] ))
)

# data frame for intracellular non bicarbonate buffering effect from Jensen et al 1998
nbb2<-data.frame(
  ph = c(a[4,3],7.1,7.1),
  hco3 = c(a[4,5] ,a[4,5] -107.8*(7.1-a[4,3] ),a[4,5] -78.5*(7.1-a[4,3] ))
)

nbb2<-data.frame(
  ph = c(a[4,3],(20-(a[4,5]-a[4,3]*(-107.8)))/-107.8,(20-(a[4,5]-a[4,3]*(-78.5)))/-78.5),
  hco3 = c(a[4,5],20,20)
)









# Plot and save Figure 2 as PDF
ggplot() + 
  geom_polygon(mapping = aes(nbb$ph, y = nbb$hco3),alpha = 0.5) +
  geom_polygon(mapping = aes(nbb2$ph, y = nbb2$hco3),alpha = 0.5) +
  #scale_x_continuous(limits = c(7.3,7.65)) +
  scale_y_continuous(limits = c(10,22)) +
  scale_shape_manual(values = c(19,17),breaks = c("a","i"),labels = c("Arterial", "Erythrocyte"))+
  geom_point(a, mapping = aes(x = pH, y = HCO3,shape = Compartment),size = 3) +
  geom_path(a, mapping = aes(x = pH, y = HCO3, group = Compartment)) +
  geom_errorbar(width = 0.015,aes(x = a[,3], ymin=a[,5]-a[,6], ymax=a[,5]+a[,6])) +
  geom_errorbarh(height = 0.6, aes(y = a[,5], xmin=a[,3]-a[,4], xmax=a[,3]+a[,4])) +
  theme_classic()+
  theme(legend.position=c(0.15,0.95),
        panel.grid.major   = element_line(colour = "grey90",size = 0.2),
        text = element_text(size=8),
        legend.background = element_rect(fill=alpha('blue', 0)),
        legend.title = element_blank(),
        axis.text = element_text(size=8),
        axis.line = element_line(size = 0.5),
        axis.title = element_text(size=8, margin = margin(t=0,r = 0,b = 0,l = 10)),
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))+
  geom_text(mapping = aes(x = a[,3], y = a[,5], label=c("Pre-dive","18 min dive","32 min dive","Pre-dive","18 min dive","32 min dive")), 
            hjust=-.1, vjust = c(-0.7,-.3,-0.5,-.7,0,-0.8),
            color="black",
            #position = position_dodge(0.9), 
            size=2.8)+
  annotate("text",size=2.8, x = c(7.65,7.62,7.46), y = c(12,22,22), vjust = 0, hjust = 0.5, label = c("10","20","30"),fontface = "bold")+
  stat_function(fun = function(pH) 10*10^(pH-(-0.0817*pH+6.7818))*aco2, linetype = "dotted", colour="black")+
  stat_function(fun = function(pH) 20*10^(pH-(-0.0817*pH+6.7818))*aco2, linetype = "dotted", colour="black")+
  stat_function(fun = function(pH) 30*10^(pH-(-0.0817*pH+6.7818))*aco2, linetype = "dotted", colour="black")+
  labs(x = expression(bold("pH")),
       y = expression(bold("[HCO"[3]^-{}*"] (mmol l"^"-1"*")")))

ggsave(filename = "./Figures/a_Figure 2.pdf",width = 9/2.54,height = 9/2.54)
ggsave(filename = "./Figures/a_Figure 2.jpeg",width = 9/2.54,height = 9/2.54,dpi = 600)
