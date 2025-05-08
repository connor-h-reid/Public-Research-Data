##### Load Packages & Import/Check Data #####

library(zoo)
library(dplyr)
library(plyr)
library(readr)
library(lubridate)

library(ggplot2)
library(ggforce)
library(ggtext)
library(egg)
library(plotly)
library(scales)

library(car)
library(lme4)
library(performance)
library(lmerTest)
library(emmeans)


# General fish data and baseline (angling) blood chemistry files
# Either verify wd or specify exact file paths before file names below (e.g., "C:/Users/.../)

#setwd([...])

master=read.csv("BC Sturgeon General Data.csv")
master$Zone2=as.factor(master$Zone2)
master$Group=as.factor(master$Group)
master$Capture=as.factor("Setline")
str(master)

baseline=read.csv("BC Baseline Blood Data.csv")
baseline$Capture=as.factor("Angled")
str(baseline)


# Accel data for full trials only
fulltrial=master[,c("Tag.ID",
                    "No.Full.Trial",
                    "Zone2",
                    "Length",
                    "Mass",
                    "Group",
                    "SampleS",
                    "Glucose",
                    "Lactate",
                    "TotalS")]
str(fulltrial)

fulltrial=fulltrial[fulltrial$No.Full.Trial==0,]

str(fulltrial)


# This is a good time to get sample sizes, summary lengths/masses
# (for included fish only)
summarystats=as.data.frame(
  fulltrial %>%
    group_by(fulltrial$Group) %>%
    summarise(count=n(),
              meanL=mean(Length),
              sdL=sd(Length),
              meanM=mean(Mass),
              sdM=sd(Mass))
)

summarystats

master %>%
  group_by(master$Group) %>%
  summarise(count=n())

fulltrial %>%
  summarise(count=n(),
            meanL=mean(Length),
            sdL=sd(Length),
            meanM=mean(Mass),
            sdM=sd(Mass))


# Test lengths, masses between groups

lengthmod=glm(log(Length)~Group,data=fulltrial,family="gaussian")
check_model(lengthmod)

par(mfrow=c(2,2))
plot(lengthmod)
par(mfrow=c(1,1))
hist(resid(lengthmod))
shapiro.test(resid(lengthmod))

massmod=glm(log(Mass)~Group,data=fulltrial,family="gaussian")
check_model(massmod)

par(mfrow=c(2,2))
plot(massmod)
par(mfrow=c(1,1))
hist(resid(massmod))
shapiro.test(resid(massmod))

kruskal.test(Length~Group,data=fulltrial)
kruskal.test(Mass~Group,data=fulltrial)


# Make sure manually trimmed accelerometer data files (folder with individual fish files,
# available on GitHub) are together in a folder; set wd as needed

setwd("[...]/CSVs for Analyses")
Myfiles = lapply(list.files(pattern="*.csv"), read.csv)
sturgeon <- data.table::rbindlist(Myfiles, fill=TRUE)
rm(Myfiles)
head(sturgeon)

table(sturgeon$Tag.ID)


##### Calculate Static and Dynamic Accelerations #####

# Make a date/time object
sturgeon$DatetimePosx = as.POSIXct(sturgeon$Timestamp, tz="UTC", format="%m-%d-%Y %H:%M:%OS")
sturgeon$DatetimeMin = strftime(sturgeon$DatetimePosx,tz="UTC",format="%m-%d-%Y %H:%M")
sturgeon$DatetimeMin = as.POSIXct(sturgeon$DatetimeMin,tz="UTC",format="%m-%d-%Y %H:%M")

# Fix temperature and "pressure" (already converted to depth in XManager) names
names(sturgeon)[names(sturgeon) == "Temp....C."] <- "Temp"
names(sturgeon)[names(sturgeon) == "Press...mBar."] <- "Depth"

head(sturgeon)
names(sturgeon)
str(sturgeon)

sturgeon$sum <- abs(sturgeon$X)+abs(sturgeon$Y)+abs(sturgeon$Z)
mean(sturgeon$sum[!is.na(sturgeon$sum)])

sturgeon$X2 <- as.numeric(sturgeon$X)
sturgeon$Y2 <- as.numeric(sturgeon$Y)
sturgeon$Z2 <- as.numeric(sturgeon$Z)


# Find, remove static acceleration from gravity
sturgeon$Xstatic <- rollmean(sturgeon$X2, 50,fill=NA)
sturgeon$Ystatic <- rollmean(sturgeon$Y2, 50,fill=NA)
sturgeon$Zstatic <- rollmean(sturgeon$Z2, 50,fill=NA)

sturgeon$Xstatic[is.na(sturgeon$Xstatic)] <- "cut"
sturgeon2 <- sturgeon[-which(sturgeon$Xstatic=="cut"),]
sturgeon2$Tag.ID=as.factor(sturgeon2$Tag.ID)
sturgeon2$Xstatic=as.numeric(sturgeon2$Xstatic)

# Calculate pitch and roll
sturgeon2$pitch <- as.numeric(sturgeon2$Xstatic)*180/3.14159265359 
sturgeon2$roll <- as.numeric(sturgeon2$Ystatic)*180/3.14159265359 

# Calculate dynamic for each axis
sturgeon2$Xdyn <- sturgeon2$X2-as.numeric(sturgeon2$Xstatic)
sturgeon2$Ydyn <- sturgeon2$Y2-as.numeric(sturgeon2$Ystatic)
sturgeon2$Zdyn <- sturgeon2$Z2-as.numeric(sturgeon2$Zstatic)

str(sturgeon2)



##### Diagnostic Plot Generation #####

### Just verify file paths and remove the "#"s before lines of 
### code in "Depth plots" and "Triaxial Plots" sections

# Depth Plots (set wd wherever works)

#pdf("Sturgeon Depth Plots.pdf", width=16, height=9)

#par(mfrow = c(1,1))
#par(mar = c(10,2,0,2))
#par(oma = c(0,5,5,6))

#for (i in 1:length(sort(unique(sturgeon2$Tag.ID))))
  
#{
  
#  sturgeon2.i <- sturgeon2[sturgeon2$Tag.ID == sort(unique(sturgeon2$Tag.ID))[i],]
  
#  print(ggplot(sturgeon2.i)+ 
#          geom_point(aes(x=DatetimePosx,y=Depth))+
#          scale_y_reverse(name="Depth (m)")+
#          scale_x_datetime(name="Time")+
#          theme_classic(base_size=25)+
#          ggtitle(label = paste("Tag ID = ", sturgeon2.i$Tag[1])) )
  
#}

#dev.off()


# Triaxial Plots
#pdf("C:/Users/CHR/Desktop/Sturgeon Movement Plots.pdf", width=16, height=9)

#par(mfrow = c(1,1))
#par(mar = c(10,2,0,2))
#par(oma = c(0,5,5,6))

#for (i in 1:length(sort(unique(sturgeon2$Tag.ID))))
  
#{
  
#  sturgeon2.i <- sturgeon2[sturgeon2$Tag.ID == sort(unique(sturgeon2$Tag.ID))[i],]
  
#  print(ggplot(sturgeon2.i)+  
#          geom_line(aes(x=DatetimePosx,y=Zdyn , colour= "Tail Beats (Z)"))+
#          geom_line(aes(x=DatetimePosx,y=Ydyn,colour="Y"))+
#          geom_line(aes(x=DatetimePosx,y=Xdyn,colour="Forward (X)"))+
#          theme_bw()+
#          ggtitle(label = paste("Tag ID = ", sturgeon2.i$Tag[1])))
  
#}

#dev.off()



##### Calculate ODBA & VeDBA #####

# Calculate ODBA, VeDBA (look these up, ODBA likely more useful)
sturgeon2$ODBA <- abs(sturgeon2$Xdyn)+abs(sturgeon2$Ydyn)+abs(sturgeon2$Zdyn)
sturgeon2$VeDBA <- sqrt(sturgeon2$Xdyn^2+sturgeon2$Ydyn^2+sturgeon2$Zdyn^2)

head(sturgeon2)

# Assign seconds:
sturgeon2$seconds <- strftime(sturgeon2$DatetimePosx, format="%m-%d-%Y %H:%M:%OS")
sturgeon2seconds <- sturgeon2 %>% group_by(Tag.ID, seconds) %>% summarise(ODBA=sum(ODBA), VeDBA=sum(VeDBA))
head(sturgeon2seconds)

# Assign other variables:
depthtemp <- sturgeon2 %>% filter(!is.na(Depth)) %>% dplyr::select(Tag.ID, seconds, Depth, Temp)
head(depthtemp)

sturgeon2seconds <- merge(sturgeon2seconds, depthtemp, by=c("Tag.ID", "seconds"), all.x=TRUE)
head(sturgeon2seconds)


# Grouping by minute 
i=1
sturgeon2seconds$SecPostRelease <- 1
for (i in 2:length(sturgeon2seconds$Tag.ID)){
  sturgeon2seconds$SecPostRelease[i] <- ifelse(sturgeon2seconds$Tag.ID[i]==sturgeon2seconds$Tag.ID[i-1],
                                               sturgeon2seconds$SecPostRelease[i-1]+1, 1)
}


# For 15 s instead of 30s, adjust #s (e.g., 16 instead of 31, 0.25 instead of 0.5, etc.)?
sturgeon2seconds$MinPostRelease <- 
                                ifelse(sturgeon2seconds$SecPostRelease<31, 0.5,
                                ifelse(sturgeon2seconds$SecPostRelease<61, 1,
                                ifelse(sturgeon2seconds$SecPostRelease<91, 1.5,
                                ifelse(sturgeon2seconds$SecPostRelease<121, 2,
                                ifelse(sturgeon2seconds$SecPostRelease<151, 2.5,
                                ifelse(sturgeon2seconds$SecPostRelease<181, 3,
                                ifelse(sturgeon2seconds$SecPostRelease<211, 3.5,
                                ifelse(sturgeon2seconds$SecPostRelease<241, 4,
                                ifelse(sturgeon2seconds$SecPostRelease<271, 4.5,
                                ifelse(sturgeon2seconds$SecPostRelease<301, 5,
                                ifelse(sturgeon2seconds$SecPostRelease<331, 5.5,
                                ifelse(sturgeon2seconds$SecPostRelease<361, 6,
                                ifelse(sturgeon2seconds$SecPostRelease<391, 6.5,
                                ifelse(sturgeon2seconds$SecPostRelease<421, 7,
                                ifelse(sturgeon2seconds$SecPostRelease<451, 7.5,
                                ifelse(sturgeon2seconds$SecPostRelease<481, 8,
                                ifelse(sturgeon2seconds$SecPostRelease<511, 8.5,
                                ifelse(sturgeon2seconds$SecPostRelease<541, 9,
                                ifelse(sturgeon2seconds$SecPostRelease<571, 9.5,
                                ifelse(sturgeon2seconds$SecPostRelease<602, 10,10.5))))))))))))))))))))


#summarise
sturgeonSum <- sturgeon2seconds %>% dplyr::group_by(Tag.ID,MinPostRelease) %>% dplyr::summarise(ODBA=mean(ODBA), VeDBA=mean(VeDBA), samples=length(Tag.ID))
head(sturgeonSum)

#add back in depth and temp:
head(sturgeon2seconds)
depthtemp2 <- sturgeon2seconds %>% filter(!is.na(Depth)) %>% dplyr::select(Tag.ID, MinPostRelease, Depth, Temp, ODBA, VeDBA) 
head(depthtemp2)  
str(depthtemp2)

#summarise mean depth, sd depth, and temperature by minute:
depthtempsum <- depthtemp2 %>% dplyr::group_by(Tag.ID, MinPostRelease) %>% dplyr::summarise(mean_depth=mean(Depth), sd_depth=sd(Depth),temp=mean(Temp))
head(depthtempsum)

sturgeonSum <- merge(sturgeonSum, depthtempsum, by=c("Tag.ID", "MinPostRelease"), all.x=TRUE)
sturgeonSum <- sturgeonSum %>% arrange(Tag.ID, MinPostRelease)


#Plot

head(sturgeonSum)

### Export to CSV if desirable, add file location
#write.csv(sturgeonSum,"ODBA Data.csv")



##### Sampling Logistics & Blood Chemistry #####

# Re-update wd if desired
# setwd("[...]")


# Sampling times

mean(master$BloodS)
sd(master$BloodS)
range(master$BloodS)

mean(master$UltraS,na.rm=T)
sd(master$UltraS,na.rm=T)
range(master$UltraS,na.rm=T)

mean(master$BiopsyS,na.rm=T)
sd(master$BiopsyS,na.rm=T)
range(master$BiopsyS,na.rm=T)


# If you *really* need a model to show that doing more things adds time
timemod=glm(SampleS~Group,data=master,family="gaussian")
par(mfrow=c(2,2))
plot(timemod)
par(mfrow=c(1,1))
hist(resid(timemod))
shapiro.test(resid(timemod))

Anova(timemod,test.statistic="F")
emmeans(timemod,pairwise~Group,type="response")

# Fit's not perfect but even log transform doesn't change the actual interpretation
timemod2=glm(log(SampleS)~Group,data=master,family="gaussian")
par(mfrow=c(2,2))
plot(timemod2)
par(mfrow=c(1,1))
hist(resid(timemod2))
shapiro.test(resid(timemod2))

Anova(timemod2,test.statistic="F")
emmeans(timemod2,pairwise~Group,type="response")

sampletimeplot=ggplot(master,aes(x=Group,y=SampleS))+
  geom_violin(scale="area",fill="grey70")+
  geom_sina(scale="area",size=1.1)+
  xlab("")+ylab("Required Sampling Time (s)\n")+
  scale_y_continuous(expand=c(0,0),limits=c(120,1200))+
  theme_classic()+
  guides(fill="none")+
  theme(axis.text.x=element_text(colour="black",size=18),
        axis.text.y=element_text(colour="black",size=18),
        axis.title.y=element_text(size=18),
        panel.background=element_blank())+
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))+
  annotate("text",label="italic(P) < 0.001",x=1.5,y=1000,parse=T,size=6)
sampletimeplot

set.seed(0)
ggsave("Sampling Time Plot.png",sampletimeplot,width=13,height=13,units="cm")


# Time spent on retrieved line before processing

plot(master$LineS~master$Group)

mean(master$LineS,na.rm=T)
sd(master$LineS,na.rm=T)
range(master$LineS,na.rm=T)


# Total time from setline retrieval to release

plot(master$TotalS~master$Group)

mean(master$TotalS)
sd(master$TotalS)
range(master$TotalS)


# Blood chemistry

range(baseline$BloodS)
range(master$BloodS)

blooddata=merge(master[,c("Tag.ID","Capture","Length","Mass","BloodS","Glucose","Lactate")],baseline[,-2],all=T)
str(blooddata)

# Can't analyse lactate as only one angled/baseline was above lower detection limit on device 
# (see "LO" values in almost all angled fish):
View(blooddata)

str(blooddata)
blooddata$Lactate[blooddata$Lactate=="LO"]=NA
blooddata$Lactate=as.numeric(blooddata$Lactate)
str(blooddata)

# Glucose and lactate plots

glumod=glm(log(Glucose)~Capture+BloodS,data=blooddata,family="gaussian")
par(mfrow=c(2,2))
plot(glumod)
par(mfrow=c(1,1))
hist(resid(glumod))
shapiro.test(resid(glumod))

Anova(glumod,type="II",test.statistic="F")

emmeans(glumod,pairwise~Capture,type="response")

gluplot=ggplot(blooddata,aes(x=Capture,y=Glucose))+
  geom_violin(scale="count",fill="grey70")+
  geom_sina(aes(fill=Capture,group=Capture),scale="count",size=1.1)+
  xlab("")+ylab("Blood Glucose (mmol/l)\n")+
  scale_y_continuous(expand=c(0,0),limits=c(0,10))+
  theme_classic()+
  guides(fill="none")+
  theme(axis.text.x=element_text(colour="black",size=18),
        axis.text.y=element_text(colour="black",size=18),
        axis.title.y=element_text(size=18),
        panel.background=element_blank())+
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))+
  annotate("text",label="italic(P) < 0.001",x=1.75,y=7.5,parse=T,size=6)+
  annotate("pointrange",x=1,y=3.98,ymin=3.60,ymax=4.41,size=1.5,linewidth=2.5,colour="black",alpha=0.5)+
  annotate("pointrange",x=2,y=1.45,ymin=1.26,ymax=1.67,size=1.5,linewidth=2.5,colour="black",alpha=0.5)
gluplot

lacplot=ggplot(blooddata,aes(x=Capture,y=Lactate))+
  geom_violin(scale="count",fill="grey70")+
  geom_sina(aes(fill=Capture,group=Capture),scale="count",size=1.1)+
  geom_rect(aes(ymin=0,ymax=0.3,xmin=0,xmax=3),fill="grey90")+
  geom_hline(aes(yintercept=0.3),lty=2)+
  geom_vline(aes(xintercept=0))+
  xlab("")+ylab("Blood Lactate (mmol/l)\n")+
  scale_y_continuous(expand=c(0,0),limits=c(0,15))+
  theme_classic()+
  guides(fill="none")+
  theme(axis.text.x=element_text(colour="black",size=18),
        axis.text.y=element_text(colour="black",size=18),
        axis.title.y=element_text(size=18),
        panel.background=element_blank())+
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
lacplot

set.seed(0)
bloodplot=ggarrange(gluplot,lacplot,ncol=2,nrow=1,labels=c("A","B"),
                    label.args=list(gp=grid::gpar(font=2,cex=1.7)))
ggsave("bloodplot.png",bloodplot,width=25,height=15,units="cm")


master %>%
  summarise(count=n(),
            mean(Glucose),
            mean(Lactate))

mean(blooddata[blooddata$Capture=="Setline",]$Lactate,na.rm=T)
range(blooddata[blooddata$Capture=="Setline",]$Lactate,na.rm=T)


plot(master$Glucose~master$LineS)
plot(master$Lactate~master$LineS)

glumod2=glm(log(Glucose)~LineS,data=master,family="gaussian")
par(mfrow=c(2,2))
plot(glumod2)
par(mfrow=c(1,1))
hist(resid(glumod2))
shapiro.test(resid(glumod2))

Anova(glumod2,type="II",test.statistic="F")

lacmod=glm(log(Lactate)~LineS,data=master,family="gaussian")
par(mfrow=c(2,2))
plot(lacmod)
par(mfrow=c(1,1))
hist(resid(lacmod))
shapiro.test(resid(lacmod))

Anova(lacmod,type="II",test.statistic="F")

##### ODBA & VeDBA Analyses #####

str(sturgeonSum)

acceldat=merge(sturgeonSum,fulltrial,by="Tag.ID")
str(acceldat)

# Confirm no fish made it through without full trial
unique(acceldat$Tag.ID[acceldat$No.Full.Trial==1])

# ODBA & VeDBA in first 30 s
releasedat=acceldat[acceldat$MinPostRelease==0.5,]

releaseODBA=glm(log(ODBA)~Group+temp+Length+SampleS,data=releasedat)

par(mfrow=c(2,2))
plot(releaseODBA)
par(mfrow=c(1,1))
hist(resid(releaseODBA))
shapiro.test(resid(releaseODBA))

Anova(releaseODBA,type="II",test.statistic="F")
plot(releasedat$ODBA~releasedat$Group)
plot(releasedat$ODBA~releasedat$temp)

releaseVeDBA=glm(log(VeDBA)~Group+temp+Length+SampleS,data=releasedat)

par(mfrow=c(2,2))
plot(releaseVeDBA)
par(mfrow=c(1,1))
hist(resid(releaseVeDBA))
shapiro.test(resid(releaseVeDBA))

Anova(releaseVeDBA,type="II",test.statistic="F")
plot(releasedat$VeDBA~releasedat$Group)
plot(releasedat$VeDBA~releasedat$temp)


# ODBA & VeDBA with time as predictor

# ODBA
# Check inclusion of individual fish as random effect
ODBAmod=lmer(log(ODBA)~Group+temp+Length+SampleS+MinPostRelease+(1|Tag.ID),data=acceldat)
ranova(ODBAmod)

check_model(ODBAmod)

# Fixed effects analysis
Anova(ODBAmod,type="II",test.statistic="F")
summary(ODBAmod)
ODBArandvariance=0.05343/(0.05343+0.16907)
ODBArandvariance*100

plot(acceldat$ODBA~acceldat$temp)

# VeDBA
# Check inclusion of individual fish as random effect
VeDBAmod=lmer(log(VeDBA)~Group+temp+Length+SampleS+MinPostRelease+(1|Tag.ID),data=acceldat)
ranova(VeDBAmod)

check_model(VeDBAmod)

# Fixed effects analysis
Anova(VeDBAmod,type="II",test.statistic="F")
summary(VeDBAmod)
VeDBArandvariance=0.04498/(0.04498+0.16188)
VeDBArandvariance*100

# Plot
ODBAplot=ggplot(acceldat,aes(x=MinPostRelease,y=ODBA,alpha=0.5))+
  geom_point()+facet_wrap("Group")+
  geom_line(aes(x=MinPostRelease,y=ODBA,group=Tag.ID))+
  scale_y_continuous(expand=c(0,0),limits=c(0,20.5))+
  ylab("ODBA (*g*)")+xlab("")+
  theme_classic()+
  theme(legend.position="none",
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
  theme(axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title.x=element_text(size=14),
        axis.title.y=ggtext::element_markdown(size=14),
        strip.text=element_text(size=14))
ODBAplot


VeDBAplot=ggplot(acceldat,aes(x=MinPostRelease,y=VeDBA,alpha=0.5))+
  geom_point()+facet_wrap("Group")+
  geom_line(aes(x=MinPostRelease,y=VeDBA,group=Tag.ID))+
  scale_y_continuous(expand=c(0,0),limits=c(0,15))+
  ylab("VeDBA (*g*)\n")+xlab("\nMinutes Post-Release")+
  theme_classic()+
  theme(legend.position="none",
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
  theme(axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title.x=element_text(size=14),
        axis.title.y=ggtext::element_markdown(size=14),
        strip.text=element_text(size=14))
VeDBAplot

# SAVE THIS PLOT

OVeDBAplot=ggarrange(ODBAplot,VeDBAplot,ncol=1,nrow=2,labels=c("A","B"),
                     label.args=list(gp=grid::gpar(font=2,cex=1.7)))
ggsave("OVeDBAplot.png",OVeDBAplot,height=20,width=30,units="cm")


# Depth & temperature (not really interesting or important I think)

depthmod=lmer(mean_depth~Group+Length+SampleS+MinPostRelease+(1|Tag.ID),data=acceldat,REML=F)
check_model(depthmod)

depthmod2=lm(mean_depth~Group+Length+SampleS+MinPostRelease,data=acceldat)
anova(depthmod,depthmod2)

depthmod=lmer(mean_depth~Group+Length+SampleS+MinPostRelease+(1|Tag.ID),data=acceldat)

Anova(depthmod,type="II",test.statistic="F")
summary(depthmod)

Depthrandvariance=13.068/(13.068+2.701)
Depthrandvariance*100


depthplot=ggplot(acceldat,aes(x=MinPostRelease,y=mean_depth,alpha=0.5))+
  geom_point()+facet_wrap("Group")+
  geom_line(aes(x=MinPostRelease,y=mean_depth,group=Tag.ID))+
  scale_y_continuous(expand=c(0,0),limits=c(20,0),trans="reverse")+
  ylab("Mean Depth (m)\n")+xlab("\nMinutes Post-Release")+
  theme_classic()+
  theme(legend.position="none",
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
  theme(axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        strip.text=element_text(size=14))
depthplot

ggsave("Sturgeon Depth Plot.png",depthplot,height=10,width=20,units="cm")
