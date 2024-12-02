
##### Notes #####

# Part II:
# Discard fish 2 (MS-222) because dose was too low (~stage III only)
# Discard fish 13 (clove oil) from location/orientation behaviour analyses (too many missing values)
# Discard blood chemistry data for fish BP1, BP2, BP3, and 2-9 because of sampling technique errors

##### Data/Housekeeping #####

library(car)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(AER)
library(factoextra)
library(plotly)
library(dplyr)
library(nnet)

# Data for Part 1 (Immobilisation Pilot)

pilot=read.csv("C:/Users/CHR/Documents/FECPL/2024/Manitoba/Data/MB Part 1 Data.csv")
pilot$Treatment=as.factor(pilot$Treatment)
levels(pilot$Treatment)
levels(pilot$Treatment)=c("Control","cDC","TENS")
str(pilot)

mean(pilot$FL)
sd(pilot$FL)
range(pilot$FL)

mean(pilot$Mass)
sd(pilot$Mass)
range(pilot$Mass)

# Part 1 - Fork lengths and masses vs. treatments

plot(pilot$FL~pilot$Treatment)

lengthmod=glm(FL~Treatment,data=pilot,family="gaussian")
par(mfrow=c(2,2))
plot(lengthmod)
par(mfrow=c(1,1))

Anova(lengthmod,type="II",test.statistic="F")

massmod=glm(Mass~Treatment,data=pilot,family="gaussian")
par(mfrow=c(2,2))
plot(massmod)
par(mfrow=c(1,1))

Anova(massmod,type="II",test.statistic="F")


# Data for Part 2 (Main Part)

surgerydat=read.csv("C:/Users/CHR/Documents/FECPL/2024/Manitoba/Data/MB Part 2 Surgery Data.csv")
surgerydat$Treatment=as.factor(surgerydat$Treatment)
surgerydat=surgerydat[surgerydat$FishID!="2",]
str(surgerydat)

behavdat1=read.csv("C:/Users/CHR/Documents/FECPL/2024/Manitoba/Data/MB Part 2 Single Behaviour Scores.csv")
behavdat1$Fish=as.factor(behavdat1$Fish)
behavdat1$Treatment=as.factor(behavdat1$Treatment)
behavdat1$Tank=as.factor(behavdat1$Tank)
behavdat1$Time=as.factor(behavdat1$Time)
behavdat1$Rheo1=as.factor(behavdat1$Rheo1)
behavdat1=behavdat1[behavdat1$Fish!="2",]
levels(behavdat1$Treatment)=c("Baseline","Clove Oil","Handling","MS-222","TENS")
str(behavdat1)

treatorder <- c("Baseline","Handling","TENS","MS-222","Clove Oil") 


# Part 2 - Fork lengths and masses vs. treatments

part2sizedat=read.csv("C:/Users/CHR/Documents/FECPL/2024/Manitoba/Data/MB Part 2 Size Data.csv")
str(part2sizedat)

mean(part2sizedat$FL,na.rm=T)
sd(part2sizedat$FL,na.rm=T)
range(part2sizedat$FL,na.rm=T)

mean(part2sizedat$Mass,na.rm=T)
sd(part2sizedat$Mass,na.rm=T)
range(part2sizedat$Mass,na.rm=T)

lengthmod2=glm(FL~Treatment,data=part2sizedat,family="gaussian")
par(mfrow=c(2,2))
plot(lengthmod2)
par(mfrow=c(1,1))

Anova(lengthmod2,type="II",test.statistic="F")

massmod2=glm(Mass~Treatment,data=part2sizedat,family="gaussian")
par(mfrow=c(2,2))
plot(massmod2)
par(mfrow=c(1,1))

Anova(massmod2,type="II",test.statistic="F")


# Part 2 - Treatment/tank distributions

set.seed(2024)
chisq.test(behavdat1$Treatment,behavdat1$Tank,simulate.p.value=T,B=5000)

table(behavdat1$Treatment,behavdat1$Tank)


##### Part I: Pilot #####

pilot$K=100*(pilot$Mass/(pilot$FL/10)^3)

pilot$TotalCounts=pilot$Count1+pilot$Count2
pilot$TotalSeqs=pilot$Sequences1+pilot$Sequences2

countmod=glm(TotalCounts~Treatment+FL+K,data=pilot,family="poisson")
dispersiontest(countmod)
countmod=glm(TotalCounts~Treatment+FL+K,data=pilot,family="quasipoisson")

Anova(countmod,type="II",test.statistic="LR")
plot(pilot$TotalCounts~pilot$Treatment)

seqmod=glm(TotalSeqs~Treatment+FL+K,data=pilot,family="poisson")
dispersiontest(seqmod)

Anova(seqmod,type="II",test.statistic="LR")
plot(pilot$TotalSeqs~pilot$Treatment)

emmeans(seqmod,pairwise~Treatment,type="response")

mean(pilot$mAPSU,na.rm=T)
range(pilot$mAPSU,na.rm=T)

mean(pilot$TENS,na.rm=T)
range(pilot$TENS,na.rm=T)

table(pilot$TotalCounts,pilot$Treatment)
table(pilot$TotalSeqs,pilot$Treatment)

moveplot1=ggplot(pilot,aes(x=Treatment,y=TotalCounts))+
  geom_dotplot(binaxis="y",stackdir="center",dotsize=0.5)+
  theme_classic()+
  ylab("Total Individual Movement Counts\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(-0.25,15))+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12))+
  annotate(geom="text",label=c("n=16","n=13","n=16"),x=c(1,2,3),y=14.25,size=5)
moveplot1

moveplot2=ggplot(pilot,aes(x=Treatment,y=TotalSeqs))+
  geom_dotplot(binaxis="y",stackdir="center",dotsize=0.5)+
  theme_classic()+
  ylab("Total Movement Event Counts\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(-0.1,5.5))+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12))+
  annotate(geom="text",label=c("n=15","n=11","n=14"),x=c(1,2,3),y=5.25,size=5)
moveplot2

totalmoveplot=ggarrange(moveplot1,moveplot2,labels="AUTO",widths=c(1,1))
totalmoveplot
ggsave("Individual Movements and Sequences Plot.png",totalmoveplot,dpi=300,width=18,height=10,units="cm")


##### Part II: Induction and Recovery Times, Surgery Movements #####

#  Induction times:

inductmod=glm(Tinduction~Treatment+Mass,data=surgerydat,family="gaussian")
par(mfrow=c(2,2))
plot(inductmod)
par(mfrow=c(1,1))
shapiro.test(resid(inductmod))

Anova(inductmod,type="II",test.statistic="F")
emmeans(inductmod,pairwise~Treatment)


# Voluntary movements:

movemod=glm(Movements~Treatment,data=surgerydat[surgerydat$Treatment!="MS-222",],family="poisson")
dispersiontest(movemod)

Anova(movemod,type="II",test.statistic="LR")
emmeans(movemod,pairwise~Treatment,type="response")


# Surgery times:

timemod=glm(Tsurgery~Treatment,data=surgerydat,family="gaussian")
par(mfrow=c(2,2))
plot(timemod)
par(mfrow=c(1,1))

Anova(timemod,type="II",test.statistic="F")
emmeans(timemod,pairwise~Treatment)


# RAMP scores:

surgerydat$RAMP=surgerydat$Tail+surgerydat$Flex+surgerydat$Equil+surgerydat$Eye

RAMPmod=glm(RAMP~Treatment,data=surgerydat,family="poisson")
dispersiontest(RAMPmod)

Anova(RAMPmod,type="II",test.statistic="LR")
emmeans(RAMPmod,pairwise~Treatment,type="response")

table(surgerydat$RAMP,surgerydat$Treatment)


# Recovery time:

recovmod=glm(TotalRecov~Treatment,data=surgerydat,family="gaussian")
par(mfrow=c(2,2))
plot(recovmod)
par(mfrow=c(1,1))

Anova(recovmod,type="II",test.statistic="F")
emmeans(recovmod,pairwise~Treatment)


##### Part II: Single Behaviour Scores #####

# Swimming behaviour: moving upstream, downstream, or slacking off

clusterdat=behavdat1[complete.cases(behavdat1[,c(12,16,17)]),]
clusterdat$Downstreams=abs(clusterdat$Downstreams)

clusters=clusterdat[,c(12,16,17)]
means=apply(clusters,2,mean)
stdevs=apply(clusters,2,sd)
normals=scale(clusters,center=means,scale=stdevs)

normals=as.data.frame(normals)
names(normals)[names(normals)=="SwimTime"]="NorSwimTime"
names(normals)[names(normals)=="Upstreams"]="NorUpstreams"
names(normals)[names(normals)=="Downstreams"]="NorDownstreams"
str(normals)

clusterdat=cbind(clusterdat,normals)
str(clusterdat)

fviz_nbclust(normals, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")
fviz_nbclust(normals, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

clusterdat$clusters=factor(kmeans(clusterdat[,c(26:28)],3)$cluster)
str(clusterdat)

levels(clusterdat$clusters)
levels(clusterdat$clusters)=c("Downstream","Inert","Upstream")

cluster.plot=plot_ly(clusterdat,x=~NorUpstreams,y=~NorDownstreams, 
                  z=~NorSwimTime,color=~clusters,colors=c("grey80","grey45","black"),
                  marker=list(size=6,
                              line=list(color="black",
                                          width = 1))) %>%
  add_markers() %>%
  layout(showlegend=TRUE,legend=list(itemsizing="trace"))
cluster.plot


plot(clusterdat$SwimTime~clusterdat$clusters)
c1mod=glm(SwimTime~clusters,data=clusterdat,family="gaussian")
par(mfrow=c(2,2))
plot(c1mod)
par(mfrow=c(1,1))
Anova(c1mod,type="II",test.statistic="F")
emmeans(c1mod,pairwise~clusters)

plot(clusterdat$Upstreams~clusterdat$clusters)
c2mod=glm(Upstreams~clusters,data=clusterdat,family="poisson")
dispersiontest(c2mod)
c2mod=glm(Upstreams~clusters,data=clusterdat,family="quasipoisson")
Anova(c2mod,type="II",test.statistic="F")
emmeans(c2mod,pairwise~clusters,type="response")

plot(clusterdat$Downstreams~clusterdat$clusters)
c3mod=glm(Downstreams~clusters,data=clusterdat,family="poisson")
dispersiontest(c3mod)
c3mod=glm(Downstreams~clusters,data=clusterdat,family="quasipoisson")
Anova(c3mod,type="II",test.statistic="F")
emmeans(c3mod,pairwise~clusters,type="response")


multmod=multinom(clusters~Treatment,data=clusterdat)
clusterprobs=as.data.frame(emmeans(multmod,~clusters|Treatment,type="prob"))
clusterprobs
plot(clusterprobs)
str(clusterprobs)

Anova(multmod)

trtmeans=emmeans(multmod,pairwise~clusters|Treatment,mode="prob")
trtmeans
trtmeansdat=as.data.frame(trtmeans$emmeans)
trtmeansdat$labelY=trtmeansdat$upper.CL+0.1
trtmeansdat$label=c("a","b","a",
                    "b","a","b",
                    "a","a","b",
                    "b","a","b",
                    "a","a","a")
trtmeansdat

clustermeans=emmeans(multmod,pairwise~Treatment|clusters,mode="prob")
clustermeans
clustermeansdat=as.data.frame(clustermeans$emmeans)
clustermeansdat$labelY=clustermeansdat$upper.CL+0.06
clustermeansdat$label=c("ab","b","a","ab","ab",
                        "c","a","b","ab","b",
                        "a","b","b","b","ab")

clustermeansdat

trtclusters=ggplot(clusterprobs,aes(x=clusters,y=prob))+
  geom_point(cex=3)+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0.3,lwd=1)+
  facet_wrap(~factor(Treatment,levels=c("Baseline","Handling","TENS","MS-222","Clove Oil")))+
  theme_bw()+scale_y_continuous(expand=c(0,0),limits=c(-0.03,1.15),breaks=c(0,0.25,0.5,0.75,1.0))+
  theme(strip.text.x=element_text(size=12),
        axis.text=element_text(size=12,colour="black"),
        axis.title.y=element_text(size=12))+
  theme(panel.grid=element_blank())+
  ylab("LS Mean Proportion of Fish\n")+
  xlab("")+
  geom_text(data=trtmeansdat,aes(x=clusters,y=labelY,label=label),size=5)
trtclusters
ggsave("Proportion of Fish in Clusters Within Treatments.png",trtclusters,width=20,height=13,units="cm",dpi=300)

clustertrts=ggplot(clusterprobs,aes(x=factor(Treatment,levels=c("Baseline","Handling","TENS","MS-222","Clove Oil")),y=prob))+
  geom_point(cex=3)+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0.3,lwd=1)+
  facet_wrap(~clusters)+
  theme_bw()+scale_y_continuous(expand=c(0,0),limits=c(-0.03,1.11))+
  theme(strip.text.x=element_text(size=12),
        axis.text=element_text(size=12,colour="black"),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(angle=35,hjust=1,vjust=1))+
  theme(panel.grid=element_blank())+
  ylab("LS Mean Proportion of Fish\n")+
  xlab("")+
  geom_text(data=clustermeansdat,aes(x=Treatment,y=labelY,label=label),size=5)
clustertrts
ggsave("Proportion of Fish in Clusters Across Treatments.png",clustertrts,width=20,height=10,units="cm",dpi=300)

totalclusterplot=ggarrange(trtclusters,clustertrts,labels="AUTO",widths=c(1,1),heights=c(1.1,1),ncol=1,nrow=2)
totalclusterplot
ggsave("Total Clusters Plot.png",totalclusterplot,width=20,height=25,units="cm")

length(clusterdat$NorSwimTime[clusterdat$Treatment=="Baseline" & !is.na(clusterdat$NorSwimTime)])
length(clusterdat$NorSwimTime[clusterdat$Treatment=="Handling" & !is.na(clusterdat$NorSwimTime)])
length(clusterdat$NorSwimTime[clusterdat$Treatment=="TENS" & !is.na(clusterdat$NorSwimTime)])
length(clusterdat$NorSwimTime[clusterdat$Treatment=="MS-222" & !is.na(clusterdat$NorSwimTime)])
length(clusterdat$NorSwimTime[clusterdat$Treatment=="Clove Oil" & !is.na(clusterdat$NorSwimTime)])



# Swimming time

behavdat1$K=100*(behavdat1$Mass/(behavdat1$FL/10)^3)

swimmod=glm(SwimTime~Treatment+K,data=behavdat1,family="gaussian")
par(mfrow=c(2,2))
plot(swimmod)
par(mfrow=c(1,1))

Anova(swimmod,type="II",test.statistic="F")
summary(swimmod)
swimmeans=emmeans(swimmod,pairwise~Treatment)
swimmeans

swimemm=as.data.frame(swimmeans$emmeans)
levels(swimemm$Treatment)

ggplot(swimemm,aes(x=factor(Treatment,level=treatorder),y=emmean))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Time Spent Swimming (s)\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,320))+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title=element_text(size=12))+
  annotate("text",x=1,y=314.90,label="a",size=6)+
  annotate("text",x=2,y=183.84,label="b",size=6)+
  annotate("text",x=3,y=191.43,label="b",size=6)+
  annotate("text",x=4,y=141.51,label="bc",size=6)+
  annotate("text",x=5,y=104.02,label="c",size=6)

Tapply(behavdat1$SwimTime~behavdat1$Treatment,fun="max",na.rm=T)

swimplot=ggplot(behavdat1,aes(x=factor(Treatment,level=treatorder),y=SwimTime))+
  geom_violin(fill="grey80",scale="count")+
  geom_point()+
  geom_point(data=swimemm,aes(x=factor(Treatment,level=treatorder),y=emmean),size=5,shape=1)+
  geom_errorbar(data=swimemm,aes(x=Treatment,y=emmean,ymin=lower.CL,ymax=upper.CL),
                width=0,linewidth=3,alpha=0.3)+
  ylab("Time Spent Swimming (s)\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,320))+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title=element_text(size=12))+
  annotate("text",x=1,y=312,label="a",size=6)+
  annotate("text",x=2,y=299,label="b",size=6)+
  annotate("text",x=3,y=312,label="b",size=6)+
  annotate("text",x=4,y=233,label="bc",size=6)+
  annotate("text",x=5,y=235,label="c",size=6)
swimplot
ggsave("Swimming Time.png",swimplot,dpi=300,height=12,width=12,units="cm")


# Swimming initiations

startmod=glm(SwimStarts~Treatment+K,data=behavdat1,family="poisson")
dispersiontest(startmod)
startmod=glm(SwimStarts~Treatment+K,data=behavdat1,family="quasipoisson")

Anova(startmod,type="II",test.statistic="LR")
startmeans=emmeans(startmod,pairwise~Treatment,type="response")
startmeans

startemm=as.data.frame(startmeans$emmeans)

ggplot(startemm,aes(x=factor(Treatment,level=treatorder),y=rate))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Swim Starts\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,12))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))


# 180º turns

turnmod=glm(Turns~Treatment+K,data=behavdat1,family="poisson")
dispersiontest(turnmod)
turnmod=glm(Turns~Treatment+K,data=behavdat1,family="quasipoisson")

Anova(turnmod,type="II",test.statistic="LR")
turnmeans=emmeans(turnmod,pairwise~Treatment,type="response")
turnmeans

turnemm=as.data.frame(turnmeans$emmeans)

turnplot=ggplot(turnemm,aes(x=factor(Treatment,level=treatorder),y=rate))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Turns\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,20))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=17.93,label="a",size=6)+
  annotate("text",x=2,y=10.93,label="a",size=6)+
  annotate("text",x=3,y=12.54,label="a",size=6)+
  annotate("text",x=4,y=3.10,label="b",size=6)+
  annotate("text",x=5,y=2.53,label="b",size=6)
turnplot
ggsave("Number of Turns.png",turnplot,dpi=300,height=12,width=12,units="cm")

behavplot=ggarrange(swimplot,turnplot,labels="AUTO",widths=c(1,1))
behavplot
ggsave("Swimming and Turning Plot.png",behavplot,width=25,height=13,units="cm")

length(behavdat1$SwimTime[behavdat1$Treatment=="Baseline" & !is.na(behavdat1$SwimTime)])
length(behavdat1$SwimTime[behavdat1$Treatment=="TENS" & !is.na(behavdat1$SwimTime)])
length(behavdat1$SwimTime[behavdat1$Treatment=="Handling" & !is.na(behavdat1$SwimTime)])
length(behavdat1$SwimTime[behavdat1$Treatment=="MS-222" & !is.na(behavdat1$SwimTime)])
length(behavdat1$SwimTime[behavdat1$Treatment=="Clove Oil" & !is.na(behavdat1$SwimTime)])

length(behavdat1$Turns[behavdat1$Treatment=="Baseline" & !is.na(behavdat1$Turns)])
length(behavdat1$Turns[behavdat1$Treatment=="TENS" & !is.na(behavdat1$Turns)])
length(behavdat1$Turns[behavdat1$Treatment=="Handling" & !is.na(behavdat1$Turns)])
length(behavdat1$Turns[behavdat1$Treatment=="MS-222" & !is.na(behavdat1$Turns)])
length(behavdat1$Turns[behavdat1$Treatment=="Clove Oil" & !is.na(behavdat1$Turns)])



# Initial rheotaxis

table(behavdat1$Rheo1,behavdat1$Treatment)

rheodat=behavdat1[behavdat1$Rheo1!="Neutral",]
table(rheodat$Rheo1,rheodat$Treatment)

chisq.test(behavdat1$Treatment,behavdat1$Rheo1,simulate.p.value=T,B=10000)

rheomod=multinom(behavdat1$Rheo1~behavdat1$FL)
Anova(rheomod)


# Edge, outer, inner proportions

behavdat1$EdgeProp=behavdat1$EdgeCount/(behavdat1$EdgeCount+behavdat1$OuterCount+behavdat1$InnerCount)
behavdat1$OuterProp=behavdat1$OuterCount/(behavdat1$EdgeCount+behavdat1$OuterCount+behavdat1$InnerCount)
behavdat1$InnerProp=behavdat1$InnerCount/(behavdat1$EdgeCount+behavdat1$OuterCount+behavdat1$InnerCount)
str(behavdat1)

edgemod=glm(EdgeProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(edgemod)
edgemeans=emmeans(edgemod,pairwise~Treatment,type="response",adjust="bonferroni")
edgemeans

edgeemm=as.data.frame(edgemeans$emmeans)
edgeemm

edgeplot=ggplot(edgeemm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Edge Zone\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=1.050,label="a",size=6)+
  annotate("text",x=2,y=0.899,label="ab",size=6)+
  annotate("text",x=3,y=0.990,label="a",size=6)+
  annotate("text",x=4,y=0.769,label="bc",size=6)+
  annotate("text",x=5,y=0.678,label="c",size=6)
edgeplot

outermod=glm(OuterProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(outermod)

outermeans=emmeans(outermod,pairwise~Treatment,type="response")
outermeans

outeremm=as.data.frame(outermeans$emmeans)
outeremm

outerplot=ggplot(outeremm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Outer Zone\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=0.132,label="a",size=6)+
  annotate("text",x=2,y=0.209,label="ab",size=6)+
  annotate("text",x=3,y=0.200,label="ab",size=6)+
  annotate("text",x=4,y=0.320,label="bc",size=6)+
  annotate("text",x=5,y=0.363,label="c",size=6)
outerplot

innermod=glm(InnerProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(innermod)
innermeans=emmeans(innermod,pairwise~Treatment,type="response")
innermeans

inneremm=as.data.frame(innermeans$emmeans)
inneremm

#Only 1 inner count in 1 BB fish, so remove CIs for this fish
#As conservative estimate, assign "significance" vs. treatments with CIs not going below 0.01667?
  # (that'd be all but TENS)

inneremm$asymp.UCL=replace(inneremm$asymp.UCL,inneremm$asymp.UCL>0.75,0.005)

innerplot=ggplot(inneremm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Inner Zone\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=0.060,label="",size=6)+
  annotate("text",x=2,y=0.289,label="ab",size=6)+
  annotate("text",x=3,y=0.153,label="a",size=6)+
  annotate("text",x=4,y=0.354,label="ab",size=6)+
  annotate("text",x=5,y=0.418,label="b",size=6)
innerplot

edgeemm$zone="Edge"
outeremm$zone="Outer"
inneremm$zone="Inner"

zoneemm=rbind(edgeemm,outeremm,inneremm)
zoneemm$Zone=as.factor(zoneemm$zone)
zoneemm$Zone=factor(zoneemm$Zone,levels=c("Edge","Outer","Inner"))
str(zoneemm)

zoneplot=ggplot(zoneemm,aes(x=factor(Treatment,level=treatorder),y=prob,fill=Zone))+
  geom_bar(stat="identity",width=0.5)+
  scale_fill_manual(values=c("grey80","grey50","grey20"))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.01))+
  xlab("")+ylab("LS Mean Proportion of Time per Tank Zone\n")+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=18,vjust=1,hjust=1,angle=45),
      axis.text.y=element_text(colour="black",size=16),
      axis.title.y=element_text(size=18),
      legend.title=element_text(size=16),
      legend.text=element_text(size=16))+
  annotate(geom="text",label=c("a","ab","a","bc","c"),
           x=c(1.4,2.42,3.4,4.42,5.4),y=0.97,size=7,colour="grey60")+
  annotate(geom="text",label=c("a","ab","ab","bc","c"),
           x=c(1.4,2.42,3.42,4.42,5.4),y=c(0.05,0.24,0.12,0.39,0.49),
           colour="grey40",size=7)+
  annotate(geom="text",label=c("–","ab","a","ab","b"),
           x=c(1.4,2.42,3.4,4.42,5.4),y=0.02,
           colour="grey10",size=7)
zoneplot
ggsave("Zoneplot.png",height=15,width=18,units="cm")


# Positive, neutral, negative rheotaxis proportions

behavdat1$PosProp=behavdat1$PosCount/(behavdat1$PosCount+behavdat1$NeuCount+behavdat1$NegCount)
behavdat1$NeuProp=behavdat1$NeuCount/(behavdat1$PosCount+behavdat1$NeuCount+behavdat1$NegCount)
behavdat1$NegProp=behavdat1$NegCount/(behavdat1$PosCount+behavdat1$NeuCount+behavdat1$NegCount)
str(behavdat1)

posmod=glm(PosProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(posmod)
posmeans=emmeans(posmod,pairwise~Treatment,type="response")
posmeans

posemm=as.data.frame(posmeans$emmeans)
posemm

ggplot(posemm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Positive Rheotaxis\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))

neumod=glm(NeuProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(neumod)
neumeans=emmeans(neumod,pairwise~Treatment,type="response",adjust="bonferroni")
neumeans

neuemm=as.data.frame(neumeans$emmeans)
neuemm

ggplot(neuemm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Neutral Rheotaxis\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))

negmod=glm(NegProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(negmod)
negmeans=emmeans(negmod,pairwise~Treatment,type="response",adjust="bonferroni")
negmeans

negemm=as.data.frame(negmeans$emmeans)
negemm

ggplot(negemm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time in Negative Rheotaxis\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=0.153,label="a",size=6)+
  annotate("text",x=2,y=0.126,label="a",size=6)+
  annotate("text",x=3,y=0.125,label="a",size=6)+
  annotate("text",x=4,y=0.078,label="b",size=6)+
  annotate("text",x=5,y=0.064,label="b",size=6)

table(behavdat1$NegCount,behavdat1$Treatment)


posemm$Orientation="+"
neuemm$Orientation="o"
negemm$Orientation="–"

directemm=rbind(posemm,neuemm,negemm)
directemm$Orientation=as.factor(directemm$Orientation)
directemm$Orientation=factor(directemm$Orientation,levels=c("+","o","–"))
str(directemm)

directplot=ggplot(directemm,aes(x=factor(Treatment,level=treatorder),y=prob,fill=Orientation))+
  geom_bar(stat="identity",width=0.5)+
  scale_fill_manual(values=c("grey80","grey50","grey20"))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.01))+
  xlab("")+ylab("LS Mean Proportion of Time in Each Orientation\n")+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=18,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(colour="black",size=16),
        axis.title.y=element_text(size=18),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16))+
  annotate(geom="text",label=c("a","a","a","b","b"),
           x=c(1.4,2.4,3.4,4.4,5.4),y=0.03,
           colour="grey10",size=7)
directplot

ggarrange(zoneplot,directplot,nrow=1,labels=c("A","B"),font.label=list(size=20))
ggsave("Zone & Orientation Plot 2.png",height=19,width=36,units="cm")



# Surfacing/near surfacing counts

behavdat1$SurfProp=behavdat1$Surfaces/60
str(behavdat1)

surfmod=glm(SurfProp~Treatment,data=behavdat1,family="quasibinomial")
Anova(surfmod,type="II",test.statistic="LR")
surfmeans=emmeans(surfmod,pairwise~Treatment,type="response",adjust="bonferroni")
surfmeans

surfemm=as.data.frame(surfmeans$emmeans)
surfemm

ggplot(surfemm,aes(x=factor(Treatment,level=treatorder),y=prob))+
  geom_bar(stat="identity",position="dodge")+
  ylab("LS Mean Proportion of Time At/Near Surface\n")+
  xlab("")+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12))+
  annotate("text",x=1,y=0.483,label="a",size=6)+
  annotate("text",x=2,y=0.324,label="a",size=6)+
  annotate("text",x=3,y=0.323,label="a",size=6)+
  annotate("text",x=4,y=0.174,label="b",size=6)+
  annotate("text",x=5,y=0.135,label="b",size=6)


##### Part II: Blood Data #####

blood=read.csv("C:/Users/CHR/Documents/FECPL/2024/Manitoba/Data/MB Part 2 Blood Data.csv")
blood$Treatment=as.factor(blood$Treatment)
blood$Time=as.factor(blood$Time)
blood$Plate=as.factor(blood$Plate)
str(blood)

hist(blood$CV)
mean(blood$CV)
min(blood$CV)
max(blood$CV)

plot(blood$Cortisol~blood$Glucose)
cor.test(blood$Cortisol,blood$Glucose,method="kendall")

glumod=glm(Glucose~Treatment*Time,data=blood[blood$Treatment!="PBaseline",],family="gaussian")
par(mfrow=c(2,2))
plot(glumod)
par(mfrow=c(1,1))
hist(resid(glumod))
shapiro.test(resid(glumod))

Anova(glumod,type="II",test.statistic="F")
emmeans(glumod,pairwise~Treatment|Time)
emmeans(glumod,pairwise~Time|Treatment)

Tapply(Glucose~Treatment+Time,data=blood,fun="min",na.action="na.omit")
Tapply(Glucose~Treatment+Time,data=blood,fun="max",na.action="na.omit")

set.seed(1)
ggplot(blood[blood$Treatment!="PBaseline",],aes(x=Treatment,y=Glucose,fill=Time))+
  geom_violin(scale="count")+
  geom_point(position=position_jitterdodge(dodge.width=0.9),size=3)+
  scale_fill_manual(values=c("grey25","grey50","grey80"))+
  xlab("")+ylab("Blood Glucose (mmol/l)\n")+
  scale_y_continuous(expand=c(0,0),limits=c(2,6))+
  theme_classic()+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=14))

BaseGluc=mean(blood$Glucose[blood$Treatment=="PBaseline"],na.rm=T)
BaseGlucStdev=sqrt(var(blood$Glucose[blood$Treatment=="PBaseline"],na.rm=T))
GlucUpperSE=BaseGluc+BaseGlucStdev
GlucLowerSE=BaseGluc-BaseGlucStdev
BaseGluc
GlucUpperSE
GlucLowerSE

max(blood$Glucose[blood$Treatment=="PBaseline"],na.rm=T)
min(blood$Glucose[blood$Treatment=="PBaseline"],na.rm=T)


lacmod=glm(LacDetect~Treatment*Time,data=blood[blood$Treatment!="PBaseline",],family="binomial")

Anova(lacmod,type="II",test.statistic="LR")
emmeans(lacmod,pairwise~Treatment|Time,type="response")
emmeans(lacmod,pairwise~Time|Treatment,type="response")

table(blood$LacDetect,blood$Treatment,blood$Time)


Tapply(Cortisol~Treatment+Time,data=blood,fun=function(x) length(unique(x)),na.action="na.omit")
Tapply(Cortisol~Treatment+Time,data=blood,fun="mean",na.action="na.omit")

#cortmod=glm(Cortisol~Treatment*Time,data=blood[blood$Treatment!="PBaseline",],family="gaussian")
#par(mfrow=c(2,2))
#plot(cortmod)
#par(mfrow=c(1,1))
#hist(resid(cortmod))
#shapiro.test(resid(cortmod))

cortmod2=glm(log10(Cortisol)~Treatment*Time,data=blood[blood$Treatment!="PBaseline",],family="gaussian")
par(mfrow=c(2,2))
plot(cortmod2)
par(mfrow=c(1,1))
hist(resid(cortmod2))
shapiro.test(resid(cortmod2))

Anova(cortmod2,type="II",test.statistic="F")

TrtOrder=c("Handling","TENS","Clove oil","MS-222")

# Raw values, then log below:
emmeans(cortmod2,pairwise~Time|Treatment,type="response")
emmeans(cortmod2,pairwise~Treatment|Time,type="response")

cmeans=emmeans(cortmod2,pairwise~Treatment|Time,type="response")
cmeans
cmeansdat=as.data.frame(cmeans$emmeans)
cmeansdat

BaseCort=mean(blood$Cortisol[blood$Treatment=="PBaseline"],na.rm=T)
BaseCortStdev=sqrt(var(blood$Cortisol[blood$Treatment=="PBaseline"],na.rm=T))
CortUpperSE=BaseCort+BaseCortStdev
CortLowerSE=BaseCort-BaseCortStdev
BaseCort
CortUpperSE
CortLowerSE
CortLowerSE=0
#Note that baseline cortisol lower cut-off is truncated at 0 (rather than showing -1.34) for graphing purposes

mean(blood$Cortisol[blood$Treatment=="PBaseline"],na.rm=T)
min(blood$Cortisol[blood$Treatment=="PBaseline"],na.rm=T)
max(blood$Cortisol[blood$Treatment=="PBaseline"],na.rm=T)

logcmeans=emmeans(cortmod2,pairwise~Treatment|Time)
logcmeans
logcmeansdat=as.data.frame(logcmeans$emmeans)
logcmeansdat

logBaseCort=mean(log10(blood$Cortisol[blood$Treatment=="PBaseline"]),na.rm=T)
logBaseCortStdev=sqrt(var(log10(blood$Cortisol[blood$Treatment=="PBaseline"]),na.rm=T))
logCortUpperSE=logBaseCort+logBaseCortStdev
logCortLowerSE=logBaseCort-logBaseCortStdev
logBaseCort
logCortUpperSE
logCortLowerSE

emmeans(cortmod2,pairwise~Treatment|Time)
emmeans(cortmod2,pairwise~Time|Treatment)

logcortplot=ggplot(logcmeansdat,aes(y=emmean,x=Treatment,fill=Time))+
  geom_bar(stat="identity",position="dodge",aes(x=factor(Treatment,level=TrtOrder)))+
  scale_fill_manual(values=c(rep(c("grey50","grey70","grey90"),times=4)))+
  geom_rect(ymin=logBaseCort-0.01,ymax=logBaseCort+0.01,xmin=0,xmax=5,
            fill="blue",alpha=0.03)+
  geom_errorbar(aes(x=Treatment,ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(width=0.9),width=0.5,lwd=0.8)+
  theme_classic()+scale_y_continuous(expand=c(0,0),limits=c(-1,3))+
  ylab("LS Mean log10 Plasma Cortisol Levels (ng/ml)")+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12))+
#  annotate(geom="text",x=c(0.7,1,1.3),y=c(1.982,1.260,0.335),label=c("a","b","b"),size=8)+
#  annotate(geom="text",x=c(1.7,2,2.3),y=c(2.198,1.215,0.748),label=c("a","b","b"),size=8)+
#  annotate(geom="text",x=c(2.7,3,3.3),y=c(1.573,0.819,0.372),label=c("a","b","b"),size=8)+
#  annotate(geom="text",x=c(3.7,4,4.3),y=c(1.369,0.732,0.137),label=c("a","b","c"),size=8)+
  annotate(geom="text",x=c(0.7,1.7,2.7,3.7),y=c(2.5,2.5,2.5,2.5),label=c("A","A","AB","B"),size=8)+
  annotate(geom="text",x=c(1,2,3,4),y=c(2.3,2.3,2.3,2.3),label=c("A","A","A","A"),size=8)+
  annotate(geom="text",x=c(1.3,2.3,3.3,4.3),y=c(2.1,2.1,2.1,2.1),label=c("A","A","A","A"),size=8)
logcortplot




cortplot=ggplot(blood[blood$Treatment!="PBaseline",],aes(y=Cortisol,x=Treatment,fill=Time))+
  geom_boxplot()+
  scale_fill_manual(values=c(rep(c("grey50","grey70","grey90"),times=4)))+
  geom_rect(ymin=CortLowerSE,ymax=CortUpperSE,xmin=0,xmax=6,
            fill="lightgreen",alpha=0.01)+
  theme_classic()+scale_y_continuous(expand=c(0,0),limits=c(0,130))+
  ylab("Plasma Cortisol Concentrations (ng/ml)")+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12))
cortplot


Tapply(Cortisol~Treatment+Time,data=blood,fun="min",na.action="na.omit")
Tapply(Cortisol~Treatment+Time,data=blood,fun="max",na.action="na.omit")
Tapply(Cortisol~Treatment+Time,data=blood,fun="mean",na.action="na.omit")

emmeans(cortmod2,pairwise~Time|Treatment,type="response")
emmeans(cortmod2,pairwise~Treatment|Time,type="response")

cortplot2=ggplot(blood[blood$Treatment!="PBaseline",],aes(y=Cortisol,x=factor(Treatment,level=TrtOrder),fill=Time))+
  geom_rect(ymin=CortLowerSE,ymax=CortUpperSE,xmin=0,xmax=6,
            fill="grey90",alpha=1)+
  geom_violin(scale="width",linewidth=0.7)+
  geom_point(position=position_jitterdodge(dodge.width=0.9,seed=2,jitter.width=0.4),size=1.5)+
  scale_fill_manual(values=c(rep(c("grey40","grey60","grey80"),times=4)))+
  theme_classic()+scale_y_continuous(expand=c(0,0),limits=c(0,111))+
  ylab("Plasma Cortisol Concentrations (ng/ml)\n")+
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(colour="black",size=12),
        axis.title.y=element_text(size=14))+
  xlab("")+
  annotate(geom="text",x=c(0.7,1.7,2.7,3.7),y=c(107, 107, 107, 107),label=c("AB","B","A","A"),size=8)+
  annotate(geom="text",x=c(0.7,1.7,2.7,3.7),y=c(74.70, 98.15, 37.76, 42.40),label=c("a","a","a","a"),size=8)+
  annotate(geom="text",x=c(1,2,3,4),y=c(20.72, 17.04, 17.28, 12.82),label=c("b","b","b","b"),size=8)+
  annotate(geom="text",x=c(1.3,2.3,3.3,4.3),y=c(30.34, 40.18, 11.65, 9.81),label=c("c","b","b","b"),size=8)+
  geom_point(data=cmeansdat,aes(x=Treatment,y=response,fill=Time),
             size=4,pch=1,stroke=2,position=position_dodge(width=0.9),
             alpha=0.4)
cortplot2
ggsave("Raw Plasma Cortisol Plot.png",cortplot2,dpi=300,width=20,height=15,units="cm")


baselinecort=data.frame(blood$Treatment[blood$Treatment=="PBaseline"],blood$Cortisol[blood$Treatment=="PBaseline"])
colnames(baselinecort)=c("Treatment","Cortisol")
time1cort=data.frame(blood$Treatment[blood$Time=="0.5"],blood$Cortisol[blood$Time=="0.5"])
colnames(time1cort)=c("Treatment","Cortisol")
time2cort=data.frame(blood$Treatment[blood$Time=="2"],blood$Cortisol[blood$Time=="2"])
colnames(time2cort)=c("Treatment","Cortisol")
time3cort=data.frame(blood$Treatment[blood$Time=="4"],blood$Cortisol[blood$Time=="4"])
colnames(time3cort)=c("Treatment","Cortisol")

time1cort=rbind(time1cort,baselinecort)
time2cort=rbind(time2cort,baselinecort)
time3cort=rbind(time3cort,baselinecort)

time1cort
time2cort
time3cort

time1mod=glm(log10(Cortisol)~Treatment,data=time1cort,family="gaussian")
par(mfrow=c(2,2))
plot(time1mod)
par(mfrow=c(1,1))
hist(resid(time1mod))
shapiro.test(resid(time1mod))

Anova(time1mod,type="II",test.statistic="F")
emmeans(time1mod,pairwise~Treatment,type="response")


time2mod=glm(log10(Cortisol)~Treatment,data=time2cort,family="gaussian")
par(mfrow=c(2,2))
plot(time2mod)
par(mfrow=c(1,1))
hist(resid(time2mod))
shapiro.test(resid(time2mod))

Anova(time2mod,type="II",test.statistic="F")
emmeans(time2mod,pairwise~Treatment,type="response")


time3mod=glm(log10(Cortisol)~Treatment,data=time3cort,family="gaussian")
par(mfrow=c(2,2))
plot(time3mod)
par(mfrow=c(1,1))
hist(resid(time3mod))
shapiro.test(resid(time3mod))

Anova(time3mod,type="II",test.statistic="F")
emmeans(time3mod,pairwise~Treatment,type="response")
