library(gdata)
library(epicalc)
library(grDevices)
library(lattice)
library(grid)
library(scatterplot3d)
library(arm)
library(Hmisc)
library(RColorBrewer) # creates nice color schemes
library(classInt)     # finds class intervals for continuous variables
library(gplots)
library(latticeExtra)
library(VIM)
library(reshape)
library(gmodels)
library(GGally)
library(gvlma)
library(MVA)
library(IC2)
library(class)
library(reshape2)
library(plyr)
library(tables)


#### Cargando bases de datos de tuberculosis
tb_world<-read.csv("world_tb_dataset.csv")
tab1(tb_world$Indicator)
tb_americas<-subset(tb_world,tb_world$WHO.region=="Americas")
data1=subset(tb_americas,Indicator %in% c("Incidence of tuberculosis (per 100 000 population per year)"))
data1<-data1[,c(4,2,5,6,7)]
data1 <- rename(data1, c(Numeric="ir_tb"))
data1 <- rename(data1, c(Low="ir_tb_low"))
data1 <- rename(data1, c(High="ir_tb_high"))
data2=subset(tb_americas,Indicator %in% c("Number of incident tuberculosis cases"))
data2<-data2[,c(4,2,5,6,7)]
data2 <- rename(data2, c(Numeric="tb_cases"))
data2 <- rename(data2, c(Low="tb_cases_low"))
data2 <- rename(data2, c(High="tb_cases_high"))
tb_lac <- merge(data1,data2,by=c("country","year"),all=TRUE)
tb_lac$country<-factor(tb_lac$country, label = c("Antigua and Barbuda","Argentina", "Bahamas", "Barbados","Belize", "Bolivia", "Brazil","Canada",
                                                 "Chile","Colombia","Costa Rica","Cuba","Dominica","Dominican Republic","Ecuador","El Salvador","Grenada",
                                                 "Guatemala","Guyana","Haiti","Honduras","Jamaica","Mexico","Nicaragua","Panama","Paraguay",
                                                 "Peru","Saint Kitts and Nevis","Saint Lucia","Saint Vincent and the Grenadines","Suriname","Trinidad and Tobago",
                                                 "United States of America","Uruguay","Venezuela"))
tb_lac_f<-subset(tb_lac,country %in% c("Argentina", "Bolivia", "Brazil","Chile","Colombia","Costa Rica","Cuba",
                                       "Dominican Republic","Ecuador","El Salvador","Guatemala","Guyana","Haiti","Honduras",
                                       "Mexico","Nicaragua","Panama","Paraguay","Peru",
                                       "Uruguay","Venezuela"))
tb_lac_f<-as.data.frame(subset(tb_lac_f, tb_lac_f$year==2013))


###############################################################################################
##################### Cargando bases de datos de diabetes #####################################
dm_world1<-read.csv("6th-Edition-Estimates_Update_2014.csv")
dm_world2<-read.csv("6th-Edition-Estimates_Update_2014_1.csv")
#dm_world1$IDF.Region[dm_world1$Country.territory=="Mexico"]<-"SACA"
#dm_world2$IDF.region[dm_world2$Country.territory=="Mexico"]<-"SACA"
#dm_world1$IDF.Region[dm_world1$Country.territory=="Haiti"]<-"SACA"
#dm_world2$IDF.region[dm_world2$Country.territory=="Haiti"]<-"SACA"
dm_americas1a<-subset(dm_world1,dm_world1$IDF.Region==c("SACA"))
dm_americas2a<-subset(dm_world2,dm_world2$IDF.region==c("SACA"))
dm_americas1b<-subset(dm_world1,dm_world1$IDF.Region==c("NAC"))
dm_americas2b<-subset(dm_world2,dm_world2$IDF.region==c("NAC"))
dm_lac1a<-dm_americas1a[,c(1,3,4,5,6,7,8,9,10)]
dm_lac2a<-dm_americas2a[,c(1,3,4,5,6,7,8)]
dm_lac1b<-dm_americas1b[,c(1,3,4,5,6,7,8,9,10)]
dm_lac2b<-dm_americas2b[,c(1,3,4,5,6,7,8)]
dm_lac1a$Country.territory <- revalue(dm_lac1a$Country.territory, c("Bolivia (Plurinational State of)"="Bolivia", "Venezuela (Bolivarian Republic of)"="Venezuela"))
dm_lac2a$Country.territory <- revalue(dm_lac2a$Country.territory, c("Bolivia (Plurinational State of)"="Bolivia", "Venezuela (Bolivarian Republic of)"="Venezuela"))
dm_lac1b$Country.territory <- revalue(dm_lac1b$Country.territory, c("Bolivia (Plurinational State of)"="Bolivia", "Venezuela (Bolivarian Republic of)"="Venezuela"))
dm_lac2b$Country.territory <- revalue(dm_lac2b$Country.territory, c("Bolivia (Plurinational State of)"="Bolivia", "Venezuela (Bolivarian Republic of)"="Venezuela"))
dm_lac_f1<-merge(dm_lac1a,dm_lac2a,by=c("Country.territory"),all=TRUE)
dm_lac_f2<-merge(dm_lac1b,dm_lac2b,by=c("Country.territory"),all=TRUE)
dm_lac_f<-as.data.frame(rbind(dm_lac_f1,dm_lac_f2))
dm_lac_f <- rename(dm_lac_f, c(Country.territory="country"))
dm_lac_f<-subset(dm_lac_f,country %in% c("Argentina", "Bolivia", "Brazil","Chile","Colombia","Costa Rica","Cuba",
                                         "Dominican Republic","Ecuador","El Salvador","Guatemala","Guyana","Haiti","Honduras",
                                         "Mexico","Nicaragua","Panama","Paraguay","Peru",
                                         "Uruguay","Venezuela"))
dm_lac_ff<-data.frame(country=dm_lac_f$country,prev_dm=dm_lac_f$Diabetes.national.prevalence....,db_cases=dm_lac_f$Diabetes.cases..20.79..in.1000s.x,dm_expenditure=dm_lac_f$Mean.diabetes.related.expenditure.per.person.with.diabetes..USD.)
tb_dm_lac<-merge(tb_lac_f,dm_lac_ff,by=c("country"),all=TRUE)


### Calculo del RAP ###
Pe<-(tb_dm_lac$prev_dm)/100
Pe
RRe<-matrix(c(3.11,2.27,4.26 ),3)
RRe
RRe[1]
ORe<-matrix(c(1.2,7.8),2)
ORe
ORe[1]
f<-function(Pe,RRe) {
  (Pe*(RRe-1)/(1 + Pe*(RRe-1)))*100
}
tb_dm_lac$RAP<-round(f(Pe,RRe[1]),2)
list(tb_dm_lac$country,tb_dm_lac$RAP)
tb_dm_lac$RAP_low<-round(f(Pe,RRe[2]),2)
list(tb_dm_lac$country,tb_dm_lac$RAP_low)
tb_dm_lac$RAP_high<-round(f(Pe,RRe[3]),2)
list(tb_dm_lac$country,tb_dm_lac$RAP_high)
tb_dm_lac$RAP_OR_low<-round(f(Pe,ORe[1]),2)
list(tb_dm_lac$country,tb_dm_lac$RAP_OR_low)
tb_dm_lac$RAP_OR_high<-round(f(Pe,ORe[2]),2)
list(tb_dm_lac$country,tb_dm_lac$RAP_OR_high)
tb_dm_lac$tb_cases_atr_dm<-round((tb_dm_lac$tb_cases*tb_dm_lac$RAP)/100,0)
tb_dm_lac$tb_cases_atr_dm_low<-round((tb_dm_lac$tb_cases*tb_dm_lac$RAP_low)/100,0)
tb_dm_lac$tb_cases_atr_dm_high<-round((tb_dm_lac$tb_cases*tb_dm_lac$RAP_high)/100,0)
plot(tb_dm_lac$RAP,tb_dm_lac$ir_tb)
plot(tb_dm_lac$ir_tb,tb_dm_lac$RAP)
summ(tb_dm_lac$RAP,graph=FALSE)
summ(tb_dm_lac$RAP_low,graph=FALSE)
summ(tb_dm_lac$RAP_high,graph=FALSE)
tb_dm_lac_sort<-tb_dm_lac[order(tb_dm_lac$prev_dm),]

quartz(width=10, height=6, pointsize=10)
plot(tb_dm_lac_sort$prev_dm, tb_dm_lac_sort$RAP, ylim = c(0,30), type = "l", xlab="Prevalencia de Diabetes (%)", ylab="Fracción Atribuible Poblacional (%)")
#make polygon where coordinates start with lower limit and
# then upper limit in reverse order
polygon(c(tb_dm_lac_sort$prev_dm,rev(tb_dm_lac_sort$prev_dm)),c(tb_dm_lac_sort$RAP_low,rev(tb_dm_lac_sort$RAP_high)),col = "grey75", border = FALSE)
lines(tb_dm_lac_sort$prev_dm, tb_dm_lac_sort$RAP, lwd = 2)
#add red lines on borders of p`olygon
lines(tb_dm_lac_sort$prev_dm, tb_dm_lac_sort$RAP_high, col="red",lty=2)
lines(tb_dm_lac_sort$prev_dm, tb_dm_lac_sort$RAP_low, col="red",lty=2)



data.table(Country=tb_dm_lac$country,TB.incidence=tb_dm_lac$ir_tb,TB.cases=tb_dm_lac$tb_cases,
           Prev.DM=tb_dm_lac$prev_dm,DM.cases=tb_dm_lac$db_cases,PAF=tb_dm_lac$RAP,PAF.low=tb_dm_lac$RAP_low,
           PAF.high=tb_dm_lac$RAP_high,TB.case.atrb.DM=tb_dm_lac$tb_cases_atr_dm,TB.case.atrb.DM.low=tb_dm_lac$tb_cases_atr_dm_low,TB.case.atrb.DM.high=tb_dm_lac$tb_cases_atr_dm_high)




