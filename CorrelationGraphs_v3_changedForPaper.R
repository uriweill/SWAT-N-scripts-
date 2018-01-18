########################################################################################
#    Program name:     CorrelationGraphs
#    Program Owner:    Nofar Harpaz
#    Date created:     2017-Jan-11
#    Description:      This program creates corelation graphs comparing correlation between 
#                      the generic full genome (Nop1-N terminal), 
#                      native (Seamless) and cherry (Tef2-Nterminal-Cherry) libraries.
#
#                      This script is based on Silvia's program:
#                      CompareOldDatavs5halftimeScaleSwapV2_2015_DifferentAutoForSP
#
#    Updated by:
#    Date of update:
#    Reason for update:
#
########################################################################################

###################### Set the working enviorment #######################

#clear the workspace
rm(list=ls())

# Define variables:
#myDir='D:/Nofar/Maya'
#myDir='E:/Nofar/Scripts/SWAT_URI'
myDir='C:/MyComputerAtMSlabBackup/Scripts/SWAT_URI'
InputFile='SWAT_GFP FG Vs dtabases V4.csv';
OutputFile='3LibGraphs_v3_changed.pdf';

#set directory
setwd(myDir)

#insall packages - activate below comments if packages weren't installed in this R version

#install.packages('plotrix')
#install.packages('RColorBrewer')
#install.packages('gplots')
#install.packages('car')
#install.packages('sp')
#install.packages('psych')
#install.packages(('corrplot'))


library('plotrix')
library('RColorBrewer')
library('gplots')
#library('car')
library('sp')
library('psych')



###################### Load the data ####################################

data1<-read.csv(InputFile,header=TRUE, na.strings=c('#N/A','NA'))


###################### Manipulate the data ##############################

## Preprocess the data:

#Adding a plate number column
data1$Plate=as.factor(as.character(strsplit(as.character(data1$X384.Plate_Row_Col),"[A-Z]\\d+", fixed=FALSE, perl=FALSE)))
#Change 0 score of intensity to NA
data1$Generic.GfpIntensity.FG.[which(data1$Generic.GfpIntensity.FG.==0)]=NA
data1$Native.GfpIntensity.FG.[which(data1$Native.GfpIntensity.FG.==0)]=NA
data1$Tef2pr.mCherryIntensity.FG.[which(data1$Tef2pr.mCherryIntensity.FG.==0)]=NA
#Delete all data1$Generic.GfpIntensity.FG.==NA from table
data1=data1[which(!is.na(data1$Generic.GfpIntensity.FG.)),]

#########################################################################
#                         Local function                                #
#########################################################################
#func name  - reduce_AF
#func des   - this function calculates the mean of the autofluorescent wells per libraray per plate, 
#             and reduce it from library score. For libraray Tef2pr.mCherryIntensity.FG. scores, the mean is 
#             calculated differently: the mean of plate 16 is reduced from this plate only, and the mean caculated 
#             from plates 1-4, is reduced from all plates except #16 (upon uri's request)
#function arguments:
#DATA       - the name of your dataframe
#COL        - the name of your column (the libraray scores), e.g: COL='Generic.GfpIntensity.FG.'
#CONTROL    - a char vector with the names identify of the control wells as they appear in the ORF column, 
#             eg. c('Control-his3d-GFPdC','URA+')
#
###########################################################################

reduce_AF=function(DATA,COL,CONTROL) {
  
  #producing a smaller table, of the requested controls
  A=DATA[which(DATA$ORF %in% CONTROL),]
  #find the index of your column of interest
  ind=which(colnames(A)==COL)
  if (COL!='Tef2pr.mCherryIntensity.FG.') {
    #calculate mean per plate
    AF_mean<-as.matrix(tapply(A[,ind], A$Plate, mean,na.rm=TRUE))
    #merge
    data.link=merge(x=DATA, y=AF_mean, by.x='Plate', by.y='row.names',sort=FALSE)
    #find the index of your column of interest in the merged data
    ind2=which(colnames(data.link)==COL)
    #substract AF mean from Intensity
    data.link[,ind2]=(data.link[,ind2])-(data.link$V1)
    #remove th extra column
    data.link=data.link[,which(colnames(data.link)!='V1')]
    return (data.link)
  }
  else {
    #calculate mean per plate 1,2,3,4, and reduce it from all plates except plate 16
    DATA[which(DATA$Plate!=16),ind]=DATA[which(DATA$Plate!=16),ind]-mean(A[which(A$Plate %in% c(1,2,3,4)),ind], na.rm=TRUE)
    #calculate mean per plate 16 and reduce if from only plate 16
    DATA[which(DATA$Plate==16),ind]=DATA[which(DATA$Plate==16),ind]-mean(A[which(A$Plate==16),ind], na.rm=TRUE)
    return (DATA)

  }
  
}
######################################################################################

## Reduce autofluorescence mean from libraries scores

data2=reduce_AF(DATA=data1,COL='Generic.GfpIntensity.FG.',CONTROL=c('Control-his3d-GFPdC','URA+'))
data2=reduce_AF(DATA=data2,COL='Native.GfpIntensity.FG.',CONTROL=c('Control-his3d-GFPdC','his3?1::SWAT / SWAT-SP','his3?1::SWAT-MTS'))  
data2=reduce_AF(DATA=data2,COL='Tef2pr.mCherryIntensity.FG.',CONTROL=c('Control-his3d-GFPdC'))  

#Change <1 score of intensity to NA
data2$Generic.GfpIntensity.FG.[which(data2$Generic.GfpIntensity.FG.<=1)]=NA
data2$Native.GfpIntensity.FG.[which(data2$Native.GfpIntensity.FG.<=1)]=NA
data2$Tef2pr.mCherryIntensity.FG.[which(data2$Tef2pr.mCherryIntensity.FG.<=1)]=NA
#Delete all data1$Generic.GfpIntensity.FG.==NA from table
data2=data2[which(!is.na(data2$Generic.GfpIntensity.FG.)),]


## Extract only rows where ORF column value starts with a 'Y'

data3=data2[which(substr(data2$ORF,1,1)=='Y'),]
  

## sort data, ascending values of Generic.GfpIntensity.FG.
order_ind<-order(data3$Generic.GfpIntensity.FG.) 
data4_s=data3[order_ind,]





#########################################################################
#                         Local functions for plots                     #
#########################################################################
#func name  - add_text_to_point
#func des   - this function add a "|" sign on a plot, at the lowest value that is bigger than the input value.
#             it also add a text specifying the input value at the next to the point on plot
#function arguments:
#DATA.COL   - the name of your vector (could be data$col, or a vector vraiable). The function expect
#             to recieve a vector of numeric values,sorted in ascending order
#VAL        - the value of which you want to compare
#DELTA      - when printing the text on the plot, the function uses this delat to define the hieght of the text 
#             comparing to the y value of the point
#
###########################################################################


add_text_to_point=function(DATA.COL,VAL,DELTA){
  ind=which(DATA.COL>VAL)
  value=DATA.COL[ind[1]]
  text(ind[1],value,labels='|', font=2)
  text(ind[1],value+DELTA, labels=VAL, font=2)
  print(c(ind[1],value))
}

#########################################################################
#func name  - change_ticks_to_int
#func des   - this function add x and y ticks to log axes, where the values are integers. 
#             This function should by called after a plot command, where log='xy'. In the plot command,
#             don't forget to erase the x and y ticks (xaxt = 'n', yaxt='n').
#             
#
#This function has no arguments.
#
#See here example of how to call the function:
#plot(data4_s$Native.GfpIntensity.FG., data4_s$Generic.GfpIntensity.FG., log='xy',xlim=c(1e-02,1e+03), ylim=c(1e-02,1e+03), xaxt = 'n', yaxt='n')
#change_ticks_to_int()
#
###########################################################################
change_ticks_to_int=function(){
  myTicks1 = axTicks(1)
  myTicks2 = axTicks(2)
  axis(1, at = myTicks1, labels = formatC(myTicks1, format = 'd'))
  axis(2, at = myTicks2, labels = formatC(myTicks2, format = 'd'))
}

change_ticks_to_dec=function(){
  myTicks1 = axTicks(1)
  myTicks2 = axTicks(2)
  axis(1, at = myTicks1, labels = formatC(myTicks1, format = 'f'))
  axis(2, at = myTicks2, labels = formatC(myTicks2, format = 'f'))
}

#########################################################################
#func name  - find_log_scale_lim
#func des   - this function gets the numeric vectors used as x and y in your plot (when a log scale is used but the data is
#             not log), and returns the value of your min and max. If your min is <0, it'll return the lowest >0 value.
#             This is because log2 of negative number are irrational. 
#             
#
#This function has 2 arguments: x and y - which are numeric vectors used as x and y in your plot
#
#
###########################################################################
find_log_scale_lim=function(x,y){
  x1=x[x>0]
  y1=y[y>0]
  abs_min=min(c(min(x1, na.rm=T), min(y1, na.rm=T)))
  abs_max=max(c(max(x1, na.rm=T), max(y1, na.rm=T)))
  return (c(abs_min,abs_max))
}

#########################################################################
#func name  - find_log_data_lim
#func des   - Same as find_log_scale_lim. Only this time x and y are the values after log2() was performed.
#             In this case, it is possible that negative values will exist.
#             
#
#This function has 2 arguments: x and y - which are numeric vectors used as x and y in your plot
#
###########################################################################

find_log_data_lim=function(x,y){
  x1=x
  y1=y
  abs_min=min(c(min(x1, na.rm=T), min(y1, na.rm=T)))
  abs_max=max(c(max(x1, na.rm=T), max(y1, na.rm=T)))
  return (c(abs_min,abs_max))
}

#########################################################################
#func name  - Cal_R_native_q / Cal_R_generic_q / Cal_R_tef2_q
#func des   - This function calculate the R sperman correlation between two columns for : the entire data (data4_s) 
#             and Native based q1-q4
#             
#
#This function has 2 arguments: x and y - which are strings of the coulmns names for which you wish to calculate R 
#
###########################################################################

Cal_R_native_q=function(x,y){
  indx=which(colnames(data4_s)==x)
  indy=which(colnames(data4_s)==y)
  #testing correllation spearman, based on Native quantile
  CORC.N.=cor.test(data4_s[,indx],data4_s[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.2=round(CORC.N.$estimate, 2)
  
  indx=which(colnames(data_Nat_q1)==x)
  indy=which(colnames(data_Nat_q1)==y)
  CORC.N.Q1=cor.test(data_Nat_q1[,indx],data_Nat_q1[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q1.2=round(CORC.N.Q1$estimate, 2)
  CORC.N.Q2=cor.test(data_Nat_q2[,indx],data_Nat_q2[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q2.2=round(CORC.N.Q2$estimate, 2)
  CORC.N.Q3=cor.test(data_Nat_q3[,indx],data_Nat_q3[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q3.2=round(CORC.N.Q3$estimate, 2)
  CORC.N.Q4=cor.test(data_Nat_q4[,indx],data_Nat_q4[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q4.2=round(CORC.N.Q4$estimate, 2)
  return (c(CORC.N.2,CORC.N.Q1.2,CORC.N.Q2.2,CORC.N.Q3.2,CORC.N.Q4.2))
  
}


Cal_R_generic_q=function(x,y){
  indx=which(colnames(data4_s)==x)
  indy=which(colnames(data4_s)==y)
  #testing correllation spearman, based on Native quantile
  CORC.N.=cor.test(data4_s[,indx],data4_s[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.2=round(CORC.N.$estimate, 2)
  
  indx=which(colnames(data_Gen_q1)==x)
  indy=which(colnames(data_Gen_q1)==y)
  CORC.N.Q1=cor.test(data_Gen_q1[,indx],data_Gen_q1[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q1.2=round(CORC.N.Q1$estimate, 2)
  CORC.N.Q2=cor.test(data_Gen_q2[,indx],data_Gen_q2[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q2.2=round(CORC.N.Q2$estimate, 2)
  CORC.N.Q3=cor.test(data_Gen_q3[,indx],data_Gen_q3[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q3.2=round(CORC.N.Q3$estimate, 2)
  CORC.N.Q4=cor.test(data_Gen_q4[,indx],data_Gen_q4[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q4.2=round(CORC.N.Q4$estimate, 2)
  return (c(CORC.N.2,CORC.N.Q1.2,CORC.N.Q2.2,CORC.N.Q3.2,CORC.N.Q4.2))
  
} 

Cal_R_tef2_q=function(x,y){
  indx=which(colnames(data4_s)==x)
  indy=which(colnames(data4_s)==y)
  #testing correllation spearman, based on Native quantile
  CORC.N.=cor.test(data4_s[,indx],data4_s[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.2=round(CORC.N.$estimate, 2)
  
  indx=which(colnames(data_Tef_q1)==x)
  indy=which(colnames(data_Tef_q1)==y)
  CORC.N.Q1=cor.test(data_Tef_q1[,indx],data_Tef_q1[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q1.2=round(CORC.N.Q1$estimate, 2)
  CORC.N.Q2=cor.test(data_Tef_q2[,indx],data_Tef_q2[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q2.2=round(CORC.N.Q2$estimate, 2)
  CORC.N.Q3=cor.test(data_Tef_q3[,indx],data_Tef_q3[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q3.2=round(CORC.N.Q3$estimate, 2)
  CORC.N.Q4=cor.test(data_Tef_q4[,indx],data_Tef_q4[,indy],method="spearman",use="pairwise.complete.obs")
  CORC.N.Q4.2=round(CORC.N.Q4$estimate, 2)
  return (c(CORC.N.2,CORC.N.Q1.2,CORC.N.Q2.2,CORC.N.Q3.2,CORC.N.Q4.2))
  
} 

###################### PDF definitions ##################################


# windows()
# 
e<-paste(myDir,"/",OutputFile,sep="")

graphics.off()  
pdf(e,width=7,height=5.8,family='Helvetica', paper='A4')
# 
# par(mar=c(4.8,4,3,1))
# par(mfrow=c(2,3))
# par(family="Helvetica")
####################### Plot ###########################################

#General
#storing my old par()
old_par=par()
#the defult is par(pty="m")
par(pty="s") #defining a square graph

#creating datasets of quantile according to each library 
#creating data of native q1 - q4
order_ind<-order(data3$Native.GfpIntensity.FG. ) 
data5_s_N=data3[order_ind,]
Nat_stat=summary(data5_s_N$Native.GfpIntensity.FG. )
Nat_quantile=quantile(data5_s_N$Native.GfpIntensity.FG.,na.rm=T, names=T)
data_Nat_q1=data5_s_N[which(!is.na(data5_s_N$Native.GfpIntensity.FG.) & data5_s_N$Native.GfpIntensity.FG.<=Nat_quantile[2]),]
data_Nat_q2=data5_s_N[which(!is.na(data5_s_N$Native.GfpIntensity.FG.) & Nat_quantile[2]<=data5_s_N$Native.GfpIntensity.FG. & data5_s_N$Native.GfpIntensity.FG.<=Nat_quantile[3]),]
data_Nat_q3=data5_s_N[which(!is.na(data5_s_N$Native.GfpIntensity.FG.) & Nat_quantile[3]<=data5_s_N$Native.GfpIntensity.FG. & data5_s_N$Native.GfpIntensity.FG.<=Nat_quantile[4]),]
data_Nat_q4=data5_s_N[which(!is.na(data5_s_N$Native.GfpIntensity.FG.) & data5_s_N$Native.GfpIntensity.FG.>=Nat_quantile[4]),]
#creating data of Generic q1 - q4
Gen_stat=summary(data4_s$Generic.GfpIntensity.FG.  )
Gen_quantile=quantile(data4_s$Generic.GfpIntensity.FG.,na.rm=T, names=T)
data_Gen_q1=data5_s_N[which(!is.na(data4_s$Generic.GfpIntensity.FG.) & data4_s$Generic.GfpIntensity.FG.<=Gen_quantile[2]),]
data_Gen_q2=data5_s_N[which(!is.na(data4_s$Generic.GfpIntensity.FG.) & Gen_quantile[2]<=data4_s$Generic.GfpIntensity.FG. & data4_s$Generic.GfpIntensity.FG.<=Gen_quantile[3]),]
data_Gen_q3=data5_s_N[which(!is.na(data4_s$Generic.GfpIntensity.FG.) & Gen_quantile[3]<=data4_s$Generic.GfpIntensity.FG. & data4_s$Generic.GfpIntensity.FG.<=Gen_quantile[4]),]
data_Gen_q4=data5_s_N[which(!is.na(data4_s$Generic.GfpIntensity.FG.) & data4_s$Generic.GfpIntensity.FG.>=Gen_quantile[4]),]
#creating data of Tef2 q1 - q4
order_ind<-order(data3$Tef2pr.mCherryIntensity.FG. ) 
data5_s_T=data3[order_ind,]
Tef_stat=summary(data5_s_T$Tef2pr.mCherryIntensity.FG. )
Tef_quantile=quantile(data5_s_T$Tef2pr.mCherryIntensity.FG. ,na.rm=T, names=T)
data_Tef_q1=data5_s_T[which(!is.na(data5_s_T$Tef2pr.mCherryIntensity.FG.) & data5_s_T$Tef2pr.mCherryIntensity.FG.<=Tef_quantile[2]),]
data_Tef_q2=data5_s_T[which(!is.na(data5_s_T$Tef2pr.mCherryIntensity.FG.) & Tef_quantile[2]<=data5_s_T$Tef2pr.mCherryIntensity.FG. & data5_s_T$Tef2pr.mCherryIntensity.FG.<=Tef_quantile[3]),]
data_Tef_q3=data5_s_T[which(!is.na(data5_s_T$Tef2pr.mCherryIntensity.FG.) & Tef_quantile[3]<=data5_s_T$Tef2pr.mCherryIntensity.FG. & data5_s_T$Tef2pr.mCherryIntensity.FG.<=Tef_quantile[4]),]
data_Tef_q4=data5_s_T[which(!is.na(data5_s_T$Tef2pr.mCherryIntensity.FG.) & data5_s_T$Tef2pr.mCherryIntensity.FG.>=Tef_quantile[4]),]



#plot 1 :SWAT libraries signal intensity
x=seq(from =1,to= length(data4_s$Generic.GfpIntensity.FG. ),by=1)
plot(x,data4_s$Generic.GfpIntensity.FG.,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="SWAT signal intensity (a.u.)" , font.lab=2, main="SWAT libraries protein abundance")
lines(x,data4_s$Tef2pr.mCherryIntensity.FG. ,type="p",col="aquamarine2",pch=19,cex=0.35)
lines(x,data4_s$Native.GfpIntensity.FG. ,type="p",col="darkgray",pch=19,cex=0.35)
lines(x,data4_s$Generic.GfpIntensity.FG.,type="p",col="darkmagenta",pch=19,cex=0.35)
legend("topleft",c("NOP1pr-GFP","Native-GFP", "TEF2pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c("darkmagenta","darkgray","aquamarine2"),bty="n",text.font=1) # 


#plot 2: NOP1-pr-GFP library signal intensity
plot(x,data4_s$Generic.GfpIntensity.FG.,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="NOP1pr-GFP signal intensity (a.u.)" , font.lab=2, main="SWAT-NOP1pr-GFP library protein abundance")
#add_text_to_point(data4_s$Generic.GfpIntensity.FG.,0,100)
add_text_to_point(data4_s$Generic.GfpIntensity.FG.,10,100) 
add_text_to_point(data4_s$Generic.GfpIntensity.FG.,100,100) 

y=log2(data4_s$Generic.GfpIntensity.FG.)
plot(x,y,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="NOP1pr-GFP signal intensity (a.u.)" , font.lab=2, main="SWAT-NOP1pr-GFP library protein abundance", yaxt='n')
myTicks=axTicks(2)
# labels <- sapply(myTicks,function(i)
#   as.expression(bquote(2^ .(i)))
#   )
labels <- sapply(myTicks,function(i) 2^ i)
axis(2,at=myTicks,labels=labels)

#plot 3: Native-pr-GFP library signal intensity
order_ind<-order(data3$Native.GfpIntensity.FG.) 
Native_sorted=data3[order_ind,colnames(data3)=='Native.GfpIntensity.FG.']
Native_sorted=Native_sorted[which(!is.na(Native_sorted))]
x=seq(from =1,to= length(Native_sorted ),by=1)
plot(x,Native_sorted,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="Native-GFP signal intensity (a.u.)" , font.lab=2, main="SWAT-Native-GFP library protein abundance")
#add_text_to_point(Native_sorted,0,200)
add_text_to_point(Native_sorted,10,200)
add_text_to_point(Native_sorted,100,200)

y=log2(Native_sorted)
plot(x,y,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="Native-GFP signal intensity (a.u.)" , font.lab=2, main="SWAT-Native-GFP library protein abundance", yaxt='n')
myTicks=axTicks(2)
labels <- sapply(myTicks,function(i) 2^ i)
axis(2,at=myTicks,labels=labels)

#plot 4: TEF2-pr-mCherry library signal intensity
order_ind<-order(data3$Tef2pr.mCherryIntensity.FG.) 
Tef2_sorted=data3[order_ind,colnames(data3)=='Tef2pr.mCherryIntensity.FG.']
Tef2_sorted=Tef2_sorted[which(!is.na(Tef2_sorted))]
x=seq(from =1,to= length(Tef2_sorted ),by=1)
plot(x,Tef2_sorted,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="TEF2pr-mCherry signal intensity (a.u.)" , font.lab=2, main="SWAT-TEF2-pr-mCherry library protein abundance")
#add_text_to_point(Tef2_sorted,0,100)
add_text_to_point(Tef2_sorted,10,100)
add_text_to_point(Tef2_sorted,100,100)

y=log2(Tef2_sorted)
plot(x,y,type="p",col="darkmagenta",pch=19,cex=0.35,xlab="Strains" ,ylab="TEF2pr-mCherry signal intensity (a.u.)" , font.lab=2, main="SWAT-TEF2-pr-mCherry library protein abundance", yaxt='n')
myTicks=axTicks(2)
labels <- sapply(myTicks,function(i) 2^ i)
axis(2,at=myTicks,labels=labels)

#Labels for all axes in the following plots (log2)

lab_nop_log=expression("log"[2]*"(SWAT-NOP1pr-GFP signal intensity (a.u.))")
lab_nat_log=expression("log"[2]*"(SWAT-Native-GFP signal intensity (a.u.))")
lab_tef_log=expression("log"[2]*"(SWAT-TEF2pr-mCherry signal intensity (a.u.))")
lab_all_3=expression("log"[2]*"(signal intensity (a.u.))")
lab_cGFP_log=expression("log"[2]*"(C' GFP signal intensity (a.u.))")
lab_ribo_log=expression("log"[2]*"(translation rate)")
lab_mas_log=expression("log"[2]*"(protein abundance)")
lab_mRNA_log=expression("log"[2]*"(mRNA abundance)")
lab_FACS_log=expression("log"[2]*"(C' GFP signal intensity (a.u.))")
lab_half_log=expression("log[2]"*"(protein half-life (min))")

#plot 5: NOP1-pr-GFP vs Native-pr-GFP -all
#testing correllation spearman, based on Native quantile
myRs=Cal_R_native_q('Generic.GfpIntensity.FG.','Native.GfpIntensity.FG.')
#plot and add text
#calculate axes max and min
x=log2(data4_s$Native.GfpIntensity.FG.)
y=log2(data4_s$Generic.GfpIntensity.FG.)
axes_lim=find_log_data_lim(x,y)
#plot(data4_s$Native.GfpIntensity.FG., data4_s$Generic.GfpIntensity.FG., type="p",pch=19,cex=0.35,ylab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-Native-pr-GFP", log='xy')
#plot(data4_s$Native.GfpIntensity.FG., data4_s$Generic.GfpIntensity.FG., type="p",pch=19,cex=0.35,ylab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-Native-pr-GFP", log='xy',xlim=c(1e-02,1e+03), ylim=c(1e-02,1e+03))
#plot(data4_s$Native.GfpIntensity.FG., data4_s$Generic.GfpIntensity.FG., type="p",pch=19,cex=0.35,ylab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-Native-pr-GFP", log='xy',xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
plot(x, y, type="p",pch=19,cex=0.35,ylab='SWAT-NOP1pr-GFP signal intensity (a.u.)' ,xlab='SWAT-Native-GFP signal intensity (a.u.)' , font.lab=2, main="Comparison of protein abundance between SWAT-NOP1pr-GFP \nvs SWAT-Native-GFP",xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]), yaxt='n', xaxt='n')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,11.5,labels=rall, font=2)
# text(10,2.6,labels=rq1, font=2)
# text(10,2.2,labels=rq2, font=2)
# text(10,1.8,labels=rq3, font=2)
# text(9.9,1.4,labels=rq4, font=2)
myTicks=axTicks(2)
labels <- sapply(myTicks,function(i) 2^ i)
axis(2,at=myTicks,labels=labels)
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)

#plot 6: NOP1-pr-GFP vs TEF2-pr-mCherry - all
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('Generic.GfpIntensity.FG.','Tef2pr.mCherryIntensity.FG.')
#plot and add text
x=log2(data4_s$Generic.GfpIntensity.FG.)
y=log2(data4_s$Tef2pr.mCherryIntensity.FG.)
axes_lim=find_log_data_lim(x,y)
#plot(data4_s$Generic.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG. ,type="p",pch=19,cex=0.35,xlab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,ylab="SWAT-TEF2-pr-mCherry protein abundance \n(mCherry intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-TEF2-pr-mCherry",log='xy')
#plot(data4_s$Generic.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG. ,type="p",pch=19,cex=0.35,xlab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,ylab="SWAT-TEF2-pr-mCherry protein abundance \n(mCherry intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-TEF2-pr-mCherry",log='xy', xlim=c(1e-03,1e+03))
#plot(data4_s$Generic.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG. ,type="p",pch=19,cex=0.35,xlab="SWAT-NOP1pr-GFP protein abundance (GFP intensity (a.u.))" ,ylab="SWAT-TEF2-pr-mCherry protein abundance \n(mCherry intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-NOP1-pr-GFP vs SWAT-TEF2-pr-mCherry",log='xy', xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
plot(x,y ,type="p",pch=19,cex=0.35,xlab=lab_nop_log ,ylab=lab_tef_log , font.lab=2, main="Comparison of protein abundance between SWAT-NOP1pr-GFP \nvs SWAT-TEF2pr-mCherry", xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))  
text(1,10,labels=rall, font=2)
# text(1.5,9.6,labels=rq1, font=2)
# text(1.5,9.2,labels=rq2, font=2)
# text(1.5,8.8,labels=rq3, font=2)
# text(1.5,8.4,labels=rq4, font=2)




#plot 7: TEF2-pr-mCherry vs Native-pr-GFP -all
#testing correllation spearman, based on Native quantile
myRs=Cal_R_native_q('Tef2pr.mCherryIntensity.FG.','Native.GfpIntensity.FG.')
#plot and add text
x=log2(data4_s$Native.GfpIntensity.FG.)
y=log2(data4_s$Tef2pr.mCherryIntensity.FG.)
axes_lim=find_log_data_lim(x,y)
#plot(data4_s$Native.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG.,type="p",pch=19,cex=0.35,ylab="SWAT-TEF2-pr-mCherry protein abundance (mCherry intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-TEF2-pr-mCherry vs SWAT-Native-pr-GFP", log='xy')
#plot(data4_s$Native.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG.,type="p",pch=19,cex=0.35,ylab="SWAT-TEF2-pr-mCherry protein abundance (mCherry intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-TEF2-pr-mCherry vs SWAT-Native-pr-GFP", log='xy', xlim=c(1e-03,1e+03))
#plot(data4_s$Native.GfpIntensity.FG.,data4_s$Tef2pr.mCherryIntensity.FG.,type="p",pch=19,cex=0.35,ylab="SWAT-TEF2-pr-mCherry protein abundance (mCherry intensity (a.u.))" ,xlab="SWAT-Native-GFP protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of SWAT-TEF2-pr-mCherry vs SWAT-Native-pr-GFP", log='xy', xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
plot(x,y,type="p",pch=19,cex=0.35,ylab=lab_tef_log ,xlab=lab_nat_log , font.lab=2, main="Comparison of protein abundance between SWAT-TEF2pr-mCherry \nvs SWAT-Nativepr-GFP",  xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5])) 
text(1,11.5,labels=rall, font=2)
# text(10,2.6,labels=rq1, font=2)
# text(10,2.2,labels=rq2, font=2)
# text(9.85,1.8,labels=rq3, font=2)
# text(10,1.4,labels=rq4, font=2)


#plot8: Histogram of Native
h1=hist(log2(data4_s$Native.GfpIntensity.FG.), breaks=100,main="Histogram of SWAT-Native-GFP" ,xlab="SWAT-Native-GFP signal intensity (a.u.)", ylab="Number of strains", col=rgb(0,0,1,1/4), xlim=c(0,12), ylim=c(0,100), xaxt='n')
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)


#plot9: Histogram of Generic
h2=hist(log2(data4_s$Generic.GfpIntensity.FG.) , breaks=100, main="Histogram of SWAT-NOP1pr-GFP",xlab="SWAT-NOP1pr-GFP signal intensity (a.u.)", ylab="Number of strains" , xlim=c(0,12), ylim=c(0,200), col=rgb(0,1,0,1/4), xaxt='n')
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)

#plot10: Histogram of Tef2
h3=hist(log2(data4_s$Tef2pr.mCherryIntensity.FG.) , breaks=100,main="Histogram of SWAT-TEF2pr-mCherry",xlab="SWAT-TEF2pr-mCherry signal intensity (a.u.)", ylab="Number of strains", col=rgb(1,0,0,1/4), xlim=c(0,12), ylim=c(0,140), xaxt='n')
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)

#plot 11: Multiple histogram
plot( h2, col=rgb(0,1,0,1/4), main="Histogram of the three SWAT libraries", xlab="signal intensity (a.u.)", ylab="Number of Strains", xlim=c(0,12), ylim=c(0,200),xaxt='n')  # first histogram
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)
plot( h1, col=rgb(0,0,1,1/4), add=T ) 
plot( h3, col=rgb(1,0,0,1/4), add=T )
legend("topleft",c("Native-GFP","NOP1pr-GFP", "TEF2pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1,1/4),rgb(0,1,0,1/4),rgb(1,0,0,1/4)),bty="n",text.font=1) # 


#plot 12: C-terminal vs Native (log)
C_log=log2(data4_s$C..GFP.intensity)
Nat_log=log2(data4_s$Native.GfpIntensity.FG.) #for negative values, log2 calculation isn't possible. NaN is given. 
#axes_lim=find_log_data_lim(Nat_log,C_log)
plot(Nat_log ,C_log,type="p",pch=19,cex=0.35,ylab=lab_cGFP_log ,xlab=lab_nat_log , font.lab=2, main="Comparison of protein abundance between C' GFP proteins \nvs SWAT-Native-GFP")
#plot(Nat_log ,C_log,type="p",pch=19,cex=0.35,ylab="log2(C' GFP protein abundance (GFP intensity (a.u.)))" ,xlab="log2(SWAT-Native-GFP protein abundance (GFP intensity (a.u.)))" , font.lab=2, main="Signal intensity of C'-tagged GFP proteins vs SWAT-Native-pr-GFP", xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
#testing correllation spearman, based on Native quantile
myRs=Cal_R_native_q('C..GFP.intensity','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,11.5,labels=rall, font=2)
# text(1,11.1,labels=rq1, font=2)
# text(1,10.7,labels=rq2, font=2)
# text(1,10.3,labels=rq3, font=2)
# text(1,9.9,labels=rq4, font=2)


#plot 12b: C-terminal vs Native (log)
AF_lim=20
higher_AF=data4_s[which(data4_s$C..GFP.intensity>AF_lim),]
C_log2=log2(higher_AF$C..GFP.intensity)
Nat_log2=log2(higher_AF$Native.GfpIntensity.FG.) #for negative values, log2 calculation isn't possible. NaN is given.
plot(Nat_log2 ,C_log2,type="p",pch=19,cex=0.35,ylab="C' GFP signal intensity (a.u.)" ,xlab="SWAT-Native-GFP signal intensity (a.u.)" , font.lab=2, main="Comparison of protein abundance between C' GFP proteins \nvs SWAT-Native-GFP", col='black', ylim=c(3.5,12), yaxt='n', xaxt='n')
myTicks=axTicks(2)
labels <- sapply(myTicks,function(i) 2^ i)
axis(2,at=myTicks,labels=labels)
myTicks=axTicks(1)
labels <- sapply(myTicks,function(i) 2^ i)
axis(1,at=myTicks,labels=labels)
Rhigh=cor.test(higher_AF$C..GFP.intensity,higher_AF$Native.GfpIntensity.FG.,method="spearman",use="pairwise.complete.obs")
Rhigh2=round(Rhigh$estimate, 2)
rall=bquote(R ~ " = " ~ .(Rhigh2))  
text(1,11.5,labels=rall, font=2, col='black')

lower_AF=data4_s[which(data4_s$C..GFP.intensity<=AF_lim),]
C_log3=log2(lower_AF$C..GFP.intensity)
Nat_log3=log2(lower_AF$Native.GfpIntensity.FG.) #for negative values, log2 calculation isn't possible. NaN is given.
lines(Nat_log3 ,C_log3,type="p",pch=19,cex=0.35, font.lab=2, col='blue')
Rlow=cor.test(lower_AF$C..GFP.intensity,lower_AF$Native.GfpIntensity.FG.,method="spearman",use="pairwise.complete.obs")
Rlow2=round(Rlow$estimate, 2)
rall=bquote(R ~ " = " ~ .(Rlow2))  
text(1,11,labels=rall, font=2, col='blue')

#what are the genes that were under AF limit in the C terminal, but had value for the native
# BlueGenes=lower_AF[which(!is.na(lower_AF$Native.GfpIntensity.FG.)),c('ORF','Gene','Native.GfpIntensity.FG.','C..GFP.intensity')]
# print('Genes under AF limit in the C terminal, that had value for the native')
# print(BlueGenes)

#plot 13: C-terminal vs Generic (log)
Gen_log=log2(data4_s$Generic.GfpIntensity.FG. ) #for negative values, log2 calculation isn't possible. NaN is given. 
#axes_lim=find_log_data_lim(Gen_log,C_log)
plot(Gen_log ,C_log,type="p",pch=19,cex=0.35,ylab=lab_cGFP_log,xlab=lab_nop_log , font.lab=2, main="Comparison of protein abundance between C' GFP proteins \nvs SWAT-NOP1pr-GFP")
#plot(Gen_log ,C_log,type="p",pch=19,cex=0.35,ylab="log2(C' GFP protein abundance (GFP intensity (a.u.)))" ,xlab="log2(SWAT-NOP1-pr-GFP protein abundance (GFP intensity (a.u.)))" , font.lab=2, main="Signal intensity of C'-tagged GFP proteins vs SWAT-NOP1-pr-GFP",xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('C..GFP.intensity','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,11.5,labels=rall, font=2)
# text(1,11.1,labels=rq1, font=2)
# text(1,10.7,labels=rq2, font=2)
# text(1,10.3,labels=rq3, font=2)
# text(1,9.9,labels=rq4, font=2)

#plot 14: C-terminal vs Tef2 (log)
Tef_log=log2(data4_s$Tef2pr.mCherryIntensity.FG.) #for negative values, log2 calculation isn't possible. NaN is given. 
#axes_lim=find_log_data_lim(Tef_log,C_log)
plot(Tef_log ,C_log,type="p",pch=19,cex=0.35,ylab=lab_cGFP_log ,xlab=lab_tef_log , font.lab=2, main="Comparison of protein abundance between C' GFP proteins \nvs SWAT-TEF2pr-mCherry")
#plot(Tef_log ,C_log,type="p",pch=19,cex=0.35,ylab="C' GFP protein abundance (GFP intensity (a.u.))" ,xlab="SWAT-TEF2-pr-mCherry protein abundance (GFP intensity (a.u.))" , font.lab=2, main="Signal intensity of C'-tagged GFP proteins vs SWAT-TEF2-pr-mCherry",xlim=c(axes_lim[1],axes_lim[2]), ylim=c(axes_lim[1],axes_lim[2]))
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('C..GFP.intensity','Tef2pr.mCherryIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,11.5,labels=rall, font=2)
# text(1,11.1,labels=rq1, font=2)
# text(1,10.7,labels=rq2, font=2)
# text(1,10.3,labels=rq3, font=2)
# text(1,9.9,labels=rq4, font=2)

#plot 12-14 together
plot(Nat_log ,C_log,type="p",pch=20,cex=0.35,ylab=lab_cGFP_log ,xlab=lab_all_3 , font.lab=2, main="Comparison of protein abundance between C' GFP proteins \nvs SWAT proteins", col=rgb(0,0,1))
lines(Gen_log ,C_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,C_log,type="p",pch=20,cex=0.35, font.lab=2,col=rgb(1,0,0))
legend("topleft",c("Native-GFP","NOP1pr-GFP", "TEF2pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 


#plot 15: Ribosome vs Native (log)
Rib_log=log2(data4_s$Ribosome ) 
plot(Nat_log ,Rib_log,type="p",pch=19,cex=0.35,ylab=lab_ribo_log ,xlab=lab_nat_log, font.lab=2, main="Comparison of translation rates as measured by ribosome profiling \nvs SWAT-Native-GFP protein abundance")
#testing correllation spearman, based on Native quantile
myRs=Cal_R_native_q('Ribosome','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,13,labels=rall, font=2)
# text(10,2.2,labels=rq1, font=2)
# text(10,1.4,labels=rq2, font=2)
# text(10,0.6,labels=rq3, font=2)
# text(10,-0.2,labels=rq4, font=2)


#plot 16: Ribosome vs Generic (log)
plot(Gen_log ,Rib_log,type="p",pch=19,cex=0.35,ylab=lab_ribo_log ,xlab=lab_nop_log , font.lab=2, main="Comparison of translation rates as measured by ribosome profiling \nvs SWAT-NOP1pr-GFP protein abundance")
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('Ribosome','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,13,labels=rall, font=2)
# text(1,11.2,labels=rq1, font=2)
# text(1,10.4,labels=rq2, font=2)
# text(1,9.6,labels=rq3, font=2)
# text(1,8.8,labels=rq4, font=2)

#plot 17: Ribosome vs Tef2 (log)
plot(Tef_log ,Rib_log,type="p",pch=19,cex=0.35,ylab=lab_ribo_log ,xlab=lab_tef_log , font.lab=2, main="Comparison of translation rates as measured by ribosome profiling \nvs SWAT-TEF2pr-mCherry protein abundance")
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('Ribosome','Tef2pr.mCherryIntensity.FG.')
# rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,13,labels=rall, font=2)
# text(1,11.2,labels=rq1, font=2)
# text(1,10.4,labels=rq2, font=2)
# text(1,9.6,labels=rq3, font=2)
# text(1,8.8,labels=rq4, font=2)

#plot 15-17 together
plot(Nat_log ,Rib_log,type="p",pch=20,cex=0.35,ylab=lab_ribo_log ,xlab=lab_all_3 , font.lab=2, main="Comparison of translation rates as measured by ribosome profiling \nvs SWAT protein abundance", col=rgb(0,0,1))
lines(Gen_log ,Rib_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,Rib_log,type="p",pch=20,cex=0.35, font.lab=2, col=rgb(1,0,0))
legend("topleft",c("Native-GFP","NOP1pr-GFP", "TEF2pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 

#plot 18: Massspec vs Native (log)
Mass_log=log2(data4_s$Mass.Spectroscopy ) 
plot(Nat_log ,Mass_log,type="p",pch=19,cex=0.35,ylab=lab_mas_log ,xlab=lab_nat_log , font.lab=2, main="Comparison of abundance of proteins as measured by \nmass spectrometry vs SWAT-Native-GFP protein abundance")
#testing correllation spearman, based on native quantile
myRs=Cal_R_native_q('Mass.Spectroscopy','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,40,labels=rall, font=2)
# text(10,23,labels=rq1, font=2)
# text(10,22,labels=rq2, font=2)
# text(10,21,labels=rq3, font=2)
# text(10,20,labels=rq4, font=2)

#plot 19: Massspec vs Generic (log)
plot(Gen_log ,Mass_log,type="p",pch=19,cex=0.35,ylab=lab_mas_log ,xlab=lab_nop_log , font.lab=2, main="Comparison of abundance of proteins as measured by \nmass spectrometry vs SWAT-NOP1pr-GFP protein abundance")
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('Mass.Spectroscopy','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,40,labels=rall, font=2)
# text(1,38,labels=rq1, font=2)
# text(1,37,labels=rq2, font=2)
# text(1,36,labels=rq3, font=2)
# text(1,35,labels=rq4, font=2)


#plot 20: Massspec vs Tef2 (log)
plot(Tef_log ,Mass_log,type="p",pch=19,cex=0.35,ylab=lab_mas_log ,xlab=lab_tef_log , font.lab=2, main="Comparison of abundance of proteins as measured by \nmass spectrometry vs SWAT-TEF2pr-mCherry protein abundance")
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('Mass.Spectroscopy','Tef2pr.mCherryIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,40,labels=rall, font=2)
# text(9,20,labels=rq1, font=2)
# text(9,19,labels=rq2, font=2)
# text(9,18,labels=rq3, font=2)
# text(9,17,labels=rq4, font=2)

#plot 18-20 together
plot(Nat_log ,Mass_log,type="p",pch=20,cex=0.35,ylab=lab_mas_log ,xlab=lab_all_3 , font.lab=2, main="Comparison of abundance of proteins as measured by \nmass spectrometry vs SWAT protein abundance", col=rgb(0,0,1))
lines(Gen_log ,Mass_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,Mass_log,type="p",pch=20,cex=0.35, font.lab=2,col=rgb(1,0,0))
legend("topleft",c("Native-pr-GFP","NOP1-pr-GFP", "TEF2-pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 


#plot 21: mRNA vs Native (log)
mRNA_log=log2(data4_s$mRNA ) 
plot(Nat_log ,mRNA_log,type="p",pch=19,cex=0.35,ylab=lab_mRNA_log ,xlab=lab_nat_log , font.lab=2, main="Comparison of mRNA abundance \nvs SWAT-Native-GFP protein abundance ")
#testing correllation spearman, based on native quantile
myRs=Cal_R_native_q('mRNA','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,14,labels=rall, font=2)
# text(-6,12.8,labels=rq1, font=2)
# text(-6,12.1,labels=rq2, font=2)
# text(-6,11.4,labels=rq3, font=2)
# text(-6,10.7,labels=rq4, font=2)

#plot 22: mRNA vs Generic (log)
plot(Gen_log ,mRNA_log,type="p",pch=19,cex=0.35,ylab=lab_mRNA_log ,xlab=lab_nop_log , font.lab=2, main="Comparison of mRNA abundance \nvs SWAT-NOP1pr-GFP protein abundance ")
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('mRNA','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,14,labels=rall, font=2)
# text(-2.6,12.8,labels=rq1, font=2)
# text(-2.6,12.1,labels=rq2, font=2)
# text(-2.6,11.4,labels=rq3, font=2)
# text(-2.78,10.7,labels=rq4, font=2)


#plot 23: mRNA vs Tef2 (log)
plot(Tef_log ,mRNA_log,type="p",pch=19,cex=0.35,ylab=lab_mRNA_log ,xlab=lab_tef_log , font.lab=2, main="Comparison of mRNA abundance \nvs SWAT-TEF2pr-mCherry protein abundance ")
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('mRNA','Tef2pr.mCherryIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,14,labels=rall, font=2)
# text(-9,12.8,labels=rq1, font=2)
# text(-9,12.1,labels=rq2, font=2)
# text(-9,11.4,labels=rq3, font=2)
# text(-9,10.7,labels=rq4, font=2)

#plot 21-23 together
plot(Nat_log ,mRNA_log,type="p",pch=20,cex=0.35,ylab=lab_mRNA_log ,xlab=lab_all_3 , font.lab=2, main="Comparison of mRNA abundance \nvs SWAT protein abundance", col=rgb(0,0,1))
lines(Gen_log ,mRNA_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,mRNA_log,type="p",pch=20,cex=0.35, font.lab=2,col=rgb(1,0,0))
legend("topleft",c("Native-pr-GFP","NOP1-pr-GFP", "TEF2-pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 

#plot 24: FACS vs Native (log)
FACS_log=log2(data4_s$FACS ) 
plot(Nat_log ,FACS_log,type="p",pch=19,cex=0.35,ylab=lab_FACS_log ,xlab=lab_nat_log , font.lab=2, main="Abundance of proteins as measured by flow cytometry analysis \nof C'-tagged GFP proteins")
#testing correllation spearman, based on native quantile
myRs=Cal_R_native_q('FACS','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,15.5,labels=rall, font=2)
# text(-5.8,15,labels=rq1, font=2)
# text(-5.8,14.5,labels=rq2, font=2)
# text(-5.8,14,labels=rq3, font=2)
# text(-6,13.5,labels=rq4, font=2)

#plot 25: FACS vs Generic (log)
plot(Gen_log ,FACS_log,type="p",pch=19,cex=0.35,ylab=lab_FACS_log ,xlab=lab_nop_log , font.lab=2, main="Abundance of proteins as measured by flow cytometry analysis \nof C'-tagged GFP protein")
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('FACS','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,15.5,labels=rall, font=2)
# text(-2.6,15,labels=rq1, font=2)
# text(-2.6,14.5,labels=rq2, font=2)
# text(-2.6,14,labels=rq3, font=2)
# text(-2.6,13.5,labels=rq4, font=2)


#plot 26: FACS vs Tef2 (log)
plot(Tef_log ,FACS_log,type="p",pch=19,cex=0.35,ylab=lab_FACS_log ,xlab=lab_tef_log , font.lab=2, main="Abundance of proteins as measured by flow cytometry analysis of \nC'-tagged GFP protein")
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('FACS','Tef2pr.mCherryIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,15.5,labels=rall, font=2)
# text(-9,15,labels=rq1, font=2)
# text(-9,14.5,labels=rq2, font=2)
# text(-9,14,labels=rq3, font=2)
# text(-9,13.5,labels=rq4, font=2)

#plot 24-26 together
plot(Nat_log ,FACS_log,type="p",pch=20,cex=0.35,ylab=lab_FACS_log ,xlab=lab_all_3 , font.lab=2, main="Abundance of proteins as measured by flow cytometry analysis \nof C'-tagged GFP protein", col=rgb(0,0,1))
lines(Gen_log ,FACS_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,FACS_log,type="p",pch=20,cex=0.35, font.lab=2,col=rgb(1,0,0))
legend("topleft",c("Native-pr-GFP","NOP1-pr-GFP", "TEF2-pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 


#plot 27: Half life vs Native (log)
Half_log=log2(data4_s$Raw.Half.life ) 
plot(Nat_log ,Half_log,type="p",pch=19,cex=0.35,ylab=lab_half_log ,xlab=lab_nat_log , font.lab=2, main="Comparison of protein half-life \nvs SWAT-Native-GFP protein abundance")
#testing correllation spearman, based on native quantile
myRs=Cal_R_native_q('Raw.Half.life','Native.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,16,labels=rall, font=2)
# text(-6.8,14.8,labels=rq1, font=2)
# text(-7,14.1,labels=rq2, font=2)
# text(-7,13.4,labels=rq3, font=2)
# text(-7,12.8,labels=rq4, font=2)

#plot 28: Half life vs Generic (log)
plot(Gen_log ,Half_log,type="p",pch=19,cex=0.35,ylab=lab_half_log ,xlab=lab_nop_log , font.lab=2, main="Comparison of protein half-life \nvs SWAT-NOP1-GFP protein abundance")
#testing correllation spearman, based on Generic quantile
myRs=Cal_R_generic_q('Raw.Half.life','Generic.GfpIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,16,labels=rall, font=2)
# text(-2.6,14.8,labels=rq1, font=2)
# text(-2.6,14.1,labels=rq2, font=2)
# text(-2.6,13.4,labels=rq3, font=2)
# text(-2.6,12.8,labels=rq4, font=2)


#plot 29: Half life vs Tef2 (log)
plot(Tef_log ,Half_log,type="p",pch=19,cex=0.35,ylab=lab_half_log ,xlab=lab_tef_log , font.lab=2, main="Comparison of protein half-life \nvs SWAT-TEF2-mCherry protein abundance")
#testing correllation spearman, based on Tef2 quantile
myRs=Cal_R_tef2_q('Raw.Half.life','Tef2pr.mCherryIntensity.FG.')
rall=bquote(R ~ " = " ~ .(myRs[1]))  
# rq1=bquote(R[Q1] ~ " = " ~ .(myRs[2]))    
# rq2=bquote(R[Q2] ~ " = " ~ .(myRs[3]))    
# rq3=bquote(R[Q3] ~ " = " ~ .(myRs[4]))    
# rq4=bquote(R[Q4] ~ " = " ~ .(myRs[5]))    
text(1,16,labels=rall, font=2)
# text(-9,14.8,labels=rq1, font=2)
# text(-9,14.1,labels=rq2, font=2)
# text(-9,13.4,labels=rq3, font=2)
# text(-9.2,12.8,labels=rq4, font=2)

#plot 24-26 together
plot(Nat_log ,Half_log,type="p",pch=20,cex=0.35,ylab=lab_half_log ,xlab=lab_all_3 , font.lab=2, main="Comparison of protein half-life \nvs SWAT protein abundance", col=rgb(0,0,1))
lines(Gen_log ,Half_log,type="p",pch=20,cex=0.35 , font.lab=2, col=rgb(0,1,0))
lines(Tef_log ,Half_log,type="p",pch=20,cex=0.35, font.lab=2,col=rgb(1,0,0))
legend("topleft",c("Native-pr-GFP","NOP1-pr-GFP", "TEF2-pr-mCherry"),cex=1.1,pt.cex=c(1.1,1.1),fill=c(rgb(0,0,1),rgb(0,1,0),rgb(1,0,0)),bty="n",text.font=1) # 

# OutputFile='Correlation_map'
# e<-paste(myDir,"/",OutputFile,sep="")
# 
# graphics.off()  
# pdf(e,width=7,height=5.8,family='Helvetica', paper='A4')

#Plot 27 - correlation plot
library(corrplot)
PlotData=data4_s[,c(10,8,7,9,15,16,12,13,14)]
colnames(PlotData)=c("C' GFP","Native-GFP","NOP1pr-GFP","TEF2pr-mCherry","mRNA abundance","mRNA Half Life","Translation rate","Protein abundance","Protein half-life")
M = cor(PlotData,PlotData,method="spearman",use="pairwise.complete.obs")
#I manually change the values of the principal diagonal, to allow changing the range of the color scheme from 0.07 to 0.87
diag(M)=rep(0.5,each=9)
corrplot(M, type='upper',method="color", cl.lim = c(0.07,0.87), tl.col=1,  addCoef.col = "black", diag=FALSE,is.corr = FALSE, tl.srt=45)

PlotData2=PlotData[,c(2:9)]
M2 = cor(PlotData2,PlotData2,method="spearman",use="pairwise.complete.obs")
#I manually change the values of the principal diagonal, to allow changing the range of the color scheme from 0.07 to 0.87
diag(M2)=rep(0.5,each=)
corrplot(M2, type='upper',method="color", cl.lim = c(0.15,0.87), tl.col=1,  addCoef.col = "black", diag=FALSE,is.corr = FALSE, tl.srt=45)

dev.off()

# 
# labels=c(2^3,2^4,2^5,2^6,2^7,2^8,2^9)
# xat<-c(3,4,5,6,7,8,9)
# 
# axis(1,at=xat,lwd=1,labels,cex.axis=0.65)












