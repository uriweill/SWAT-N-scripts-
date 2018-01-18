### Peroxisomal SWAT DHFR PCA screen 

#### PART 1: GET THE PROTEIN IDs, PLATE ORGANIZATION AND IMAGE READING ORDER TO CREATE
####  A MATRIX WITH ALL PROTEIN PAIRS IN THE READING ORDER
#### PART 2: COMPLETE THE TABLE WITH MEASUREMENTS AND REMOVE FREQUENT FLYERS
#### PART 3: NORMALIZE WITH DAY 0 VALUES AND ROW-COLUMN MEDIAN VALUES
#### PART 4: INCLUDE EXTERNAL INTERACTION DATA
#### PART 5: CALCULATE Z SCORES FOR RAW VALUES AND NORMALIZED VALUES (AFTER DAY 0 NORMALIZATION
####  AND AFTER DAY 0 + ROW-COLUMN NORMALIZATION) AND ADDRESS RECALL SCORES
#### PART 6: CORRECT FOR PROTEIN ABUNDANCE 
#### PART 7: CREATION OF HEAT MAP FOR ARTICLE


#### PART 1: GET THE PROTEIN IDs, PLATE ORGANIZATION AND IMAGE READING ORDER TO CREATE
#### A MATRIX WITH ALL PROTEIN PAIRS IN THE READING ORDER

PATH = "C:/Users/Bram/Documents/Project SWAT DHFR PCA/SWAT general info and data/"

# Peroxisomal protein IDs
Peroxi = read.csv(paste(PATH,"Peroxi-csv.csv",sep=""),sep=";",header=F)

# Bait plate organization (3 unique plates)
Baits = read.csv(paste(PATH,"Peroxi-baits.csv",sep=""),sep=";",header=F)
baits = list(Baits[1:32,],Baits[33:64,],Baits[65:96,]) # Split baits (DHFR F[1,2]) into 3 plates

# Prey plate organization (2 unique plates)
Preys = read.csv(paste(PATH,"Peroxi-preys.csv",sep=""),sep=";",header=F)
preys = list(Preys[1:32,],Preys[33:64,]) # Split preys (DHFR F[3]) into 2 plates

Peroxi[,2] = as.numeric(as.character(Peroxi[,2])) # turn Peroxisomal protein IDs in numeric values

## Create row and column coordinates in the order they are read during image analysis
# The order of measuring colony growth in ImageJ is very special:
# 1  2  X  X  5  6  X  X  9  10 ... 45 46
# 3  4  X  X  7  8  X  X  11 12     47 48
# X  X  X  X  X  X  X  X  X  X
# X  X  X  X  X  X  X  X  X  X
# 49 50 X  X 53 54  X  X  57 58 ... 93 94
# 51 52 X  X 55 56  X  X  59 60     95 96
# After 384 measurements, the following 384 are taken as followed:
# 1  2  385 386  5  6  389  390  9  10 ...
# 3  4  387 388  7  8  391  392  11 12
# X  X  X   X    X  X  X    X    X  X
# X  X  X   X    X  X  X    X    X  X
# 49 50 433 434 53 54  437  438  57 58 ...
# 51 52 435 436 55 56  439  440  59 60
# Then the next 384:
# 1   2   385 386  5  6  389  390    9   10 
# 3   4   387 388  7  8  391  392    11  12
# 769 770  X   X  773 774  X    X   777 778
# 771 772  X   X  775 776  X    X   779 780 ...
# 49  50  433 434 53  54  437  438  57  58 
# 51  52  435 436 55  56  439  440  59  60
# The final 384 follows the same principles...

# The following lines make sure the row-column coordinates can be obtained from the image analysis order:
roworder = c()
colorder = c()
for (i in 0:1) {
  for (j in 0:1) {
    for (k in 0:7) {
      for (l in 0:11) {
        for (m in 1:2) {
          for (n in 1:2) {
            roworder = c(roworder,m+k*4+i*2)
            colorder = c(colorder,n+l*4+j*2)
          }
        }
      }
    }
  }
} 

# Now we combine the image analysis order with the IDs of the proteins in plate format (the baits and preys lists) 
# to get the list of protein pairs in order of reading.
baits1ordered = c()
baits2ordered = c()
baits3ordered = c()
preys1ordered = c()
preys2ordered = c()
# The following lines match the raw data value analyzed through image analysis with a corresponding protein ID:
for (i in 1:1536) {
  if (is.element(baits[[1]][roworder[i],colorder[i]],Peroxi[,2])) {
    baits1ordered = c(baits1ordered,which(Peroxi[,2]==baits[[1]][roworder[i],colorder[i]]))
  } else { baits1ordered = c(baits1ordered,NA) }
  
  if (is.element(baits[[2]][roworder[i],colorder[i]],Peroxi[,2])) {
    baits2ordered = c(baits2ordered,which(Peroxi[,2]==baits[[2]][roworder[i],colorder[i]]))
  } else { baits2ordered = c(baits2ordered,NA) }
  
  if (is.element(baits[[3]][roworder[i],colorder[i]],Peroxi[,2])) {
    baits3ordered = c(baits3ordered,which(Peroxi[,2]==baits[[3]][roworder[i],colorder[i]]))
  } else { baits3ordered = c(baits3ordered,NA) }
  
  if (is.element(preys[[1]][roworder[i],colorder[i]],Peroxi[,2])) {
    preys1ordered = c(preys1ordered,which(Peroxi[,2]==preys[[1]][roworder[i],colorder[i]]))
  } else { preys1ordered = c(preys1ordered,NA) }
  
  if (is.element(preys[[2]][roworder[i],colorder[i]],Peroxi[,2])) {
    preys2ordered = c(preys2ordered,which(Peroxi[,2]==preys[[2]][roworder[i],colorder[i]]))
  } else { preys2ordered = c(preys2ordered,NA) }
}

# Now we prepare the dataframe used during the remainder of the analysis
# The order at which the plates were analyzed is as followed (with each condition as duplicates):
# C-term DHFR F[1,2] plate 1 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 1 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 1 vs. N-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 1 vs. N-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 2 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 2 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 2 vs. N-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 2 vs. N-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 3 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 3 vs. C-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 3 vs. N-term DHFR F[3] plate 1
# C-term DHFR F[1,2] plate 3 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 1 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 1 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 1 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 1 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 2 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 2 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 2 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 2 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 3 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 3 vs. C-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 3 vs. N-term DHFR F[3] plate 1
# N-term DHFR F[1,2] plate 3 vs. N-term DHFR F[3] plate 1
BaitID = rep(c(rep(baits1ordered,4),rep(baits2ordered,4),rep(baits3ordered,4)),2) # Provide Bait IDs in image reading order
BaitORF = as.character(Peroxi[BaitID,1]) # Provide Bait ORF names
Bait = as.character(Peroxi[BaitID,3]) # Provide Bait protein names
PreyID = rep(c(preys1ordered,preys2ordered),12) # Provide Prey IDs in image reading order
PreyORF = as.character(Peroxi[PreyID,1]) # Provide Prey ORF names
Prey = as.character(Peroxi[PreyID,3]) # Provide Prey protein names
Pl = 24 # Unique plates
Si = 1536 # Plate size
Re = 2 # Replicates
Plate = rep(1:Pl, each = Si)
Interaction_ID = 1:(Si*Pl) # Unique ID for each tested interaction
base_name = "PEROXI_1-24_Mtx1_"
sourcedir = "C:/Users/Bram/Documents/Project SWAT DHFR PCA/SWAT screens/SWAT peroxi full screen - images and raw data/"
days = c("day0","day4") # Days of measurement

PexPCA = data.frame("ID" = Interaction_ID, "Plate" = Plate, 
                    "BaitID" = BaitID, "BaitORF" = BaitORF,
                    "Bait" = Bait, "Baittag" = c(rep("C",1536*12),rep("N",1536*12)),
                    "PreyID" = PreyID, "PreyORF" = PreyORF,
                    "Prey" = Prey, "Preytag" = rep(c(rep("C",1536*2),rep("N",1536*2)),6),
                    "Baitrow" = rep(roworder,24), "Preycol" = rep(colorder,24),
                    "Replicate1" = rep(NA,1536*24),
                    "Replicate2" = rep(NA,1536*24)
                    )

# Change the type of data to numeric or character:
PexPCA$BaitID = as.numeric(as.character(PexPCA$BaitID))
PexPCA$PreyID = as.numeric(as.character(PexPCA$PreyID))
PexPCA$BaitORF = as.character(PexPCA$BaitORF)
PexPCA$PreyORF = as.character(PexPCA$PreyORF)
PexPCA$Bait = as.character(PexPCA$Bait)
PexPCA$Prey = as.character(PexPCA$Prey)
Peroxi[,1] = as.character(Peroxi[,1])
Peroxi[,3] = as.character(Peroxi[,3])

#### PART 2: COMPLETE THE TABLE WITH MEASUREMENTS AND REMOVE FREQUENT FLYERS

# Read raw data
day = "day4"
for (plates in 1:(Re*Pl)) {
    tmp = read.table(paste(sourcedir,day,"_Peroxi_Full_MTX1_xls/",
                           base_name,"d",strsplit(day,split="")[[1]][nchar(day)],"_",sprintf("%03d",plates),".xls",sep=""),sep="\t",header=F)
    if ((plates %% 2) == 0) {
      PexPCA$Replicate1[((ceiling(plates/2)-1)*1536+1):((ceiling(plates/2)-1)*1536+1536)] = as.numeric(tmp[,3])
    } else {
      PexPCA$Replicate2[((ceiling(plates/2)-1)*1536+1):((ceiling(plates/2)-1)*1536+1536)] = as.numeric(tmp[,3])
    }
}
# Convert values to log2 scale
PexPCA$Replicate1 = log2(PexPCA$Replicate1+1)
PexPCA$Replicate2 = log2(PexPCA$Replicate2+1)

toRemove = c()
for (i in 1:nrow(PexPCA)) {
  if (is.na(PexPCA$BaitID[i]) | is.na(PexPCA$PreyID[i])) {
    toRemove = c(toRemove,i)
  }
}
PexPCA = PexPCA[-toRemove,]

genes = unique(PexPCA$Bait)
S12 = c(); S3 = c(); T12 = c(); T3 = c()
for (i in 1:length(genes)) {
  tmp = PexPCA[which(PexPCA$Bait == genes[i]),]
  tmp2 = tmp[tmp$Baittag == "N",]
  S12 = c(S12, median(tmp2$Replicate1,na.rm=T))
  tmp2 = tmp[tmp$Baittag == "C",]
  T12 = c(T12, median(tmp2$Replicate1,na.rm=T))
  tmp = PexPCA[which(PexPCA$Prey == genes[i]),]
  tmp2 = tmp[tmp$Preytag == "N",]
  S3 = c(S3, median(tmp2$Replicate1,na.rm=T))
  tmp2 = tmp[tmp$Preytag == "C",]
  T3 = c(T3, median(tmp2$Replicate1,na.rm=T))
}
# Obtain day 0 colony size for later normalizaton
PexPCA0 = PexPCA
day = "day0"
for (plates in 1:(Re*Pl)) {
  tmp = read.table(paste(sourcedir,day,"_Peroxi_Full_MTX1_xls/",
                         base_name,"d",strsplit(day,split="")[[1]][nchar(day)],"_",sprintf("%03d",plates),".xls",sep=""),sep="\t",header=F)
  if ((plates %% 2) == 0) {
    PexPCA0$Replicate1[((ceiling(plates/2)-1)*1536+1):((ceiling(plates/2)-1)*1536+1536)] = as.numeric(tmp[,3])
  } else {
    PexPCA0$Replicate2[((ceiling(plates/2)-1)*1536+1):((ceiling(plates/2)-1)*1536+1536)] = as.numeric(tmp[,3])
  }
}
PexPCA0$Replicate1 = log2(PexPCA0$Replicate1+1)
PexPCA0$Replicate2 = log2(PexPCA0$Replicate2+1)



# Remove values below 13.5 (this is a shortcut way to remove positions without actual colonies)
PexPCA02 = PexPCA0[PexPCA$Replicate1 > 13.5 & PexPCA$Replicate2 > 13.5,]
PexPCA2 = PexPCA[PexPCA$Replicate1 > 13.5 & PexPCA$Replicate2 > 13.5,]

# For good measure, remove values if there is not supposed to be a value
PexPCA02 = PexPCA02[!is.na(PexPCA2$BaitID) & !is.na(PexPCA2$PreyID),]
PexPCA2 = PexPCA2[!is.na(PexPCA2$BaitID) & !is.na(PexPCA2$PreyID),]

# Remove frequent flyers (sticky highly abundant proteins that generally interact non-specifically with
# many other proteins):
FF = read.csv(paste(PATH,"DHFR PCA frequent flyers.csv",sep=""),sep=";",header=F)
FFB = rep(1,nrow(PexPCA2)); FFP = FFB
for (i in 1:nrow(PexPCA2)) {
  if (is.element(as.character(PexPCA2$Bait[i]),as.character(FF[,1]))) {
    FFB[i] = 0
  }
  if (is.element(as.character(PexPCA2$Prey[i]),as.character(FF[,1]))) {
    FFP[i] = 0
  }
}
PexPCA2 = PexPCA2[FFB == 1 & FFP == 1,]
PexPCA02 = PexPCA02[FFB == 1 & FFP == 1,]

### PART 3: NORMALIZE WITH DAY 0 VALUES AND ROW-COLUMN MEDIAN VALUES

day40 = data.frame("Colony_size_day0" = PexPCA02$Replicate1,
                   "Colony_size_day4" = PexPCA2$Replicate1)

# Linear regression of day 4 data against day 0 data
library(affy)
m = lm(PexPCA2$Replicate1~PexPCA02$Replicate1)
intercept = m[[1]][1]; slope = m[[1]][2]

library(ggplot2)
ggplot(subset(day40,Colony_size_day0 > 10.5), aes(x=Colony_size_day0, y = Colony_size_day4)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope=slope,intercept=intercept,lwd=2,col="blue") +
  xlim(c(12,16)) +
  ylim(c(12,19)) +
  xlab("Colony size (day 0)") +
  ylab("Colony size (day 4)") +
  theme_bw()

# Normalization step (we keep the original raw values separately):
PexPCA2$Repl1orig = PexPCA2$Replicate1
PexPCA2$Repl2orig = PexPCA2$Replicate2
PexPCA2$Replicate1 = PexPCA2$Replicate1 - PexPCA02$Replicate1*slope - intercept
m = lm(PexPCA2$Replicate2~PexPCA02$Replicate2)
intercept = m[[1]][1]; slope = m[[1]][2]
PexPCA2$Replicate2 = PexPCA2$Replicate2 - PexPCA02$Replicate2*slope - intercept

# Normalization of the data according to the median row and column values:
PexPCA2 = PexPCA2[order(PexPCA2$ID),]
norm1Tot = c()
norm2Tot = c()
for (plID in 1:Pl) {
  pltmp = PexPCA2[PexPCA2$Plate == plID,]
  colID = unique(pltmp$Preycol)
  rowID = unique(pltmp$Baitrow)
  coltmp1 = c()
  coltmp2 = c()
  for (i in 1:length(colID)) {
    coltmp1 = c(coltmp1,median(pltmp$Replicate1[pltmp$Preycol == colID[i]],na.rm=T))
    coltmp2 = c(coltmp2,median(pltmp$Replicate2[pltmp$Preycol == colID[i]],na.rm=T))
  }
  rowtmp1 = c()
  rowtmp2 = c()
  for (i in 1:length(rowID)) {
    rowtmp1 = c(rowtmp1,median(pltmp$Replicate1[pltmp$Baitrow == rowID[i]],na.rm=T))
    rowtmp2 = c(rowtmp2,median(pltmp$Replicate2[pltmp$Baitrow == rowID[i]],na.rm=T))
  }
  rowref1 = rep(0,32); colref1 = rep(0,48); rowref1[rowID] = rowtmp1; colref1[colID] = coltmp1
  rowref2 = rep(0,32); colref2 = rep(0,48); rowref2[rowID] = rowtmp2; colref2[colID] = coltmp2
  norm1 = c()
  norm2 = c()
  for (i in 1:nrow(pltmp)) {
    norm1 = c(norm1,pltmp$Replicate1[i] - mean(rowref1[pltmp$Baitrow[i]],colref1[pltmp$Preycol[i]]))
    norm2 = c(norm2,pltmp$Replicate2[i] - mean(rowref2[pltmp$Baitrow[i]],colref2[pltmp$Preycol[i]]))
  }
  norm1Tot = c(norm1Tot,norm1)
  norm2Tot = c(norm2Tot,norm2)
}

PexPCA2$RepNorm1 = norm1Tot
PexPCA2$RepNorm2 = norm2Tot

## Shapiro test for normal distributions:
NormRepl1 = c(); NormRepl2 = c()
for (i in 1:Screen_size) {
  NormRepl1 = c(NormRepl1,shapiro.test(PexPCA2$RepNorm1[PexPCA2$Plate == i])[[2]])
  NormRepl2 = c(NormRepl2,shapiro.test(PexPCA2$RepNorm2[PexPCA2$Plate == i])[[2]])
}  

### PART 4: INCLUDE EXTERNAL INTERACTION DATA

# Get all physical protein-protein interaction data of yeast and exclude the data
# from our own screen (Tarassov 2008)
IntData = read.csv(paste(PATH,"Yeast-BioGRID-4-2017.txt",sep=""),sep="\t",header=T)
IntData = data.frame(IntData)
IntData = IntData[as.character(IntData[,6])=="physical",]
Tarassov2 = which(as.character(IntData[,7])=="Tarassov K (2008)")
Tarassov = IntData[Tarassov2,]
IntData = IntData[-Tarassov2,]
PexPCA2$Tarassov = rep(0,nrow(PexPCA2))
PexPCA2$Int_in_literature = rep(0,nrow(PexPCA2))

# We narrow down the data to interactions tested in this experiment:
IntDataPeroxi = IntData
PeroxiPPI = c()
for (i in 1:nrow(IntDataPeroxi)) {
  if (is.element(IntData[i,1],Peroxi[,1]) & is.element(IntData[i,2],Peroxi[,1])) {
    PeroxiPPI = c(PeroxiPPI,i)
    print(i)
  }
}
IntDataPeroxi = IntDataPeroxi[PeroxiPPI,]

# We label interactons that were identified in the literature in our dataframe:
for (i in 1:nrow(PexPCA2)) {
  tmp1 = intersect(c(which(IntDataPeroxi[,1]==as.character(PexPCA2$BaitORF[i]))),
                   c(which(IntDataPeroxi[,2]==as.character(PexPCA2$PreyORF[i]))))
  tmp2 = intersect(c(which(IntDataPeroxi[,2]==as.character(PexPCA2$BaitORF[i]))),
                   c(which(IntDataPeroxi[,1]==as.character(PexPCA2$PreyORF[i]))))
  tmp3 = intersect(c(which(Tarassov[,1]==as.character(PexPCA2$BaitORF[i]))),
                   c(which(Tarassov[,2]==as.character(PexPCA2$PreyORF[i]))))
  tmp4 = intersect(c(which(Tarassov[,2]==as.character(PexPCA2$BaitORF[i]))),
                   c(which(Tarassov[,1]==as.character(PexPCA2$PreyORF[i]))))
  if ((length(tmp1)+length(tmp2)) > 0) {PexPCA2$Int_in_literature[i] = length(tmp1)+length(tmp2); print(paste("General",i))}
  if ((length(tmp3)+length(tmp4)) > 0) {PexPCA2$Tarassov[i] = length(tmp3)+length(tmp4); print(paste("Tarassov",i))}
}

### PART 5: CALCULATE Z SCORES FOR RAW VALUES AND NORMALIZED VALUES (AFTER DAY 0 NORMALIZATION
### AND AFTER DAY 0 + ROW-COLUMN NORMALIZATION) AND ADDRESS RECALL SCORES
PexPCA2$Zscore1 = rep(NA,nrow(PexPCA2))
PexPCA2$Zscore2 = rep(NA,nrow(PexPCA2))
PexPCA2$Zscore1orig = rep(NA,nrow(PexPCA2))
PexPCA2$Zscore2orig = rep(NA,nrow(PexPCA2))
PexPCA2$Zscore1N = rep(NA,nrow(PexPCA2))
PexPCA2$Zscore2N = rep(NA,nrow(PexPCA2))
PexPCA2$ZscoreMin = rep(NA,nrow(PexPCA2))
PexPCA2$ZscoreMinorig = rep(NA,nrow(PexPCA2))
PexPCA2$ZscoreMinN = rep(NA,nrow(PexPCA2))

for (i in 1:Pl) {
  for (j in 13:18) {
    tmp = PexPCA2[,j][PexPCA2$Plate == i]
    tmp = (tmp-mean(tmp,na.rm=T))/sd(tmp,na.rm=T)
    PexPCA2[,j+8][PexPCA2$Plate == i] = tmp
  }
}

PexPCA2$ZscoreMin = apply(PexPCA2[,c(21,22)],1,min)
PexPCA2$ZscoreMinorig = apply(PexPCA2[,c(23,24)],1,min)
PexPCA2$ZscoreMinN = apply(PexPCA2[,c(25,26)],1,min)

# Focus on recall scores to observe improvements in curve:
Recall_sel1 = rep(list(vector(length=nrow(PexPCA2))),3)

count = 1  
for (i in c(27,28,29)) {
  PexPCA2 = PexPCA2[order(PexPCA2[,i],decreasing=T),]
  tmp = PexPCA2$Int_in_literature
  tmp[tmp!=0] = 1
  for (j in 1:length(tmp)) {
    Recall_sel1[[count]][[j]] = sum(tmp[1:j])/sum(tmp)
  }
  count = count + 1
}
names(Recall_sel1) = c("Day 0 normalization","Raw data","Day 0 and column-row normalization")

ReValues = c()
ID = c()
for (i in 1:3) {
  ReValues = c(ReValues,Recall_sel1[[i]])
}
Recalldf = data.frame("rank"= rep(1:(length(ReValues)/3),3), "Recall" =  ReValues,
                      "ID" = rep(c("Day 0 normalizaton","Raw data","Day 0 and column-row normalization"),
                                   each = length(ReValues)/3))

colchoice = c("dark blue","blue","light blue")
q = ggplot(Recalldf[Recalldf$rank<max(Recalldf$rank),],aes(x=rank,y=Recall,col=ID)) +
  geom_line(size=1,lwd=1) +
  #geom_abline(slope=3/nrow(Recalldf),intercept = 0,lwd=1) +
  scale_color_manual(values = colchoice) +
  #xlim(c(0,500)) +
  #ylim(c(0,0.2)) +
  theme_bw() +
  xlab("Rank")
ggsave("C:/Users/Bram/Documents/Project SWAT DHFR PCA/Recall.png",width = 8, height = 5,q)

write.csv(Recalldf, "C:/Users/Bram/Documents/Project SWAT DHFR PCA/Recall.csv")

### PART 6: CORRECT FOR PROTEIN ABUNDANCE 

# Get the general protein properties
PeroxiProp = read.csv(paste(PATH,"Peroxi properties.csv",sep=""),sep=";",header=T,stringsAsFactors=F)  
PeroxiProp2 = PeroxiProp[,c(1,2,7,10,11,13,18)]
colnames(PeroxiProp2) = c("ORF","Gene","Localization","AbundanceN",
                          "AbundanceC","TMDs","Topology")
PeroxiProp2$AbundanceN = as.numeric(as.character(PeroxiProp2$AbundanceN))
PeroxiProp2$AbundanceC = as.numeric(as.character(PeroxiProp2$AbundanceC))
PeroxiProp2$TMDs = as.numeric(as.character(PeroxiProp2$TMDs))
PeroxiProp2$Topology = as.character(PeroxiProp2$Topology)

# Insert the abundance values of each protein in the data frame
for (i in 1:nrow(PeroxiProp2)) {
  if (is.na(PeroxiProp2$AbundanceC[i] & !is.na(PeroxiProp2$AbundanceN[i]))) {
    PeroxiProp2$AbundanceC[i] = PeroxiProp2$AbundanceN[i]
    print(i)
  }
}

PexPCA2$MinAb = rep(NA,nrow(PexPCA2))
PexPCA2$MaxAb = rep(NA,nrow(PexPCA2))
PexPCA2$MeanAb = rep(NA,nrow(PexPCA2))

PrOrd = c()
uBaitID = unique(PexPCA2$BaitID)
uBaitORF = unique(PexPCA2$BaitORF)
for (i in 1:length(uBaitORF)) {
  PrOrd = c(PrOrd,which(as.character(PeroxiProp2$ORF) == as.character(uBaitORF[i])))
}
freq = rep(0,96)
PrOrd = PrOrd[order(uBaitID)]
uBaitID = uBaitID[order(uBaitID)]
freq[uBaitID] <- PrOrd
PrOrd2 = freq

for (i in 1:nrow(PexPCA2)) {
  if (as.character(PexPCA2$Baittag[i]) == "N") {
    B = PeroxiProp2$AbundanceN[PrOrd2[PexPCA2$BaitID[i]]]
  } else {
    B = PeroxiProp2$AbundanceC[PrOrd2[PexPCA2$BaitID[i]]]
  }
  if (as.character(PexPCA2$Preytag[i]) == "N") {
    P = PeroxiProp2$AbundanceN[PrOrd2[PexPCA2$PreyID[i]]]
  } else {
    P = PeroxiProp2$AbundanceC[PrOrd2[PexPCA2$PreyID[i]]]
  }
  PexPCA2$MinAb[i] = min(B,P)
  PexPCA2$MeanAb[i] = mean(B,P)
  PexPCA2$MaxAb[i] = max(B,P)
}

Replmean = c()
for (i in 1:nrow(PexPCA2)) { Replmean = c(Replmean,mean(PexPCA2$Replicate1[i],PexPCA2$Replicate2[i])) }
PexPCA2$Replmean = Replmean

# Focus on positive results (mininimum Zscore of normalized data > 3)
PexPCA3 = PexPCA2[PexPCA2$ZscoreMinN > 3,]
PexPCA3$IntPos = PexPCA3$Int_in_literature
PexPCA3$IntPos[PexPCA3$Int_in_literature > 0] = 1
PexPCA3$IntPos = PexPCA3$IntPos + 1
PexPCA3$In_literature = "Yes"
PexPCA3$In_literature[PexPCA3$IntPos == 1] = "No"

# Determine the Z score and abundance cutoffs for removal
# of results with high maximum protein abundance and low Z scores
ratio = c()
size = c()
Abopt = c()
Zopt = c()
for (i in seq(30,100,by=3)) {
  for (j in seq(3,5,by=0.05)) {
    tmp = PexPCA3[PexPCA3$MaxAb > i & PexPCA3$ZscoreMinN < j,]
    ratio = c(ratio,sum(tmp$IntPos-1)/nrow(tmp))
    size = c(size,nrow(tmp))
    Zopt = c(Zopt,j)
    Abopt = c(Abopt,i)
  }
  print(i)
}
Opt = data.frame("Ratio"=ratio,"Size"=size,
                  "Ab"=Abopt,"Z"=Zopt)
Opt = Opt[Opt$Ratio<0.1,]
Opt = Opt[-1,]
Opt = Opt[!is.na(Opt$Size),]
max(Opt$Size)
Opt[which(Opt$Size==max(Opt$Size)),]

vline = data.frame(x1 = 32, x2 = 32, y1 = 3, y2 = 4.85)
hline = data.frame(x1 = 32, x2 = Inf,y1 = 4.85, y2 = 4.85)
q = ggplot(data=PexPCA3,aes(x=MaxAb,y=ZscoreMinN,col=In_literature)) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "Boundary"), show_guide=F, size = 1, data = vline) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "Boundary"), show_guide=F, size = 1, data = hline) +
  geom_point() +
  #scale_colour_manual(values = c("#d95f02","#7570b3","#1b9e77")) +
  #scale_colour_manual(values = c("#d95f02","#e6e600","blue")) +
  scale_colour_manual(values = c("#d95f02","purple","#e6e600")) +
  xlab("Abundance") +
  ylab("Z score") +
  ylim(c(3,23)) +
  theme_bw()
ggsave("C:/Users/Bram/Documents/Project SWAT DHFR PCA/SWAT Abundance vs Z score.png",q)

write.csv(PexPCA3,"C:/Users/Bram/Documents/Project SWAT DHFR PCA/Literature fig data.csv")

# Label positive results
PexPCA2$Pos = rep(0,nrow(PexPCA2))
PexPCA2$Pos[PexPCA2$ZscoreMinN > 3] = 1
PexPCA2$Pos[PexPCA2$MaxAb > 30 & PexPCA2$ZscoreMinN < 4.85] = 0

# Plot Z scores against each other
q = ggplot(PexPCA2,aes(x=Zscore1N,y=Zscore2N)) +
  geom_point(color="black",alpha=0.2) +
  geom_abline(slope=1,intercept=0,color="blue",lwd=1) +
  xlab("Z score replicate 1") +
  ylab("Z score replicate 2") +
  theme_bw()
ggsave("C:/Users/Bram/Documents/Project SWAT DHFR PCA/Z scores comparison",q,width=7,height=7)

# Plot ranked protein-protein interaction scores vs literature data
toRemove = PexPCA3$ID[PexPCA3$MaxAb > 30 & PexPCA3$ZscoreMinN < 4.85]
PexPCA2bu = PexPCA2
PexPCA2fig = PexPCA2[-which(is.element(PexPCA2$ID,toRemove)),]

PexPCA2fig = PexPCA2fig[order(PexPCA2fig$ZscoreMinN,decreasing=T),]
Spread = data.frame("Zscore_rank"=which(PexPCA2fig$Int_in_literature>0))
q = ggplot(Spread,aes(x=Zscore_rank)) +
  geom_histogram(binwidth=100,col="Blue",fill="Blue") +
  ylab("Percentage interactions found in literature") +
  xlab("Z score rank") +
  theme_bw()
ggsave("C:/Users/Bram/Documents/Project SWAT DHFR PCA/Z scores vs literature.png",q,width=7,height=6)  

PexPCA2 = PexPCA2[order(PexPCA2$ZscoreMinN,decreasing=T),]

# Number of unique proteins involved in positive interactions
Sel1Pos = PexPCA2[PexPCA2$Pos == 1,]
Involved = c(Sel1Pos$Bait,Sel1Pos$Prey)
Involved = Involved[order(Involved,decreasing=T)]
table(Involved) # 49 PROTEINS

# Number of unique interactions
UniqueInt = c()
for (i in 1:nrow(Sel1Pos)) {
  tmp = c(Sel1Pos$Bait[i],Sel1Pos$Prey[i])
  tmp = tmp[order(tmp)]
  UniqueInt = c(UniqueInt,paste(tmp[1],tmp[2]))
}
hist(table(UniqueInt))
length(unique(UniqueInt)) # 109 UNIQUE INTERACTIONS FOR 243 POSITIVE RESULTS

# Count occurrences of SWAT and original library strains in positive results
TotalNN = length(intersect(which(PexPCA2$Baittag == "N"),which(PexPCA2$Preytag == "N")))
TotalNC = length(intersect(which(PexPCA2$Baittag == "N"),which(PexPCA2$Preytag == "C")))
TotalCN = length(intersect(which(PexPCA2$Baittag == "C"),which(PexPCA2$Preytag == "N")))
TotalCC = length(intersect(which(PexPCA2$Baittag == "C"),which(PexPCA2$Preytag == "C")))
PosNN = length(intersect(which(Sel1Pos$Baittag == "N"),which(Sel1Pos$Preytag == "N")))
PosNC = length(intersect(which(Sel1Pos$Baittag == "N"),which(Sel1Pos$Preytag == "C")))
PosCN = length(intersect(which(Sel1Pos$Baittag == "C"),which(Sel1Pos$Preytag == "N")))
PosCC = length(intersect(which(Sel1Pos$Baittag == "C"),which(Sel1Pos$Preytag == "C")))
PerfNN = PosNN/TotalNN
PerfNC = PosNC/TotalNC
PerfCN = PosCN/TotalCN
PerfCC = PosCC/TotalCC
# Performance NN: 0.64% of tested interactions are positive
# Performance NC: 0.94% of tested interactions are positive
# Performance CN: 0.98% of tested interactions are positive
# Performance CC: 1.57% of tested interactions are positive

### PART 7: CREATION OF HEAT MAP FOR ARTICLE

# To understand which unique interactions required inclusion of SWAT strains and
# create a square matrix for the heat map:
Bmat = c()
Pmat = c()
Npos = c()
nHits = c()
nTrials = c()
PexPCA2 = PexPCA2[PexPCA2$Pos == 1,]
RefIDs = unique(c(as.character(PexPCA2$Bait),as.character(PexPCA2$Prey)))
Zsc = c()
PexPCA2 = PexPCA2[!is.na(PexPCA2$PreyID),]
for (i in 1:length(RefIDs)) {
  for (j in i:length(RefIDs)) {
    B = RefIDs[i]
    P = RefIDs[j]
    Bmat = c(Bmat,B)
    Pmat = c(Pmat,P)
    tmp = PexPCA2[(as.character(PexPCA2$Bait)==B & as.character(PexPCA2$Prey)==P) |
                    (as.character(PexPCA2$Bait)==P & as.character(PexPCA2$Prey)==B),]
    nTrials = c(nTrials,nrow(tmp))
    tmp = tmp[tmp$ZscoreMinN > 3,]
    if (nrow(tmp) > 0) {
      
      tmp2 = c()
      for (k in 1:nrow(tmp)) {
        tmp2 = c(tmp2,paste(as.character(tmp$Baittag[k]),as.character(tmp$Preytag[k])))
      }
      if (is.element("C C",tmp2)) { Npos = c(Npos,"NoN") } else { Npos = c(Npos,"N") }
      nHits = c(nHits,nrow(tmp))
      print(i)
      Zsc = c(Zsc,max(tmp$ZscoreMinN))
      
    } else {
      Npos = c(Npos,"NoN")
      nHits = c(nHits,0)
      Zsc = c(Zsc,0)
    }
    
  }
}

# Data frame for heat map
fire = data.frame("Bait"=Bmat,"Prey"=Pmat,"nTrials"=nTrials,"nHits"=nHits,"Zscore"=Zsc,"Npos"=Npos)
fire$ref = 1:nrow(fire)

# Duplicate results:
fire3 = fire
fire3$Bait = fire3$Prey
fire3$Prey = fire$Bait
fire3 = rbind(fire,fire3)

# Focus on what matters:
fire4 = fire3[,c(1,2,4,6)]
fire4$Nspec = rep(0,nrow(fire4))
fire4$Nspec[fire4$Npos == "N"] = 2
fire4$Nspec[fire4$Npos == "NoN" & fire4$nHits > 0] = 1
fire4 = fire4[,-c(3:4)]
fire4[,3] = as.numeric(as.character(fire4[,3]))
w <- reshape(fire4, 
             timevar = "Prey",
             idvar = "Bait",
             direction = "wide")
rownames(w) = w[,1]
w = w[,-1]
wm = as.matrix(w)
wm = wm[,order(colnames(wm))]
wm = wm[order(rownames(wm)),]
for (i in 1:ncol(wm)) {
    colnames(wm)[i] = substr(colnames(wm)[i],start=7,stop=nchar(colnames(wm)[i]))}

# After reorganizations, create the heat map distinguishing results obtained
# without the need of SWAT strains and those that required SWAT strains:
palette <- colorRampPalette(c('#0000ff','#0000FF'),bias=2000)(256)
palette = c("black",palette)
palette = c("black","blue","yellow")
#palette = c("black","#258039","#F5BE41")
cairo_ps(filename = "Heat map 1.ps" ,
    width = 7, height = 7, pointsize = 12)
plot.new()
heatmap(wm, Colv = NA, Rowv = NA, symm = T, col = palette, cexRow = 0.01, cexCol = 0.01)
dev.off()
tiff(filename = "Heat map 1 with label.tif", res = 200)
plot.new()
heatmap(wm, Colv = NA, Rowv = NA, symm = T, col = palette)
dev.off()

write.csv(wm,"C:/Users/Bram/Documents/wm.csv")
tiff(filename = "Heat map 1.tif", width = 5, height = 5, res = 250)
plot.new()
heatmap(wm, Colv = NA, Rowv = NA, symm = T, col = palette, margins = c(0,0), 
        labRow = NA, labCol = NA)
dev.off()
# Left part of the map: fire5
# Right part of the map: fire4
fire5 = fire3[,c(1,2,5)]
fire5[,3] = as.numeric(as.character(fire5[,3]))
w2 <- reshape(fire5, 
             timevar = "Prey",
             idvar = "Bait",
             direction = "wide")
rownames(w2) = w2[,1]
w2 = w2[,-1]
wm2 = as.matrix(w2)
wm2 = wm2[,order(colnames(wm2))]
wm2 = wm2[order(rownames(wm2)),]
for (i in 1:ncol(wm2)) {
  colnames(wm2)[i] = substr(colnames(wm2)[i],start=8,stop=nchar(colnames(wm2)[i]))}

write.csv(wm2,"C:/Users/Bram/Documents/wm2.csv")

palette <- colorRampPalette(c('#cdffcd','#00FF00'),bias=20000000)(256)
palette <- colorRampPalette(c('white','#5f00d3'),bias=20000000)(256)
palette <- colorRampPalette(c('white','red'),bias=2000)(256)
palette = c("black",palette)
tiff(filename = "Heat map 2 green comp.tif", res = 72)
plot.new()
heatmap(wm2, Colv = NA, Rowv = NA, symm = T, col = palette, margins = c(0,0), 
        labRow = NA, labCol = NA)
dev.off()

# Save all interaction data
PexPCA2bu2 = PexPCA2bu[,c(4,5,6,8,9,10,15,16,13,14,17,18,20,25,26,29,31,34)]
PexPCA2bu2 = PexPCA2bu2[order(PexPCA2bu2$ZscoreMinN,decreasing=T),]

write.csv(PexPCA2bu2,"C:/Users/Bram/Documents/Project SWAT DHFR PCA/SWAT DHFR PCA all results 10-7-17.csv")

# Different format to show positive results

Nreq = c()
for (i in 1:nrow(PexPCA2)) {
  tmp = intersect(which(fire3$Bait == PexPCA2$Bait[i]),which(fire3$Prey == PexPCA2$Prey[i]))
  if (length(tmp) > 1) {print(fire3[tmp,])}
  Nreq = c(Nreq,fire3$Npos[tmp[1]])
}
Nreq = c(Nreq,Nreq)

Baitlist = c()
Preylist = c()
for (i in 1:nrow(PexPCA2)) {
  Baitlist = c(Baitlist,PexPCA2$Bait[i])
  if (Nreq[i] == 1) {
    Preylist = c(Preylist,paste("aaa",PexPCA2$Prey[i],sep=""))
  } else {
    Preylist = c(Preylist,PexPCA2$Prey[i])
  }
}
for (i in 1:nrow(PexPCA2)) {
  Baitlist = c(Baitlist,PexPCA2$Prey[i])
  if (Nreq[i] == 1) {
    Preylist = c(Preylist,paste("aaa",PexPCA2$Bait[i],sep=""))
  } else {
    Preylist = c(Preylist,PexPCA2$Bait[i])
  }
}

results = list()
for (i in 1:length(unique(Baitlist))) {
  tmp = which(Baitlist==unique(Baitlist)[i])
  results[[i]] = unique(Preylist[tmp])
  print(length(results[[i]]))
}

final = matrix("NA",nrow=length(unique(Baitlist)),ncol=21)
final[,1] = unique(Baitlist)
for (i in 1:length(results)) {
  for (j in 1:length(results[[i]])) {
    final[i,j+1] = results[[i]][j]
  }
}

write.csv(final,"C:/Users/Bram/Documents/Project SWAT DHFR PCA/SWAT DHFR PCA different format 10-7-17.csv")

PexPCA2newscreen = PexPCA2
