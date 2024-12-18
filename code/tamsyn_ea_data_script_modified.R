#EA (DEFRA) water sampling data analysis

#Downloaded complete set of Devon and Cornwall water sampling results for 2011-2021 and saved as csvs
#https://environment.data.gov.uk/water-quality/view/download/new

#They contain a lot of excess data, so filter and append datasets together

#read in data

{
 		# input oldest year
 		oldest <- 2013
 		# input newest year
 		newest <- 2024
 		
 		# initialise data storage list
 		ea_data <- list()
 		# step through each year
 		for (i in oldest:newest) {
 			# sequentially name and import each dataset
 			ea_data[[i]] <- read.csv(
 				paste0(
 					"raw_data/DC-", i, ".csv"
 				)
 			)
 		}
 		# use do.call to rbind each list item together
 		ea_data <- do.call(
 			rbind, ea_data
  )
}

#some serious filtering
#categories: ANY AGRICULTURAL	ANY BIOTA	ANY FISH - NOT INCLUDING FLATFISH - WHOLE ANIMAL	ANY LEACHATE	ANY NON-AQUEOUS LIQUID	ANY OIL	ANY SEWAGE	ANY SHELLFISH (MOLLUSC & CRUSTACEAN) - WHOLE ANIMAL	ANY SOLID/SEDIMENT - UNSPECIFIED	ANY TRADE EFFLUENT	ANY WATER	CANAL WATER	COASTAL / MARINE SEDIMENT	CRUDE SEWAGE	ESTUARINE WATER	ESTUARY SEDIMENT	ESTUARY SEDIMENT - SUB TIDAL	FINAL SEWAGE EFFLUENT	GROUNDWATER	GROUNDWATER - PURGED/PUMPED/REFILLED	GROUNDWATER - STATIC/UNPURGED	MINEWATER	MINEWATER (FLOWING/PUMPED)	MYTILUS EDULIS - MUSSEL - MUSCLE	MYTILUS EDULIS - MUSSEL - WHOLE ANIMAL	POND / LAKE / RESERVOIR WATER	RIVER / RUNNING SURFACE WATER	SALMO TRUTTA - BROWN TROUT - WHOLE ANIMAL	SEA WATER	STORM SEWER OVERFLOW DISCHARGE	SURFACE DRAINAGE	TRADE EFFLUENT - FRESHWATER RETURNED ABSTRACTED	TRADE EFFLUENT - SALINE WATER RETURNED ABSTRACTED
#of interest: ANY WATER ESTUARINE WATER GROUNDWATER GROUNDWATER - PURGED/PUMPED/REFILLED GROUNDWATER - STATIC/UNPURGED MINEWATER MINEWATER (FLOWING/PUMPED) POND / LAKE / RESERVOIR WATER RIVER / RUNNING SURFACE WATER SURFACE DRAINAGE


#1 by sample type: using grep allows match when there are different categories i.e. of groundwater

{
MINEWATER <- ea_data[grep("MINEWATER", ea_data$sample.sampledMaterialType.label), ]
GROUNDWATER <- ea_data[grep("GROUNDWATER", ea_data$sample.sampledMaterialType.label), ]
RIVER <- ea_data[grep("RIVER / RUNNING SURFACE WATER", ea_data$sample.sampledMaterialType.label), ]

ALL_WATER <- ea_data[grep("WATER", ea_data$sample.sampledMaterialType.label), ]
}

#2 by sampling location i.e. obtain list of sampling sites per catchment
#downloaded ids (grid referecnes) from here: https://environment.data.gov.uk/water-quality/view/explore?search=hayle&area=6-27&samplingPointType.group=&samplingPointStatus%5B%5D=open&loc=157079%2C33797&_limit=500
#these aren't grid references but EA ids, they also don't match the ids from the river basin data i.e. at https://environment.data.gov.uk/catchment-planning/WaterBody/GB108049000560 
#but the SW...IDs are the ones in the database, so need to use them

{
###1 HAYLE
Hayle <- read.csv("raw_data/sp_hayle.csv" , header=T)
HayleIDs <- Hayle$notation
RIVER_Hayle <- RIVER[grepl(paste(HayleIDs, collapse="|"), RIVER$sample.samplingPoint.notation), ]

RIVER_Hayle <- ALL_WATER[grepl(paste(HayleIDs, collapse="|"), ALL_WATER$sample.samplingPoint.notation), ]

###2 RED RIVER

RR <- read.csv("raw_data/sp_red.csv" , header=T)
RRIDs <- RR$notation
Red_RIVER <- RIVER[grepl(paste(RRIDs, collapse="|"), RIVER$sample.samplingPoint.notation), ]

Red_RIVER <- ALL_WATER[grepl(paste(RRIDs, collapse="|"), ALL_WATER$sample.samplingPoint.notation), ]

####3 CARNON RIVER

CR <- read.csv("raw_data/sp_carnon.csv" , header=T)
CRIDs <- CR$notation
Carnon_RIVER <- RIVER[grepl(paste(CRIDs, collapse="|"), RIVER$sample.samplingPoint.notation), ]

Carnon_RIVER <- ALL_WATER[grepl(paste(CRIDs, collapse="|"), ALL_WATER$sample.samplingPoint.notation), ]


####4 RIVER COBER

CoR <- read.csv("raw_data/sp_cober.csv" , header=T)
CoRIDs <- CoR$notation
Cober_RIVER <- RIVER[grepl(paste(CoRIDs, collapse="|"), RIVER$sample.samplingPoint.notation), ]

Cober_RIVER <- ALL_WATER[grepl(paste(CoRIDs, collapse="|"), ALL_WATER$sample.samplingPoint.notation), ]
}

#now all sorts of possibilities, including multivariate distance etc with different parameters, but focus on metals for now

#3 by measurment type- focusing on metals and nutrients
# made a list of interesting values from 2016 dataset in determinand.label column
{
metals_nut <- selected_params$selected_params

Hayle_met_nut <- RIVER_Hayle[grepl(paste(metals_nut, collapse="|"), RIVER_Hayle$determinand.label), ]
Red_met_nut <- Red_RIVER[grepl(paste(metals_nut, collapse="|"), Red_RIVER$determinand.label), ]
Carnon_met_nut <- Carnon_RIVER[grepl(paste(metals_nut, collapse="|"), Carnon_RIVER$determinand.label), ]
Cober_met_nut <- Cober_RIVER[grepl(paste(metals_nut, collapse="|"), Cober_RIVER$determinand.label), ]

##tidy up the data, get rid of useless columns

colnames(Hayle_met_nut)
Hayle_met_nut <- Hayle_met_nut[ -c(1,2,4,7:9,11,14:17) ]
Red_met_nut <- Red_met_nut[ -c(1,2,4,7:9,11,14:17) ]
Carnon_met_nut <- Carnon_met_nut[ -c(1,2,4,7:9,11,14:17) ]
Cober_met_nut <- Cober_met_nut[ -c(1,2,4,7:9,11,14:17) ]
}

#make a new column for a unique sample i.e. location&type&time and to add units to measuremnt

library(stringr)

#split sample time column into date,year and rest
{
Hayle_met_nut[c('Year', 'Month', 'Time')] <- str_split_fixed(Hayle_met_nut$sample.sampleDateTime, '-', 3)
Red_met_nut[c('Year', 'Month', 'Time')] <- str_split_fixed(Red_met_nut$sample.sampleDateTime, '-', 3)
Carnon_met_nut[c('Year', 'Month', 'Time')] <- str_split_fixed(Carnon_met_nut$sample.sampleDateTime, '-', 3)
Cober_met_nut[c('Year', 'Month', 'Time')] <- str_split_fixed(Cober_met_nut$sample.sampleDateTime, '-', 3)

Hayle_met_nut$Sample_Type_Time <- paste(Hayle_met_nut$sample.samplingPoint.notation,Hayle_met_nut$Year, Hayle_met_nut$Month, sep = "_")
Hayle_met_nut$Measure_Unit <- paste(Hayle_met_nut$determinand.label, Hayle_met_nut$determinand.unit.label, sep="_")

Red_met_nut$Sample_Type_Time <- paste(Red_met_nut$sample.samplingPoint.notation,Red_met_nut$Year, Red_met_nut$Month, sep = "_")
Red_met_nut$Measure_Unit <- paste(Red_met_nut$determinand.label, Red_met_nut$determinand.unit.label, sep="_")

Carnon_met_nut$Sample_Type_Time <- paste(Carnon_met_nut$sample.samplingPoint.notation,Carnon_met_nut$Year, Carnon_met_nut$Month, sep = "_")
Carnon_met_nut$Measure_Unit <- paste(Carnon_met_nut$determinand.label, Carnon_met_nut$determinand.unit.label, sep="_")

Cober_met_nut$Sample_Type_Time <- paste(Cober_met_nut$sample.samplingPoint.notation,Cober_met_nut$Year, Cober_met_nut$Month, sep = "_")
Cober_met_nut$Measure_Unit <- paste(Cober_met_nut$determinand.label, Cober_met_nut$determinand.unit.label, sep="_")

#now get rid of the individual columsn

colnames(Hayle_met_nut)
Hayle_met_nut <- Hayle_met_nut[ -c(1:3,5:9) ]
Red_met_nut <- Red_met_nut[ -c(1:3,5:9) ]
Carnon_met_nut <- Carnon_met_nut[ -c(1:3,5:9) ]
Cober_met_nut <- Cober_met_nut[ -c(1:3,5:9) ]
}


#####WIDE FORMAT

#resahpe into wide format- with all the different measured parametrs for each sampling point

wide_hayle <- reshape(Hayle_met_nut, idvar = "Sample_Type_Time", timevar = "Measure_Unit", direction = "wide")
wide_red <- reshape(Red_met_nut, idvar = "Sample_Type_Time", timevar = "Measure_Unit", direction = "wide")
wide_carnon <- reshape(Carnon_met_nut, idvar = "Sample_Type_Time", timevar = "Measure_Unit", direction = "wide")
wide_cober <- reshape(Cober_met_nut, idvar = "Sample_Type_Time", timevar = "Measure_Unit", direction = "wide")

#now split back so time is separate to sample location & site, and add them back again

library("stringr")

library("ggplot2")

new_cols <- as.data.frame(str_split_fixed(wide_hayle$Sample_Type_Time, "_", 2))
wide_hayle$Sample <-paste(new_cols$V1, new_cols$V2, sep="_")
wide_hayle$Site <- new_cols$V1
wide_hayle$Time <- new_cols$V2

new_cols1 <- as.data.frame(str_split_fixed(wide_red$Sample_Type_Time, "_", 2))
wide_red$Sample <-paste(new_cols1$V1, new_cols1$V2, sep="_")
wide_red$Site <- new_cols1$V1
wide_red$Time <- new_cols1$V2

new_cols2 <- as.data.frame(str_split_fixed(wide_carnon$Sample_Type_Time, "_", 2))
wide_carnon$Sample <-paste(new_cols2$V1, new_cols2$V2, sep="_")
wide_carnon$Site <- new_cols2$V1
wide_carnon$Time <- new_cols2$V2

new_cols3 <- as.data.frame(str_split_fixed(wide_cober$Sample_Type_Time, "_", 2))
wide_cober$Sample <-paste(new_cols3$V1, new_cols3$V2, sep="_")
wide_cober$Site <- new_cols3$V1
wide_cober$Time <- new_cols3$V2

#Now lets plot

head(wide_carnon)

#sort out headers (should do this earlier)
colnames(wide_hayle) <- gsub("\\s+", "", colnames(wide_hayle))
colnames(wide_hayle) <- gsub("-", "_", colnames(wide_hayle))
colnames(wide_hayle) <- gsub("/", ".", colnames(wide_hayle))

colnames(wide_red) <- gsub("\\s+", "", colnames(wide_red))
colnames(wide_red) <- gsub("-", "_", colnames(wide_red))
colnames(wide_red) <- gsub("/", ".", colnames(wide_red))

colnames(wide_carnon) <- gsub("\\s+", "", colnames(wide_carnon))
colnames(wide_carnon) <- gsub("-", "_", colnames(wide_carnon))
colnames(wide_carnon) <- gsub("/", ".", colnames(wide_carnon))

colnames(wide_cober) <- gsub("\\s+", "", colnames(wide_cober))
colnames(wide_cober) <- gsub("-", "_", colnames(wide_cober))
colnames(wide_cober) <- gsub("/", ".", colnames(wide_cober))

#plotting just metals- whole metals (not filtered/dissolved-these are available for all as well)
####not added cober yet

##Hayle

wide_hayle$Site <- factor(wide_hayle$Site,levels = c("SW-82221895","SW-82221891","SW-82221881","SW-82221871","SW-82220474","SW-82221863","SW-82221859","SW-82221849","SW-82221844","SW-82221813","SW-82211838", "SW-82220303"))

colnames(wide_hayle)

wide_hayle$Site

Cd <- ggplot(wide_hayle, aes(x = Site, y = result.Cadmium_Cd_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Cadmium (ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Pb <- ggplot(wide_hayle, aes(x = Site, y = result.Lead_asPb_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Lead(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Zn <- ggplot(wide_hayle, aes(x = Site, y = result.Zinc_asZn_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Zinc(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Cu <- ggplot(wide_hayle, aes(x = Site, y = result.Copper_Cu_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Copper(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

As <- ggplot(wide_hayle, aes(x = Site, y = result.Arsenic_As_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Arsenic(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Fe <- ggplot(wide_hayle, aes(x = Site, y = result.Iron_asFe_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Iron(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Al <- ggplot(wide_hayle, aes(x = Site, y = result.Aluminium_Al_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Aluminium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Ni <- ggplot(wide_hayle, aes(x = Site, y = result.Nickel_Ni_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Nickel(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

Cr <- ggplot(wide_hayle, aes(x = Site, y = result.Chromium_Cr_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Chromium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge","SW-82221881" = "Drym Farm","SW-82221871" = "Binner Bridge","SW-82220474" = "Great Work Adit","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station","SW-82211838" = "River Hayle D/S Hayle Stw", "SW-82220303" = "Millpool stream"))

library(gridExtra)

grid.arrange(Cd,Cu,Zn, nrow=1)
grid.arrange(Pb, As, Fe, nrow=1)
grid.arrange (Al, Ni, Cr, nrow=1)

##Red River

wide_red$Site <- factor(wide_red$Site,levels = c("SW-82310101","SW-82310103","SW-82310132","SW-82310160","SW-82310183","SW-82310186","SW-82310704","SW-82310021","SW-82310411","SW-82310530","SW-82310608","SW-82310471","SW-82310477","SW-82310570","SW-82310609","SW-82310656","SW-82310658","SW-82310660"))

Cd_red <- ggplot(wide_red, aes(x = Site, y = result.Cadmium_Cd_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Cadmium (ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))


Pb_red <- ggplot(wide_red, aes(x = Site, y = result.Lead_asPb_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Lead(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Zn_red <- ggplot(wide_red, aes(x = Site, y = result.Zinc_asZn_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Zinc(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Cu_red <- ggplot(wide_red, aes(x = Site, y = result.Copper_Cu_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Copper(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

As_red <- ggplot(wide_red, aes(x = Site, y = result.Arsenic_As_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Arsenic(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Fe_red <- ggplot(wide_red, aes(x = Site, y = result.Iron_asFe_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Iron(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Al_red <- ggplot(wide_red, aes(x = Site, y = result.Aluminium_Al_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Aluminium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Ni_red <- ggplot(wide_red, aes(x = Site, y = result.Nickel_Ni_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Nickel(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

Cr_red <- ggplot(wide_red, aes(x = Site, y = result.Chromium_Cr_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Chromium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310021" = "Tehidy","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310530" = "Barripper (Roseworthy)","SW-82310608" = "Reen (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310477" = "Trevoole  (Roseworthy)","SW-82310570" = "Praze (Roseworthy)","SW-82310609" = "Gibralter discharge (Reen)","SW-82310656" = "Treslothan (Reen)","SW-82310658" = "Treslothan woods (Reen)","SW-82310660" = "U/s Treslothan adit (Reen)"))

library(gridExtra)

grid.arrange(Cd_red,Cu_red,Zn_red, nrow=1)
grid.arrange(Pb_red, As_red, Fe_red, nrow=1)
grid.arrange (Al_red, Ni_red, Cr_red, nrow=1)

##Carnon River

wide_carnon$Site <- factor(wide_carnon$Site,levels = c("SW-81950522","SW-81950610","SW-81950702","SW-81950555","SW-81950570","SW-81950622","SW-81950636","SW-81951044","SW-81950650","SW-81950659","SW-81950670","SW-81950421","SW-81950814","SW-81950840","SW-81950880","SW-81950888","SW-81950890","SW-81950892","SW-81950898","SW-81950899"))

Cd_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Cadmium_Cd_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Cadmium (ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Pb_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Lead_asPb_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Lead(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Zn_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Zinc_asZn_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Zinc(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Cu_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Copper_Cu_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Copper(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

As_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Arsenic_As_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Arsenic(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Fe_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Iron_asFe_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Iron(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Al_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Aluminium_Al_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Aluminium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

  scale_x_discrete(labels=c("SW-81950622" = "Twelveheads D/S St Day Stream","SW-81950636" = "Twelveheads U/S St Day Stream","SW-81950709" = "Carnon Cons Mines Clemows Valley Dam","SW-81950594" = "Carnon D/S County Adit Culverted Stretch","SW-81950512" = "Carnon Downs (New) Stw Fe","SW-81950513" = "Carnon Downs (New) Stw Inlet","SW-81950555" = "Bissoe Br Gauging Stn","SW-81950659" = "Chacewater Viaduct","SW-81950522" = "Devoran Bridge","SW-81950640" = "Wheal Prosper","SW-81950666" = "D/S Blackwater Stw","SW-81950670" = "U/S Blackwater Stw","SW-81950650" = "U/S Chacewater Stw","SW-81950514" = "U/S Devoran Stw","SW-81950702" = "Clemmows Stream U/S Carnon River","SW-81950503" = "Point Mills","SW-81950570" = "D/S County/Wellington Adits","SW-81950610" = "U/S Grenna Bridge","SW-81951044" = "St. Day Stream U/S Carnon River", "SW-81950896" = "Glenmoor Farm Stream","SW-81950899" = "Hicks Church Row","SW-81950890" = "Comford","SW-81950898" = "Coppice Gardens","SW-81950888" = "Gwennap Church","SW-81950814" = "Hicks Mill","SW-81950894" = "Penpons Mill","SW-81950892" = "Penventon","SW-81950840" = " Trehaddle","SW-81950880" = "D/S Lanner St Day Stw","SW-81950421" = " D/S Points Mill Bridge","SW-81950803" = "U/S Bissoe G Stn","SW-81950886" = "U/S Lanner St Day Stw","SW-81950897" = "Tresavean Leat"))

Ni_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Nickel_Ni_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Nickel(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

Cr_carnon <- ggplot(wide_carnon, aes(x = Site, y = result.Chromium_Cr_ug.l, fill =)) + 
  geom_boxplot(fill= "dodgerblue")+
  theme_bw(base_size = 10)+
  theme(legend.justification=c(0.02,0.95), legend.position=c(0.01,0.98), legend.title = element_text(size=12), legend.text = element_text(size=12), legend.background = element_rect(linetype="solid", size=.3, colour="black") )+
  theme(axis.text = element_text(size = 8), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Chromium(ug/l)") + xlab(" ")+
  scale_x_discrete(labels=c("SW-81950522" = "Devoran Bridge","SW-81950610" = "Grenna Bridge","SW-81950702" = "Clemmows Stream","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950650" = "Chacewater","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))



library(gridExtra)

grid.arrange(Cd_carnon,Cu_carnon,Zn_carnon, nrow=1)
grid.arrange(Pb_carnon, Cr_carnon, Fe_carnon, nrow=1)
grid.arrange (Al_carnon, Ni_carnon, As_carnon, nrow=1)

#River comparsions

grid.arrange(Cd, Cd_red, Cd_carnon, nrow=1)
grid.arrange(Cu, Cu_red, Cu_carnon, nrow=1)
grid.arrange(Zn, Zn_red, Zn_carnon, nrow=1)
grid.arrange(Pb, Pb_red, Pb_carnon, nrow=1)
grid.arrange(As, As_red, As_carnon, nrow=1)
grid.arrange(Fe, Fe_red, Fe_carnon, nrow=1)
grid.arrange(Al, Al_red, Al_carnon, nrow=1)
grid.arrange(Ni, Ni_red, Ni_carnon, nrow=1)
grid.arrange(Cr, Cr_red, Cr_carnon, nrow=1)

#associations between metals
#HAYLE
#all
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Aluminium_Al_ug.l + result.Arsenic_As_ug.l + result.Lead_asPb_ug.l + result.Nickel_Ni_ug.l + result.Iron_asFe_ug.l + result.Cr_Filtered_ug.l + result.pH_phunits, data = wide_hayle, col = rainbow(11)[wide_hayle$Site])

#selected 
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l + result.Arsenic_As_ug.l + log(result.Iron_asFe_ug.l) + result.Lead_asPb_ug.l + result.Aluminium_Al_ug.l, data = wide_hayle, col = rainbow(11)[wide_hayle$Site])
      
#best
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l, data = wide_hayle, col = rainbow(11)[wide_hayle$Site])

#RED
#all
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Aluminium_Al_ug.l + result.Arsenic_As_ug.l + result.Lead_asPb_ug.l + result.Nickel_Ni_ug.l + result.Iron_asFe_ug.l + result.Cr_Filtered_ug.l + result.pH_phunits, data = wide_hayle, col = rainbow(11)[wide_red$Site])

#selected 
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l + result.Arsenic_As_ug.l + log(result.Iron_asFe_ug.l) + result.Lead_asPb_ug.l + result.Aluminium_Al_ug.l, data = wide_hayle, col = rainbow(11)[wide_red$Site])

#best
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l, data = wide_hayle, col = rainbow(11)[wide_red$Site])


#CARNON
#all
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Aluminium_Al_ug.l + result.Arsenic_As_ug.l + result.Lead_asPb_ug.l + result.Nickel_Ni_ug.l + result.Iron_asFe_ug.l + result.Cr_Filtered_ug.l + result.pH_phunits, data = wide_hayle, col = rainbow(11)[wide_carnon$Site])

#selected 
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l + result.Arsenic_As_ug.l + log(result.Iron_asFe_ug.l) + result.Lead_asPb_ug.l + result.Aluminium_Al_ug.l, data = wide_hayle, col = rainbow(11)[wide_carnon$Site])

#best
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l + result.Nickel_Ni_ug.l, data = wide_hayle, col = rainbow(11)[wide_carnon$Site])



#individual sites of interest
#for each catchemnet, the columns are different!
#filter out the interesting metals and nutrients (may be more to add)
#add sample_time as the rownames (ids)

rownames(wide_hayle) <- wide_hayle$Sample_Type_Time
colnames(wide_hayle)
metnut_hayle <- wide_hayle[,c(45,46,47,13,18,10,31,23,30,33,17,11,4,14,21,27,32,7,36)]

rownames(wide_red) <- wide_red$Sample_Type_Time
colnames(wide_red)
metnut_red <- wide_red[,c(51,52,53,8,7,6,5,35,16,15,3,25,17,22,2,20,23,12,21)]

rownames(wide_carnon) <- wide_carnon$Sample_Type_Time
colnames(wide_carnon)
metnut_carn <- wide_carnon[,c(44,45,46,12,3,9,20,8,2,16,34,25,27,31,37,14,7,29,23)]

rownames(wide_cober) <- wide_cober$Sample_Type_Time
colnames(wide_cober)
metnut_cob <- wide_cober[,c(47,48,49,14,7,11,36,38,5,32,26,16,10,25,24,28,35,2,18)]

colnames(metnut_carn)
colnames(metnut_hayle)
colnames(metnut_red)
colnames(metnut_cob)

#selecting sites

abcrowan <- metnut_hayle[grep("SW-82221895", metnut_hayle$Site),] 
colMeans(abcrowan[ , c(4:12)], na.rm = TRUE)

sterth <- metnut_hayle[grep("SW-82221813", metnut_hayle$Site),] 
colMeans(sterth[ , c(4:12)], na.rm = TRUE)

Relubbas <- metnut_hayle[grep("SW-82221844", metnut_hayle$Site),] 
colMeans(Relubbas[ , c(4:12)], na.rm = TRUE)

God <- metnut_hayle[grep("SW-82221863", metnut_hayle$Site),] 
colMeans(God[ , c(4:12)], na.rm = TRUE)

Drym <- metnut_hayle[grep("SW-82221871", metnut_hayle$Site),] 
colMeans(Drym[ , c(4:12)], na.rm = TRUE)

Crowan <- metnut_hayle[grep("SW-82221891", metnut_hayle$Site),] 
colMeans(Crowan[ , c(4:12)], na.rm = TRUE)

Brea <- metnut_red[grep("SW-82310186", metnut_red$Site),] 
colMeans(Brea[ , c(4:12)], na.rm = TRUE)

Gwith <- metnut_red[grep("SW-82310101", metnut_red$Site),] 
colMeans(Gwith[ , c(4:12)], na.rm = TRUE)

Tresk <- metnut_red[grep("SW-82310183", metnut_red$Site),] 
colMeans(Tresk[ , c(4:12)], na.rm = TRUE)

Ros <- metnut_red[grep("SW-82310160", metnut_red$Site),] 
colMeans(Ros[ , c(4:121)], na.rm = TRUE)

Kieve <- metnut_red[grep("SW-82310132", metnut_red$Site),] 
colMeans(Kieve[ , c(4:12)], na.rm = TRUE)

Grenna<- metnut_red[grep("SW-81950610", metnut_carn$Site),] 
colMeans(Grenna[ , c(4:12)], na.rm = TRUE)

Bissoe <- metnut_carn[grep("SW-81950555", metnut_carn$Site),] 
colMeans(Bissoe[ , c(4:12)], na.rm = TRUE)

Chace <- metnut_carn[grep("SW-81950659", metnut_carn$Site),] 
colMeans(Chace[ , c(4:12)], na.rm = TRUE)

LanSTC<- metnut_carn[grep("SW-81950880", metnut_carn$Site),] 
colMeans(LanSTC[ , c(4:12)], na.rm = TRUE)

Comm<- metnut_carn[grep("SW-81950890", metnut_carn$Site),] 
colMeans(Comm[ , c(4:12)], na.rm = TRUE)

Cop<- metnut_carn[grep("SW-81950898", metnut_carn$Site),] 
colMeans(Cop[ , c(4:12)], na.rm = TRUE)

Helston<- metnut_cob[grep("SW-82010152", metnut_cob$Site),]
colMeans(Helston[ , c(4:12)], na.rm = TRUE)

Lowertown<- metnut_cob[grep("SW-82010156", metnut_cob$Site),]
colMeans(Lowertown[ , c(4:12)], na.rm = TRUE)

Wendron<- metnut_cob[grep("SW-RSN0011", metnut_cob$Site),]
colMeans(Wendron[ , c(4:12)], na.rm = TRUE)

Trenear<- metnut_cob[grep("SW-82010187", metnut_cob$Site),]
colMeans(Trenear[ , c(4:12)], na.rm = TRUE)

Med<- metnut_cob[grep("SW-82011025", metnut_cob$Site),]
colMeans(Med[ , c(4:12)], na.rm = TRUE)


###

#to each add a colum with site name, then combine all sites for one river and add a column with river, then can combine all

abcrowan$Site2 <- 'Above Crowan'
sterth$Site2 <- 'StErth'
Relubbas$Site2 <- 'Relubbas'
God$Site2 <- 'Godolphin'
Drym$Site2 <- 'Drym'
Crowan$Site2 <- 'Crowan'

Tresk$Site2 <- 'Treskillard'
Ros$Site2 <- 'Roscroggan'
Kieve$Site2 <- 'Kieve bridge'
Gwith$Site2 <- 'Gwithian'
Brea$Site2 <- 'Brea'

Grenna$Site2 <- 'Grenna'
Bissoe$Site2 <- 'Bissoe'
Chace$Site2 <- 'Chacewater viaduct'
LanSTC$Site2 <- 'Lanner STW'
Comm$Site2 <- 'Comford'
Cop$Site2 <- 'Coppice gardens'

Helston$Site2 <-'Helston'
Lowertown$Site2 <- 'Lowertown'
Wendron$Site2 <- 'Wendron'
Trenear$Site2 <- 'Trenear' 
Med$Site2 <- 'Medlyn'

#adding sites to rivers

Hayle_sites <- rbind(Relubbas,God,Drym,Crowan,abcrowan,sterth)
Red_sites <- rbind(Tresk,Ros,Kieve,Gwith,Brea)
Carnon_sites <-rbind(Bissoe,Chace,LanSTC,Comm,Cop,Grenna)
Cober_sites <-rbind(Helston,Lowertown,Wendron,Trenear,Med)

#adding river names
Hayle_sites$River <- 4
Red_sites$River <- 1
Carnon_sites$River <- 8
Cober_sites$River <- 2

#adding all sites

All_sites <-rbind(Hayle_sites,Red_sites,Carnon_sites,Cober_sites)

#correlations between each factor, just for sites of interest and overall

pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Aluminium_Al_ug.l + result.Arsenic_As_ug.l + result.Lead_asPb_ug.l + result.Nickel_Ni_ug.l + result.Iron_asFe_ug.l + result.pH_phunits, data = All_sites, col = rainbow(12)[All_sites$Site])
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Nickel_Ni_ug.l, data = All_sites, col = rainbow(10)[All_sites$River])

pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Nickel_Ni_ug.l, data = Hayle_sites, col = rainbow(12)[Hayle_sites$Site])
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Nickel_Ni_ug.l, data = Carnon_sites, col = rainbow(20)[Carnon_sites$Site])
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Nickel_Ni_ug.l, data = Red_sites, col = rainbow(20)[Red_sites$Site])
pairs(~ result.Copper_Cu_ug.l + result.Zinc_asZn_ug.l + result.Cadmium_Cd_ug.l+ result.Nickel_Ni_ug.l, data = Red_sites, col = rainbow(6)[Cober_sites$Site])
#pca
# https://www.datacamp.com/tutorial/pca-analysis-r   

colnames(All_sites)
rownames(All_sites)

all_sites_metnut_pca <- prcomp(na.omit(All_sites[,c(4:19)]), center = TRUE, scale. = TRUE)
summary(all_sites_metnut_pca) #71, 81, 87 %

all_sites_met_pca <- prcomp(na.omit(All_sites[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(all_sites_met_pca) #68, 80, 89 %
str(all_sites_met_pca)

all_met_pca <- all_sites_met_pca$x
all_met_pca <- data.frame(sample = row.names(all_met_pca), all_met_pca)
all_met_pca[c('Site', 'Time')] <- str_split_fixed(all_met_pca$sample, '_', 2)

#library(githubinstall)
#githubinstall("vqv/ggbiplot")
library(ggbiplot)

#actually the ggbiplot does not seem great, can't find an option for effectively sclaing the arrows
#worth trying PCAtools package instead-YES
#we want to have arrows indicating strength of contribution of each factor, i.e. pH will not be so important I think

ggbiplot(all_sites_met_pca, groups = all_met_pca$Site, obs.scale = 1, var.scale = 1)


#now lets add in all the sites again, and perhaps do on a per river basis intially
#and consider removing outliers

#Hayle

hayle_met_pca <- prcomp(na.omit(metnut_hayle[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(hayle_met_pca) #51,78,90 % variance

hayle_nut_pca <- prcomp(na.omit(metnut_hayle[,c(12:17)]), center = TRUE, scale. = TRUE)
summary(hayle_nut_pca) #57,79,93 % variance

hayle_pca <- hayle_met_pca$x
hayle_pca <- data.frame(sample = row.names(hayle_pca), hayle_pca)
hayle_pca[c('Site', 'Time')] <- str_split_fixed(hayle_pca$sample, '_', 2)

hayle_npca <- hayle_nut_pca$x
hayle_npca <- data.frame(sample = row.names(hayle_npca), hayle_npca)
hayle_npca[c('Site', 'Time')] <- str_split_fixed(hayle_npca$sample, '_', 2)

ggbiplot(hayle_met_pca, groups = hayle_pca$Site, obs.scale = 1, var.scale = 1)
ggbiplot(hayle_nut_pca, groups = hayle_npca$Site, obs.scale = 1, var.scale = 1)
#("SW-82221895" = "Above Crowan","SW-82221891" = "Crowan Bridge", "SW-82221871" = "Binner Bridge","SW-82221863" = "Godolphin Bridge","SW-82221859" = "Wheal Godolphin Adit","SW-82221849" = "Trescowe Adit","SW-82221844" = "Relubbus","SW-82221813" = "St Erth Gauging Station", "SW-82220303" = "Millpool stream"))

#Red

red_met_pca <- prcomp(na.omit(metnut_red[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(red_met_pca) #41,58,73 % variance

red_nut_pca <- prcomp(na.omit(metnut_red[,c(12:17)]), center = TRUE, scale. = TRUE)
summary(red_nut_pca) #56,79,90 % variance

red_pca <- red_met_pca$x
red_pca <- data.frame(sample = row.names(red_pca), red_pca)
red_pca[c('Site', 'Time')] <- str_split_fixed(red_pca$sample, '_', 2)

red_npca <- red_nut_pca$x
red_npca <- data.frame(sample = row.names(red_npca), red_npca)
red_npca[c('Site', 'Time')] <- str_split_fixed(red_npca$sample, '_', 2)

ggbiplot(red_met_pca, groups = red_pca$Site, obs.scale = 1, var.scale = 1)
ggbiplot(red_nut_pca, groups = red_npca$Site, obs.scale = 1, var.scale = 1)
#("SW-82310101" = "Gwithian footbridge","SW-82310103" = "Gwithian Towans","SW-82310132" = "Kieve Bridge","SW-82310160" = "Roscroggan Bridge","SW-82310183" = "Treskillard","SW-82310186" = "Above Brea Tin works","SW-82310704" = "Coombe (Tehidy)","SW-82310411" = "Nancemellin (Roseworthy)","SW-82310471" = "Botetoe (Roseworthy)","SW-82310570" = "Praze (Roseworthy)")

#Carnon

carn_met_pca <- prcomp(na.omit(metnut_carn[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(carn_met_pca) #58,78,88 % variance

carn_nut_pca <- prcomp(na.omit(metnut_carn[,c(12:17)]), center = TRUE, scale. = TRUE)
summary(carn_nut_pca) #64,92,96 % variance

carn_pca <- carn_met_pca$x
carn_pca <- data.frame(sample = row.names(carn_pca), carn_pca)
carn_pca[c('Site', 'Time')] <- str_split_fixed(carn_pca$sample, '_', 2)

carn_npca <- carn_nut_pca$x
carn_npca <- data.frame(sample = row.names(carn_npca), carn_npca)
carn_npca[c('Site', 'Time')] <- str_split_fixed(carn_npca$sample, '_', 2)

ggbiplot(carn_met_pca, groups = carn_pca$Site, obs.scale = 1, var.scale = 1)
ggbiplot(carn_nut_pca, groups = carn_npca$Site, obs.scale = 1, var.scale = 1)
#("SW-81950522" = "Devoran Bridge","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))

#Cober

cob_met_pca <- prcomp(na.omit(metnut_cob[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(cob_met_pca) #58,78,88 % variance

carn_nut_pca <- prcomp(na.omit(metnut_carn[,c(12:17)]), center = TRUE, scale. = TRUE)
summary(carn_nut_pca) #64,92,96 % variance

carn_pca <- carn_met_pca$x
carn_pca <- data.frame(sample = row.names(carn_pca), carn_pca)
carn_pca[c('Site', 'Time')] <- str_split_fixed(carn_pca$sample, '_', 2)

carn_npca <- carn_nut_pca$x
carn_npca <- data.frame(sample = row.names(carn_npca), carn_npca)
carn_npca[c('Site', 'Time')] <- str_split_fixed(carn_npca$sample, '_', 2)

ggbiplot(carn_met_pca, groups = carn_pca$Site, obs.scale = 1, var.scale = 1)
ggbiplot(carn_nut_pca, groups = carn_npca$Site, obs.scale = 1, var.scale = 1)
#("SW-81950522" = "Devoran Bridge","SW-81950555" = "Bissoe Br","SW-81950570" = "D/S County/Wellington Adits","SW-81950622" = "Twelveheads D/S","SW-81950636" = "Twelveheads U/S","SW-81951044" = "St. Day Stream","SW-81950659" = "Chacewater Viaduct","SW-81950670" = "Blackwater","SW-81950421" = "Points Mill","SW-81950814" = "Hicks Mill","SW-81950840" = "Trehaddle","SW-81950880" = "Lanner Stw","SW-81950888" = "Gwennap Church","SW-81950890" = "Comford","SW-81950892" = "Penventon","SW-81950898" = "Coppice Gardens","SW-81950899" = "Hicks Church Row"))


#all sites

Every_site <-rbind(metnut_hayle,metnut_red,metnut_carn)
met_pca <- prcomp(na.omit(Every_site[,c(4:12)]), center = TRUE, scale. = TRUE)
summary(met_pca) #53,72,83 % variance

pca <- met_pca$x
pca <- data.frame(sample = row.names(pca), pca)
pca[c('Site', 'Time')] <- str_split_fixed(pca$sample, '_', 2)
pca$River <- substr(pca$Site, 0, 7)

ggplot(pca, aes(x = PC1, y = PC2, colour= Site)) + 
  geom_point(aes(shape= factor(River)))
                

ggbiplot(met_pca, groups = pca$Site, obs.scale = 1, var.scale = 1)



##compiled dataset, extracted from this script then manually compiled and edited in excel

comp <- read.csv("compiled.csv", header = T)

#for each metal, want to plot a line graph of conc v. distance to mouth

as <- ggplot(comp, aes(x=DTM, y=result.Arsenic_As_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Arsenic(ug/l)") + xlab("Distance to mouth (km)")

cop <- ggplot(comp, aes(x=DTM, y=result.Copper_Cu_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Copper(ug/l)") + xlab("Distance to mouth (km)")

Cad <- ggplot(comp, aes(x=DTM, y=result.Cadmium_Cd_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Cadmium(ug/l)") + xlab("Distance to mouth (km)")

iron <- ggplot(comp, aes(x=DTM, y=result.Iron_asFe_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  ylim(0,10000)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Iron(ug/l)") + xlab("Distance to mouth (km)")

Lead <- ggplot(comp, aes(x=DTM, y=result.Lead_asPb_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Lead(ug/l)") + xlab("Distance to mouth (km)")

Zinc <- ggplot(comp, aes(x=DTM, y=result.Zinc_asZn_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Zinc(ug/l)") + xlab("Distance to mouth (km)")

Nick <- ggplot(comp, aes(x=DTM, y=result.Nickel_Ni_ug.l, colour=Catchment)) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("dodgerblue2", "orange", "seagreen3","firebrick"))+
  theme_bw(base_size = 12)+
  theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Nickel(ug/l)") + xlab("Distance to mouth (km)")

library(gridExtra)
grid.arrange(cop,Cad,as,Zinc, nrow=2 )



