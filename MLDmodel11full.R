
## Set working directory in supercomputer or local computer ###
dir5="/home/ajith/chandanp/BD-120K/BD2005sps/"
setwd(dir5)

## Get data input files structures
statespace <- read.csv("new_mask_file.csv")
traps <- read.csv("Traps.csv")
captures <- read.csv("Capture.csv")
sex <- read.csv("Sex.csv")
#effort <- read.csv("maralionEffort1.csv")

## Extract specific information from the input files ##
Xsex <- sex[,2]
#Xeffort <- effort[,4:ncol(effort)]

## For a big file such as this effort file, convert any characters entered into numeric. But be careful if there were other errors (eg: #NUM! error).##
## 'sapply' converts these columns into arbitrary integers ##

#Xeffort <- sapply(Xeffort, as.numeric) 
#alive=matrix(1,nrow=length(unique(captures[,"ANIMAL_ID"])),ncol=ncol(Xeffort))

## Get required SCR functions from the directory ###
source("e2dist.R")
source("SCRi.fn.par1-cheetah_sex.R")
source("scrDataWOeffort.R")

## Format data of class scrobj ##

scrMaraLionData <- scrData(traps=traps, captures=captures, statespace=statespace, Xsex=Xsex, Xeff=Xeffort)

### Run analysis to estimate parameters

niter <- 120000

### Number of chains to be run ###
nchains <- 3

## Set model number ##
modelno <- "BD2005sps"

## Run the model ##

scrMaraLionAnal <- SCRi.fn.par1(scrMaraLionData, modelno=modelno, nc=nchains, ni = niter, 
                                burn = 20000, skip = 1, nz = 100,theta=1,Msigma = 1, Mb = 0, 
                                Msex=1, Msexsigma = 1, Xsex = Xsex, 
                                ss.prob=NULL, coord.scale = 1000, area.per.pixel = 0.336, 
                                thinstatespace = 1, maxNN = 40, dumprate = 1000)   


