#Run for one cow and use parameters from file
#record network flows
#shell arg is animal ID

#to run for a single run but use the default values that are specified for the OAT then choose runtype=OAT and numDivsOAT=1

rm(list=ls())
graphics.off()


library(microPop)

source('makeParamSets.R')
source('modelFunc.R')
source('dataFuncs.R')
source('loadDataFunc.R')
source('bioChemFuncs.R')
source('parameters.R')
source('FuncsForMicroPop.R') #this defines the new myRateFuncs functions
source('quickPlot.R')

useFeedData=TRUE
useGasData=TRUE

useNetworkFuncs=FALSE

paramsFromFile=TRUE

spinUpTime.hours=10

dataFolder='Data/'

CSVfiles='../CSVfiles'

#load('../mechModelOATfit1.Rdata')

id='601199'

additive='Control'
diet='Mixed'

print(id)

resultsFilepath=paste('../../Results/MechModel/singleCows/',sep='')


parset=c(
    'yield.Sug'=0.318,
    'washOut'=0.12,
    'khyd'=8,
    'Ks.H2'=1.2e-5)

#input DFs----------------
#Create microbe data frames
SugarUsers = createDF(paste0(CSVfiles,'/SugarUtilizers.csv'))
AAUsers = createDF(paste0(CSVfiles,'/AminoAcidUsers.csv'))
MethanogensH2 = createDF(paste0(CSVfiles,'/Hydrogenotrophic_methanogen.csv'))
microbeNames = c('SugarUsers','AAUsers','MethanogensH2')

#get system files
sys.res = createDF(paste0(CSVfiles,'/ResourceSysInfoTroy.csv'))
sys.bac = createDF(paste0(CSVfiles,'/microbesMixedTroy.csv'))

source('makeInputDFs.R')

#RATE FUNCS
myRateFuncs=rateFuncsDefault
myRateFuncs$entryRateFunc=entryRateFunc
myRateFuncs$removalRateFunc=removalRateFunc


out=modelFunc(
    microbeNames=microbeNames,
    resourceSysInfo=sys.res,
    microbeSysInfo=sys.bac,
    rateFuncs=myRateFuncs,
    ID=id,
    additive=additive,
    basal.diet=diet,
    params=parset,
    resultsFilepath=resultsFilepath,
    dataFolder=dataFolder,
    spinUpTime.hours=spinUpTime.hours,
    useNetworkFuncs=useNetworkFuncs,
    useFeedData=useFeedData,
    useGasData=useGasData
)

quickPlot(out)

microbeNames=out$parms$microbeNames
resourceNames=out$parms$resourceNames
time=out$solution[,1]
Bac=out$solution[,microbeNames]
methane=out$solution[,'CH4.gas']
MPR=methane*out$parms$Smats$washOut['CH4.gas']
 
if (useNetworkFuncs){

    DFs=networkDFfromMPoutput(
        max(time),
        out,
        microbeNames,
        sumOverPaths=TRUE
    )
    
    vv=getVNPlotObject(DFs$nodes,DFs$edges,
        scaleNodes= TRUE,scaleEdges = TRUE,
        mainTitle=paste('Network after',round((max(time)-min(time)-spinUpTime.hours)/24,1),'d'), 
           figWidth = 500,figHeight = 500)
    
    print(vv)
    
}

    

