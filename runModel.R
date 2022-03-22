#Run for one cow and use parameters from file
#record network flows
#shell arg is animal ID

#TIME IS IN HOURS THROUGHOUT

#to run for a single run but use the default values that are specified for the OAT then choose runtype=OAT and numDivsOAT=1

rm(list=ls())
#graphics.off()


library(microPop)

source('makeParamSets.R')
source('modelFunc.R')
source('getPolymerFrac.R')
source('loadDataFunc.R')
source('bioChemFuncs.R')
source('parameters.R') #makes paramList
source('FuncsForMicroPop.R') #this defines the new myRateFuncs functions
source('quickPlot.R')
source('plotCompareGas.R')

timeStep=1/60 #one minute
times.h=seq(0,48,by=timeStep)

useNetworkFuncs=FALSE

init.with.CH4.data=FALSE

spinUpTime.hours=3

additive='Control'
diet='Mixed'

#sensitivity analysis indicates that these 4 parameters are v influential
parset=c(
    'yield.Sug'=0.318,
    'washOut'=0.12,
    'khyd.scale'=4,#8,
    'Ks.H2'=1.2e-5)

parNames=names(parset)


#load feed data
#id=601199
#load(paste(dataFolder,'feedRate',id,'.RData',sep=''))
#feedMat=cbind('time(h)'=24*(feedRate[,1]-min(feedRate[,1])),'DMIR(kg/h)'=feedRate[,2]/24)
#save(feedMat,file='data/feedMat.Rdata')

#feed data in a matrix of time (column 1) and dry matter intake rate (kg/d) (col 2)
#load('data/feedMat.Rdata')

#load gas data
#id=601199
#load(paste(dataFolder,'gas',id,'.RData',sep=''))
#MPRmol.p.h=60*animalGasData[,'CH4.g.min.']/16.04
#gasMat=cbind('time(h)'=(animalGasData[,'JD']-min(animalGasData[,'JD']))*24,'CH4(mol/h)'=MPRmol.p.h)
#save(gasMat,file='data/gasMat.Rdata')


load('data/gasMat.Rdata')
load('data/feedMat.Rdata')
load('data/dietComposition.Rdata')


#input DFs----------------
#Create microbe data frames
load('data/SugarUsers.Rdata')# = createDF(paste0(CSVfiles,'/SugarUtilizers.csv'))
load('data/AAUsers.Rdata')# = createDF(paste0(CSVfiles,'/AminoAcidUsers.csv'))
load('data/MethanogensH2.Rdata')# = createDF(paste0(CSVfiles,'/Hydrogenotrophic_methanogen.csv'))

#save(SugarUsers,file='data/SugarUsers.Rdata')
#save(MethanogensH2,file='data/MethanogensH2.Rdata')
#save(AAUsers,file='data/AAUsers.Rdata')
microbeNames = c('SugarUsers','AAUsers','MethanogensH2')

#get system files
#sys.res = createDF(paste0(CSVfiles,'/ResourceSysInfoTroy.csv'))
#sys.bac = createDF(paste0(CSVfiles,'/microbesMixedTroy.csv'))

#save(sys.res,file='data/sys.res.Rdata')
#save(sys.bac,file='data/sys.bac.Rdata')

load('data/sys.res.Rdata')
load('data/sys.bac.Rdata')

#--update data frames with new parameter values------------------

#Note stoichiom is affected by yield for sug users and AA users
if ('yield.Sug'%in%parNames){
    source('calc.stoichiom.yields.R')
    SugarUsers['yield','Sugar']=as.numeric(parset['yield.Sug'])
    new.stoichioms=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$stoi
    SugarUsers['Rtype','H2O']=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$waterStatus
    p.names=c('Sugar','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')
    for (pn in p.names){
        SugarUsers['stoichiom',pn]=as.numeric(new.stoichioms[pn])
    }
}
if ('washOut'%in%parNames){
    sys.res['washOut',]=as.numeric(parset['washOut'])
}
if ('washOut'%in%parNames){
    sys.bac['washOut',]=as.numeric(parset['washOut'])
}
if ('khyd.scale'%in%parNames){
    paramList[['khyd']]=as.numeric(parset['khyd.scale'])*paramList[['khyd.vec']]
}
if ('Ks.H2'%in%parNames){
    MethanogensH2['halfSat','H2']=as.numeric(parset['Ks.H2'])
}
#-------------------------------------------------------------------


#RATE FUNCS
myRateFuncs=rateFuncsDefault
myRateFuncs$entryRateFunc=entryRateFunc
myRateFuncs$removalRateFunc=removalRateFunc


out=modelFunc(
    times.h=times.h,
    microbeNames=microbeNames,
    resourceSysInfo=sys.res,
    microbeSysInfo=sys.bac,
    rateFuncs=myRateFuncs,
    additive=additive,
    basal.diet=diet,
    dataFolder=dataFolder,
    spinUpTime.hours=spinUpTime.hours,
    useNetworkFuncs=useNetworkFuncs,
    feedMat=feedMat,
    gasMat=gasMat,
    init.with.CH4.data=init.with.CH4.data,
    dietCompositionMat=dietComposition,
    paramList
)

quickPlot(out)

plotCompareGas(out,spinUpTime.hours,gasMat)

time=out$solution[,1]
microbeNames=out$parms$microbeNames

 
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
           figWidth = 600,figHeight = 600)
    
    print(vv)
    
}

    

