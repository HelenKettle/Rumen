#TIME IS IN HOURS THROUGHOUT

rm(list=ls())
#graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

#source('parameters.R')

timeStep=1/60 #one minute
times.h=seq(0,48,by=timeStep)

useNetworkFuncs=FALSE

init.with.CH4.data=FALSE

spinUpTime.hours=1

additive='Control'
diet='Mixed'

#sensitivity analysis indicates that these 4 parameters are v influential
parset=c(
    'yield.Sug'=0.318,
    'washOut'=0.12,
    'khyd.scale'=16,
    'Ks.H2'=1.2e-5)

parNames=names(parset)

microbeNames = c('SugarUsers','AAUsers','MethanogensH2')

#--update data frames with new parameter values------------------

#Note stoichiom is affected by yield for sug users and AA users
if ('yield.Sug'%in%parNames){
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
    paramList=paramList
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

    

