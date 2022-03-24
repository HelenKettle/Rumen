#change yield of sugar users

rm(list=ls())
graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

times.h=seq(0,24,by=1/60)

spinUpTime.hours=2

yield.Sug=0.2

#change yield of sugarusers
SugarUsers['yield','Sugar']=yield.Sug
new.stoichioms=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$stoi
SugarUsers['Rtype','H2O']=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$waterStatus
p.names=c('Sugar','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')
for (pn in p.names){
    SugarUsers['stoichiom',pn]=as.numeric(new.stoichioms[pn])
}


out=rumenModel(
    times.h,
    spinUpTime.hours=spinUpTime.hours,
    additive='Control',
    basal.diet='Mixed',
    dietCompositionMat=dietCompositionMat,
    feedMat=feedMat,
    gasMat=gasMat,
    resourceSysInfo=sys.res,
    microbeSysInfo=sys.bac,
    microbeNames=c('SugarUsers','AAUsers','MethanogensH2'),
    useNetworkFuncs=FALSE,
    init.with.CH4.data=FALSE
)
 

quickPlot(out)

plotCompareGas(out,spinUpTime.hours,gasMat)


    

