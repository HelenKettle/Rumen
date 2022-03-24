#change transit time

rm(list=ls())
graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

times.h=seq(0,24,by=1/60)

spinUpTime.hours=2

transitTime.h=12

sys.res['washOut',]=1/transitTime.h
sys.bac['washOut',]=1/transitTime.h


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


    

