#change rate of hydrolysis

rm(list=ls())
#graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

times.h=seq(0,24,by=1/60)

spinUpTime.hours=2

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
    init.with.CH4.data=FALSE,
    paramList=list('khyd.scale'=16)#change value of scale on hydrolysis vector
)
 

quickPlot(out)

plotCompareGas(out,spinUpTime.hours,gasMat)



    

