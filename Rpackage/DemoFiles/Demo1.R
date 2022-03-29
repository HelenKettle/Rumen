#use all defaults

rm(list=ls())
graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

times.h=seq(0,24,by=1/60)

spinUp=1

#RATE FUNCS
out=rumenModel(
    times.h,
    spinUpTime.hours=spinUp,
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
dev.copy2pdf(file='microPopRumen/vignettes/demo1a.pdf')

plotCompareGas(out,spinUp,gasMat)

dev.copy2pdf(file='microPopRumen/vignettes/demo1b.pdf')

    

