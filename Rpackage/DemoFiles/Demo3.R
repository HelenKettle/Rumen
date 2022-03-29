#change start values of microbes

rm(list=ls())
graphics.off()

library(microPop)
library(devtools)

load_all('microPopRumen')

times.h=seq(0,24,by=1/60)

spinUpTime.hours=2

#change initial conditions-------------------------------------
#microbes
sys.bac['startValue',c('SugarUsers','AAUsers','MethanogensH2')]=c(10,5,1)


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
dev.copy2pdf(file='microPopRumen/vignettes/demo3a.pdf')

plotCompareGas(out,spinUpTime.hours,gasMat)
dev.copy2pdf(file='microPopRumen/vignettes/demo3b.pdf')


    

