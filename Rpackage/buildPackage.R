rm(list=ls())

library(devtools)
#library(microPop)

#load in latest microPop
#orig='/home/helen/MicroPop/Package1.5.1/microPop/'
#load_all(orig)


#create('microPopGut')

devtools::document('microPopRumen')

#test()
#test_file('tests/testthat/test_inflowFuncs.R')
#test_file('tests/testthat/testpH.R')
#test_file('tests/testthat/test_removalRate.R')
#test_file('tests/testthat/test_entryRate.R')
#test_file('tests/testthat/test_combinePaths.R')
#test_file('tests/testthat/test_GrowthLim.R')

#devtools::test('microPopGut')

#devtools::install('microPopGut')
#install('../MicroPopPackageDec21/microPop', repos = NULL, type="source")

devtools::build('microPopRumen',vignettes=FALSE) #builds a package file from package sources. You can can use it to build a binary version of your package.

devtools::build('microPopRumen')

devtools::check('microPopRumen')#,vignettes=FALSE) 
#pf='/home/helen/MicroPop/MicroPopGutPackage/microPopGut'
#load_all(pf)


#from file
#install.packages('/home/helen/MicroPop/MicroPopPackageDec21/microPop',repos=NULL,type='source')
#library(microPop)
#help(package='microPop')

#from CRAN
#install.packages('microPop')
#library(microPop)
#help(package='microPop')



#test(pf)

#test_coverage(pf) #brings up results in firefox!

runDemos=TRUE
if (runDemos){
    demoNames=c('human1','human2','human3','human4','phyto','rumen','strains','p
hages')
    for (demo in demoNames){
        source(paste0(pf,'inst/DemoFiles/',demo,'.R'))
    }
}


install(pf)

library(microPop)

#help(package='microPop')

plot(c(0,10),c(0,4),type='n',xlab='time',ylab='B')

myRateFuncs=rateFuncsDefault

myRateFuncs$removalRateFunc=function(varName, varValue, stateVarValues, time,
                                  washOut, parms) {
    

    if (parms$myPars[['compartment']]==1){
        V=3
    }
    
        if (varValue <= 0) {
            v = 0
        } else {
            v = V * varValue
        }
    
    
    return(v)
}

out=microPopModel(
    microbeNames='Bacteroides',
    times=seq(0,10,0.1),
    rateFuncs=myRateFuncs,
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeSysInfoHuman,
    plotOptions=list(plotFig=FALSE),
    myPars=c('compartment'=1)
)
   
lines(out$solution[,1],out$solution[,2],type='l',xlab='time',ylab='B')
