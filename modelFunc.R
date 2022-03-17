modelFunc=function(
    microbeNames,
    resourceSysInfo,
    microbeSysInfo,
    rateFuncs=rateFuncsDefault,
    ID='001',
    additive='Control',
    basal.diet='Mixed',
    params=parset,
    resultsFilepath=NULL,
    dataFolder=NULL,
    spinUpTime.hours=0,
    useNetworkFuncs=FALSE,
    useFeedData=FALSE,
    useGasData=FALSE
){

    parNames=names(parset)

    #REQUIRES: (these are sourced later in this file)
    #bioChemFuncs.R
    #parameters.R
    #FuncsForMicroPop.R

    new.filename=paste(resultsFilepath,'run',ID,'.Rdata',sep='')
    print('filename for saving')
    print(new.filename)

    polymer.frac.gPkg<<-getPolymerFrac(basal.diet,additive,dataFolder)
    
    print(polymer.frac.gPkg)
    timeIntHours=1/60 #1 minute

    polymer.names<<-c('NDF','NSC','Protein')
    vfa.names<<-c('Acetate','Butyrate','Propionate')
    gas.names<<-c('H2'='H2.gas','SIC'='CO2.gas','CH4'='CH4.gas')
    

    resNames<<-getAllResources(microbeNames)  

    khyd<<-khyd.orig
    vfa.absorption<<-vfa.absorption.orig


    #add in variables that aren't resources
    DF1= get(microbeNames[1])
    DF<<-cbind(DF1,'NDF'=c('X',rep(NA,6)),'NSC'=c('X',rep(NA,6)),'Protein'=c('X',rep(NA,6)))
    if ('H2'%in%resNames){
        DF<<-cbind(DF,'H2.gas'=c('X',rep(NA,6)))}
    if ('SIC'%in%resNames){
        DF<<-cbind(DF,'CO2.gas'=c('X',rep(NA,6)))}
    if ('CH4'%in%resNames){
         DF<<-cbind(DF,'CH4.gas'=c('X',rep(NA,6)))}
    assign(microbeNames[1],DF,envir = .GlobalEnv)
    

        #print(DF)

#global variables (want to output from ODE rate functions)
    pH<<-NA
    ode.times<<-NA
    Sco2.keep<<-NA
    gasTransfer.CO2<<-NA #g/l/h
    gasTransfer.CH4<<-NA
    gasTransfer.H2<<-NA



    #load gas and feed data
    TSmat<<-loadDataFunc(ID,dataFolder,spinUpTime.hours,useFeedData,useGasData)


    #CHANGE TO SHORTER TIME?
    times=TSmat[1:200,1]
    


    if (useGasData){
    #set initial value of CH4.gas based on mean of data
        sys.res['startValue','CH4.gas']<<-mean(TSmat[,'MPRmPh'])/sys.res['washOut','CH4.gas']
    }else{
        sys.res['startValue','CH4.gas']=4.5
    }


   print('starting microPopModel()')
    
    out=microPopModel(
        microbeNames=microbeNames,
        times=times,
        resourceSysInfo=sys.res,
        microbeSysInfo=sys.bac,
        pHLimit=FALSE,
        rateFuncs=myRateFuncs,
        odeFunc=derivsDefault,
        checkingOptions=list(balanceTol=1e-2,reBalanceStoichiom=FALSE,checkMassConv=FALSE,checkStoichiomBalance=TRUE,stoiTol=1),
        plotOptions=list(yLabel='quantity',xLabel='time (h)',plotFig=TRUE),
        #odeOptions=list('atol'=1e-8,'rtol'=1e-8,'method'='lsoda'),
        #numStrains=numStrains,
        #strainOptions=strainOptions,
        networkAnalysis=useNetworkFuncs
        )

    print('after')


    #save(out,ID,pH,Sco2.keep,ode.times,TSmat,sys.res,sys.bac,parset,file=new.filename)

    return(out)

    
}
