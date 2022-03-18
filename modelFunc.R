modelFunc=function(
    microbeNames,
    resourceSysInfo,
    microbeSysInfo,
    rateFuncs=rateFuncsDefault,
    ID='001',
    additive='Control',
    basal.diet='Mixed',
    params=parset,
    dataFolder=NULL,
    spinUpTime.hours=0,
    useNetworkFuncs=FALSE,
    useFeedData=FALSE,
    useGasData=FALSE,
    paramList
){

   # returns list containing 'solution' (matrix from ODE solver), parms (microPop parameters),myPars (rumen parameters)))

    myPars=paramList


    myPars[['f.X']] = c('NSC'=paramList$fch.x,'Protein'=paramList$fpro.x)
    
    #convert to rumen temperature
    KH.co2  =  paramList$KH.co2.s*exp(-(19410/(paramList$R*100))*(1/298.15-1/paramList$T.rumen))
    KH.ch4  =  paramList$KH.ch4.s*exp(-(14240/(paramList$R*100))*(1/298.15-1/paramList$T.rumen))
    KH.h2   =  paramList$KH.h2.s*exp(-(4180/(paramList$R*100))*(1/298.15-1/paramList$T.rumen))
    

    myPars[['KH']]=c('H2.gas'=KH.h2,'CO2.gas'=KH.co2,'CH4.gas'=KH.ch4)
    
# Acid-base constants mol/L
    myPars[['K.w']]  =  exp(paramList$deltaH0.Ka.w/(paramList$R*100)*(1/298.15-1/paramList$T.rumen))*1e-14
    myPars[['K.a.co2']]  =  10^(-6.35)*exp(paramList$deltaH0.Ka.co2/(paramList$R*100)*(1/298.15-1/paramList$T.rumen))
    myPars[['K.a.nh4']]  =   10^(-9.25)*exp(paramList$deltaH0.Ka.nh4/(paramList$R*100)*(1/298.15-1/paramList$T.rumen))
    myPars[['K.a.ac']]  =  10^(-4.76)
    myPars[['K.a.bu']] = 10^(-4.82)
    myPars[['K.a.pr']] = 10^(-4.88)
    myPars[['K.a.vfa']] = 10^(-4.76)#K.a.ac

    myPars[['polymer.frac.gPkg']]=getPolymerFrac(basal.diet,additive,dataFolder)
    
    timeIntHours=1/60 #1 minute

    #add in variables that aren't resources
    DF1= get(microbeNames[1])
    DF=cbind(DF1,'NDF'=c('X',rep(NA,6)),'NSC'=c('X',rep(NA,6)),'Protein'=c('X',rep(NA,6)))
   # if ('H2'%in%resNames){
        DF=cbind(DF,'H2.gas'=c('X',rep(NA,6)))
#}
   # if ('SIC'%in%resNames){
        DF=cbind(DF,'CO2.gas'=c('X',rep(NA,6)))#}
 #   if ('CH4'%in%resNames){
         DF=cbind(DF,'CH4.gas'=c('X',rep(NA,6)))#}
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
    myPars[['TSmat']]=loadDataFunc(ID,dataFolder,spinUpTime.hours,useFeedData,useGasData)


    #CHANGE TO SHORTER TIME?
    times=myPars[['TSmat']][1:1000,1]
    


    if (useGasData){
    #set initial value of CH4.gas based on mean of data
        sys.res['startValue','CH4.gas'] = mean(myPars[['TSmat']][,'MPRmPh'])/sys.res['washOut','CH4.gas']
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
        networkAnalysis=useNetworkFuncs,
        myPars=myPars

    )

    print('after')


    #save(out,ID,pH,Sco2.keep,ode.times,TSmat,sys.res,sys.bac,parset,file=new.filename)

    return(list(solution=out$solution,parms=out$parms,myPars=myPars))

    
}
