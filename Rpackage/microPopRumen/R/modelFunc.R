modelFunc=function(
    times.h,
    microbeNames,
    resourceSysInfo,
    microbeSysInfo,
    rateFuncs=rateFuncsDefault,
    additive='Control',
    basal.diet='Mixed',
    dataFolder=NULL,
    spinUpTime.hours=0,
    useNetworkFuncs=FALSE,
    feedMat=NULL,
    gasMat=NULL,
    init.with.CH4.data=FALSE,
    dietCompositionMat=NULL,
    paramList=list(

        khyd.vec =  c('NDF'=0.05,'NSC'=0.20,'Protein'=0.22), #/h hydrolysis of polymers
#khyd.orig  =  c('NDF'=1.0,'NSC'=1.0,'Protein'=1.0) #/h #NOTE after GSA have changed khyd.NSC to 0.1
        kLa  = 8.33, ##Liquid-gas transfer coef (/h)(Batstone 200 /d = 8.33 /h) #1.07 
        vfa.absorption = c('Acetate'=0.33,'Propionate'=0.51,'Butyrate'=0.46), #/h VFA absorption through rumen wall from dijkstra et al 1993
        kd  =  8.3333e-4, #death rate of microbes
        fch.x  = 0.20, #mass fraction of carbs in biomass
        fpro.x  = 0.55, #mass fraction of protein in biomass
        Z0 =  0.14, #used to convert SIC to soluble CO2 (Sco2)
# Physicochemical parameters --------------------------------
        Vol.l = 100, # Volume in the liquid phase, L (Berends cite Taweel et al 2006)
        Vol.g = 40,  # Volume in the gas phase, L (Berends - no ref)
        Ptot  =  1.01325, # System pressure, bar. (approx same at atm)
        T.rumen  =  39+273.15 ,	# Temperature, K
                                        # Henrys constant  M/bar, at T = 25C (298.15K)
        KH.co2.s  = 0.035,
        KH.ch4.s  =  0.0014,
        KH.h2.s  =  7.8e-4,
        R  =  8.314*1e-2,   # bar*L/(mol*K)
    # Equilibrium constants 
        deltaH0.Ka.w  =  55900,
        deltaH0.Ka.co2  =  7646,
        deltaH0.Ka.nh4  =  51965
    )


){

   # returns list containing 'solution' (matrix from ODE solver), parms (microPop parameters),myPars (rumen parameters)))

    myPars=paramList

    myPars[['polymer.names']]=c('NDF','NSC','Protein')
    myPars[['vfa.names']]=c('Acetate','Butyrate','Propionate')
    myPars[['gas.names']]=c('H2'='H2.gas','SIC'='CO2.gas','CH4'='CH4.gas')

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

    
    if (!is.null(dietCompositionMat)){
        myPars[['polymer.frac.gPkg']]=getPolymerFrac(basal.diet,additive,dietCompositionMat)
    }else{
        stop('please provide dietCompositionMat for calculating polymer ratios')
    }
        
    
#    timeIntHours=10/60 #1 minute

    #add in variables that aren't microbial resources
    DF1= get(microbeNames[1])
    DF=cbind(DF1,'NDF'=c('X',rep(NA,6)),'NSC'=c('X',rep(NA,6)),'Protein'=c('X',rep(NA,6)))
    DF=cbind(DF,'H2.gas'=c('X',rep(NA,6)))
    DF=cbind(DF,'CO2.gas'=c('X',rep(NA,6)))#}
    DF=cbind(DF,'CH4.gas'=c('X',rep(NA,6)))#}
    assign(microbeNames[1],DF,envir = .GlobalEnv)
    


    #interpolate feed data to times.h and add spin up time
    if (!is.null(feedMat)){
        print('here')
        myPars$useFeedData=TRUE
        tstep=mean(diff(feedMat[,1]))
        spin.up.start=times.h[1]-spinUpTime.hours
        new.times.h=seq(spin.up.start,max(times.h),by=tstep)
        DMIR.new=approx(feedMat[,1],feedMat[,2],new.times.h,rule=2)$y
        TSmat=cbind('Time'=new.times.h-min(new.times.h),'DMIR'=DMIR.new)
        myPars[['TSmat']]=TSmat
    }else{
        myPars$useFeedData=TRUE
    }


    if (init.with.CH4.data & !is.null(gasMat)){
    #set initial value of CH4.gas based on mean of data
        sys.res['startValue','CH4.gas'] = mean(gasMat[,'CH4(mol/h)'])/sys.res['washOut','CH4.gas']
    }else{
        sys.res['startValue','CH4.gas']= 0
    }

    
   print('starting microPopModel()')
    
    out=microPopModel(
        microbeNames=microbeNames,
        times=times.h,
        resourceSysInfo=sys.res,
        microbeSysInfo=sys.bac,
        pHLimit=FALSE,
        rateFuncs=myRateFuncs,
        #odeFunc=derivsDefault,
        checkingOptions=list(balanceTol=1e-2,reBalanceStoichiom=FALSE,checkMassConv=FALSE,checkStoichiomBalance=TRUE,stoiTol=1),
        plotOptions=list(plotFig=FALSE),
        #odeOptions=list('atol'=1e-8,'rtol'=1e-8,'method'='lsoda'),
        #numStrains=numStrains,
        #strainOptions=strainOptions,
        networkAnalysis=useNetworkFuncs,
        myPars=myPars

    )

    print('after')

    print(names(out))

    if (useNetworkFuncs){
        return.list=list(solution=out$solution,
            flow.uptake=out$flow.uptake,flow.production=out$flow.production,
            parms=out$parms,myPars=myPars)
    }else{
        return.list=list(solution=out$solution,parms=out$parms,myPars=myPars)
    }
        

    return(return.list)

    
}
