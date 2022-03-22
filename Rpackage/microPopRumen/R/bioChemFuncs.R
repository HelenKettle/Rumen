solverh=function(Z,Svfa,SIC,SIN,myPars){
    #Function for determination of H+ (mol) i.e. ionH used when converting SIC to Sco2
    #H+ is also known as proton or hydron
    #(Munoz-Tamayo et al 2016)
    #requires: K.a.co2  K.a.nh4  K.a.vfa K.w
    #SIC in moles (soluble inorganic carbon)
    #SIN in moles (soluble inorganic nitrogen)
    #Svfa soluble volatile fatty acids (moles)

    C=NA*seq(1,6)
    C[1] = 1
    C[2] = (myPars$K.a.co2 + myPars$K.a.nh4 + myPars$K.a.vfa + Z + SIN)
    C[3] = (myPars$K.a.co2*Z - myPars$K.w - myPars$K.a.co2*SIC - myPars$K.a.vfa*Svfa + myPars$K.a.nh4*(myPars$K.a.co2 + Z) +
            SIN*(myPars$K.a.co2 + myPars$K.a.vfa) + myPars$K.a.vfa*(myPars$K.a.co2 + myPars$K.a.nh4 + Z))
    C[4]= (myPars$K.a.vfa*(myPars$K.a.co2*Z + myPars$K.a.nh4*(myPars$K.a.co2 + Z)) - myPars$K.w*(myPars$K.a.co2 + myPars$K.a.nh4 + myPars$K.a.vfa) -
           myPars$K.a.co2*SIC*(myPars$K.a.nh4 + myPars$K.a.vfa) - myPars$K.a.vfa*Svfa*(myPars$K.a.co2 + myPars$K.a.nh4) +
           myPars$K.a.co2*myPars$K.a.nh4*Z + myPars$K.a.co2*myPars$K.a.vfa*SIN)
    C[5] = (myPars$K.a.co2*myPars$K.a.nh4*myPars$K.a.vfa*Z - myPars$K.w*(myPars$K.a.vfa*(myPars$K.a.co2 + myPars$K.a.nh4) +myPars$K.a.co2*myPars$K.a.nh4) -
            myPars$K.a.co2*myPars$K.a.nh4*myPars$K.a.vfa*SIC -myPars$K.a.co2*myPars$K.a.nh4*myPars$K.a.vfa*Svfa)
    C[6] = - myPars$K.a.co2*myPars$K.a.nh4*myPars$K.a.vfa*myPars$K.w
    Hroots = rev(polyroot(as.vector(rev(C))))
    Proots = Re(Hroots[abs(Im(Hroots))<1e-15 & Re(Hroots)>=0])
    H = min(max(Proots)[1],10^(-4))

    return(H)

}


partialPressure=function(gas.name,all.gas.names,allSubConc,mw,Ptot){
    tiny=1e-16
    #returns partial pressure (molar fraction of gases)*Ptot :. same units as Ptot
    #all.gas.names is a vector of strings eg c('Gh2','Gco2','Gch4')
    all.g.mol=allSubConc[all.gas.names]
    gas.pp=Ptot*all.g.mol[gas.name]/(sum(all.g.mol,na.rm=TRUE)+tiny)
    return(gas.pp)
}


gasTransferRateFunc=function(soluble.conc,gas.pp,K.henry,kLa){
    #returns mass transfer in moles/L/h
    #K.henry in M/bar (=mol/L/bar) at rumen temp
    #soluble.conc is in mol/L
    #gas.pp is in bar
    #kLa is gas transfer coef (/h)
    return(kLa*(soluble.conc-K.henry*gas.pp))
}

Sco2FromSIC=function(allSubConc,mw,vfa.names,Z0,myPars){
  #Gives Sco2 in g/L (and pH)
    #uses solverh to compute H+ i.e. ionH
    #allSubConc (g/L) 
    #mw molar mass
    SIC=allSubConc['SIC']/mw['SIC']
    SIN=allSubConc['NH3']/mw['NH3']
    Svfa = sum(allSubConc[vfa.names]/mw[vfa.names],na.rm=TRUE)
    ionH = solverh(Z0,Svfa,SIC,SIN,myPars)
    #print('ionH')
    #print(ionH)
    soluble.conc = allSubConc['SIC']*(1- myPars$K.a.co2/(myPars$K.a.co2+ionH))
    pH=-log10(ionH)

    return(c(max(soluble.conc,0),pH))
}

