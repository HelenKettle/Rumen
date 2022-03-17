solverh=function(Z,Svfa,SIC,SIN){
    #Function for determination of H+ (mol) i.e. ionH used when converting SIC to Sco2
    #H+ is also known as proton or hydron
    #(Munoz-Tamayo et al 2016)
    #global K.a.co2  K.a.nh4  K.a.vfa K.w
    #SIC in moles (soluble inorganic carbon)
    #SIN in moles (soluble inorganic nitrogen)
    #Svfa soluble volatile fatty acids (moles)
    C=NA*seq(1,6)
    C[1] = 1
    C[2] = (K.a.co2 + K.a.nh4 + K.a.vfa + Z + SIN)
    C[3] = (K.a.co2*Z - K.w - K.a.co2*SIC - K.a.vfa*Svfa + K.a.nh4*(K.a.co2 + Z) +
            SIN*(K.a.co2 + K.a.vfa) + K.a.vfa*(K.a.co2 + K.a.nh4 + Z))
    C[4]= (K.a.vfa*(K.a.co2*Z + K.a.nh4*(K.a.co2 + Z)) - K.w*(K.a.co2 + K.a.nh4 + K.a.vfa) -
           K.a.co2*SIC*(K.a.nh4 + K.a.vfa) - K.a.vfa*Svfa*(K.a.co2 + K.a.nh4) +
           K.a.co2*K.a.nh4*Z + K.a.co2*K.a.vfa*SIN)
    C[5] = (K.a.co2*K.a.nh4*K.a.vfa*Z - K.w*(K.a.vfa*(K.a.co2 + K.a.nh4) +K.a.co2*K.a.nh4) -
            K.a.co2*K.a.nh4*K.a.vfa*SIC -K.a.co2*K.a.nh4*K.a.vfa*Svfa)
    C[6] = - K.a.co2*K.a.nh4*K.a.vfa*K.w
    Hroots = rev(polyroot(as.vector(rev(C))))
    Proots = Re(Hroots[abs(Im(Hroots))<1e-15 & Re(Hroots)>=0])
    H = min(max(Proots)[1],10^(-4));
    return(H)
}

gasTransferRateFunc=function(soluble.conc,gas.pp,K.henry,kLa){
    #returns mass transfer in moles/L/h
    #K.henry in M/bar (=mol/L/bar) at rumen temp
    #soluble.conc is in mol/L
    #gas.pp is in bar
    #kLa is gas transfer coef (/h)
    return(kLa*(soluble.conc-K.henry*gas.pp))
}

partialPressure=function(gas.name,all.gas.names,allSubConc,mw,Ptot){
    #returns partial pressure (molar fraction of gases)*Ptot :. same units as Ptot
    #all.gas.names is a vector of strings eg c('Gh2','Gco2','Gch4')
    all.g.mol=allSubConc[all.gas.names]
    return(Ptot*all.g.mol[gas.name]/(sum(all.g.mol,na.rm=TRUE)+eps))
}

Sco2FromSIC=function(allSubConc,mw,vfa.names,Z0){
  #Gives Sco2 in g/L (and pH)
    #uses solverh to compute H+ i.e. ionH
    #allSubConc (g/L) 
    #mw molar mass
    Svfa = sum(allSubConc[vfa.names]/mw[vfa.names],na.rm=TRUE)
    ionH = solverh(Z0,Svfa,allSubConc['SIC']/mw['SIC'],allSubConc['NH3']/mw['NH3'])
    soluble.conc = allSubConc['SIC']*(1- K.a.co2/(K.a.co2+ionH))
    pH=-log10(ionH)
    return(c(max(soluble.conc,0),pH))
}

