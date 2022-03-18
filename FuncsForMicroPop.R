removalRateFunc=function(varName,varValue,stateVarValues,time,washOut,parms){
  #removal rate from rumen

#print('removal rate')
#print(varName)

  mw=parms$molarMass
  allSubConc=stateVarValues[parms$resourceNames]

    hydrolysis=0;     death=0;     gasTransfer=0;  absorption=0

  if (varName%in%parms$myPars$polymer.names){
      hydrolysis=parms$myPars$khyd[varName]
  }#hydrolysis of polymers

  gname=getGroupName(varName,microbeNames)
  if (gname%in%microbeNames){death=parms$myPars$kd}#death of microbes
  
  if (varName%in%c('SIC','H2','CH4')){#dissolved gases transfer to gas phase
    if (varName=='SIC'){
      soluble.conc=Sco2FromSIC(allSubConc,mw,
          parms$myPars$vfa.names,parms$myPars$Z0,parms$myPars)[1]
    }else{
      soluble.conc=varValue
    }

    #print(parms$myPars)
    
    pp=partialPressure(parms$myPars$gas.names[varName],
        parms$myPars$gas.names,allSubConc,mw,parms$myPars$Ptot) #bar

    #print('pp')
    #print(pp)

#    print('parms$myPars$KH')
 #       print(parms$myPars$KH)
    
    gasTransfer=mw[varName]*gasTransferRateFunc(soluble.conc/mw[varName],
        pp,parms$myPars$KH[parms$myPars$gas.names[varName]],parms$myPars$kLa) #g/L/h

  #  print('gasTransfer')
   # print(gasTransfer)
  }
    
    if (varName%in%parms$myPars$vfa.names){
        absorption=parms$myPars$vfa.absorption[varName]
    }
    
  v=(washOut[varName] + hydrolysis + death + absorption)*varValue + gasTransfer 

    #print('removal rate')
    #print(v)
  return(v)

}


entryRateFunc=function(varName,varValue,stateVarValues,time,inflowRate,parms){

    #print('varName')
    #print(varName)
 #   print(names(stateVarValues))


  #resource (or microbial strain mass) per unit time
    mw=parms$molarMass
  #allMicrobeConc=stateVarValues[parms$microbeNames]
    allMicrobeConc=stateVarValues[parms$allStrainNames]
    #print(names(parms))
    #print('allMicrobeConc')

  #entry rate from outside the system----------------------------------------------
    gname=getGroupName(varName,parms$microbeNames)

    if (gname%in%parms$microbeNames){
        v.in=inflowRate[gname]/parms$numStrains

    }else{
        
        if (varName%in%parms$myPars$polymer.names){
            if (useFeedData){
                v.in=parms$myPars$polymer.frac.gPkg[varName]*
                    approx(parms$myPars[['TSmat']][,'Time'],parms$myPars[['TSmat']][,'DMIR'],time,rule=2,ties=mean)$y/parms$myPars$Vol.l #g/l/h

          }else{
              v.in=inflowRate[varName]
          }
     # }else if (varName=='CH4.chamber'){
     #     v.in=Vg*stateVarValues['CH4.gas']
     
      }else{
          v.in=inflowRate[varName]
      #        print(v.in)
      }
  }

#    print('v.in')
#        print(v.in)

   #entry rate from inside the system-----------------------------------------------
  hydrolysis=0
  input.from.dead.cells=0
  gasTransfer=0

  if (varName=='Sugar'){#entry rate of sugar from hydrolysed polymers
      hydrolysis=sum(parms$myPars$khyd[c('NDF','NSC')]*stateVarValues[c('NDF','NSC')])#g/L/h

  }

    if (varName=='AminoAcid'){#entry rate of amino acids from hydrolysed proteins
        hydrolysis=parms$myPars$khyd['Protein']*stateVarValues['Protein']#g/L/h

    }
    
    if (varName=='NSC' | varName=='Protein'){#dead microbial cells become Znsc and Zpro
        input.from.dead.cells=parms$myPars$f.X[varName]*parms$myPars$kd*sum(allMicrobeConc) #g/L/h
#        print(input.from.dead.cells)
    }



    if (varName%in%parms$myPars$gas.names){#gas transferred from dissolved gases

        soluble.name=names(parms$myPars$gas.names[parms$myPars$gas.names==varName])

 
        if (soluble.name=='SIC'){
          
        v=Sco2FromSIC(stateVarValues,mw,parms$myPars$vfa.names,parms$myPars$Z0,parms$myPars)
        soluble.conc=v[1] #g/L

        if (length(ode.times)<20){
            pH <<-c(pH,v[2])
            ode.times<<-c(ode.times,time)
            Sco2.keep<<-c(Sco2.keep,soluble.conc)
        }else{
            time.diff=time-max(ode.times,na.rm=TRUE)
            if (time.diff>(10/60)){ #only save every 10 mins
                pH <<-c(pH,v[2])
                ode.times<<-c(ode.times,time)
                Sco2.keep<<-c(Sco2.keep,soluble.conc)
            }
            
        }


      }else{
        soluble.conc=stateVarValues[soluble.name] #g/L
      }
      

        pp=partialPressure(varName,parms$myPars$gas.names,stateVarValues,mw,parms$myPars$Ptot) #bar

      gasTransfer=parms$myPars$Vol.l*gasTransferRateFunc(soluble.conc/mw[soluble.name],pp,parms$myPars$KH[varName],parms$myPars$kLa) #mol/h


#      print('pp')
#       print(pp)
      
      if (length(ode.times)<20){
 
          if (varName=='H2.gas'){gasTransfer.H2<<-c(gasTransfer.H2,gasTransfer)}
          if (varName=='CO2.gas'){gasTransfer.CO2<<-c(gasTransfer.CO2,gasTransfer)}
          if (varName=='CH4.gas'){gasTransfer.CH4<<-c(gasTransfer.CH4,gasTransfer)}

      }else{
          
          time.diff=time-max(ode.times,na.rm=TRUE)
          
          if (time.diff>(10/60)){ #only save every 10 mins
              if (varName=='H2.gas'){gasTransfer.H2<<-c(gasTransfer.H2,gasTransfer)}
              if (varName=='CO2.gas'){gasTransfer.CO2<<-c(gasTransfer.CO2,gasTransfer)}
              if (varName=='CH4.gas'){gasTransfer.CH4<<-c(gasTransfer.CH4,gasTransfer)}
          }
      }
  }
    
  v=v.in+hydrolysis+input.from.dead.cells+gasTransfer

 #   print('v.in')
 #   print(v.in)
 #   print('hydrolysis')
 #   print(hydrolysis)
 #   print('input.from.dead.cells')
 #   print(input.from.dead.cells)
 #   print('gasTransfer')
 #   print(gasTransfer)
 #  print('entry rate')
 #  print(v)

    return(v)
}

