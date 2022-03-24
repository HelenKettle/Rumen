removalRateFuncRumen=function(varName,varValue,stateVarValues,time,washOut,parms){
  #removal rate from rumen

#print('removal rate')
#print(varName)

  mw=parms$molarMass
  allSubConc=stateVarValues[parms$resourceNames]

    hydrolysis=0;     death=0;     gasTransfer=0;  absorption=0

  if (varName%in%parms$myPars$polymer.names){
      hydrolysis=parms$myPars$khyd[varName]
  }#hydrolysis of polymers

  gname=microPop::getGroupName(varName,parms$microbeNames)
  if (gname%in%parms$microbeNames){death=parms$myPars$kd}#death of microbes


  #compute gas transfer from liquid to headspace for SIC, H2 and CH4
  #for SIC compute soluble concentration using Sco2FromSIC()
  #Compute partial pressure in headspace (gas.pp in bar)
  #Compute transfer across interface using gasTransferRateFunc (g/L/h)
  
  if (varName%in%c('SIC','H2','CH4')){#dissolved gases transfer to gas phase
      
      if (varName=='SIC'){
          soluble.conc=Sco2FromSIC(allSubConc,mw,
              parms$myPars$vfa.names,parms$myPars$Z0,parms$myPars)[1]
      }else{
          soluble.conc=varValue
      }

      gas.name=parms$myPars$gas.names[varName]
      
      gas.pp=partialPressure(gas.name,parms$myPars$gas.names,
          allSubConc,mw,parms$myPars$Ptot) #bar
    
      gasTransfer=mw[varName]*gasTransferRateFunc(soluble.conc/mw[varName],
          gas.pp,parms$myPars$KH[gas.name],parms$myPars$kLa) #g/L/h

  }
  
    
    if (varName%in%parms$myPars$vfa.names){
        absorption=parms$myPars$vfa.absorption[varName]
    }
    
  v=(washOut[varName] + hydrolysis + death + absorption)*varValue + gasTransfer 

    #print('removal rate')
    #print(v)
  return(v)

}


entryRateFuncRumen=function(varName,varValue,stateVarValues,time,inflowRate,parms){


  #resource (or microbial strain mass) per unit time
    mw=parms$molarMass

    allMicrobeConc=stateVarValues[parms$allStrainNames]

  #entry rate from outside the system----------------------------------------------
    gname=microPop::getGroupName(varName,parms$microbeNames)

    if (gname%in%parms$microbeNames){
        v.in=inflowRate[gname]/parms$numStrains
        
    }else{
        
        if (varName%in%parms$myPars$polymer.names){
            if (parms$myPars$useFeedData){
                v.in=parms$myPars$polymer.frac.gPkg[varName]*
                    stats::approx(parms$myPars[['TSmat']][,'Time'],parms$myPars[['TSmat']][,'DMIR'],time,rule=2,ties=mean)$y/parms$myPars$Vol.l #g/l/h

            }else{
                v.in=inflowRate[varName]
            }
        }else{
            v.in=inflowRate[varName]
        }
    }


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
    }



    if (varName%in%parms$myPars$gas.names){#gas transferred from dissolved gases

        soluble.name=names(parms$myPars$gas.names[parms$myPars$gas.names==varName])

 
        if (soluble.name=='SIC'){
          
            v=Sco2FromSIC(stateVarValues,mw,parms$myPars$vfa.names,parms$myPars$Z0,parms$myPars)
            soluble.conc=v[1] #g/L

        }else{
            soluble.conc=stateVarValues[soluble.name] #g/L
        }
      

        pp=partialPressure(varName,parms$myPars$gas.names,stateVarValues,mw,parms$myPars$Ptot) #bar

        gasTransfer=parms$myPars$Vol.l*gasTransferRateFunc(soluble.conc/mw[soluble.name],pp,parms$myPars$KH[varName],parms$myPars$kLa) #mol/h

      

    }

    
    v=v.in+hydrolysis+input.from.dead.cells+gasTransfer
    
    return(v)
}

