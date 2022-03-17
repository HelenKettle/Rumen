#makeInputDFs=function(parset){

    parNames=names(parset)
    
    if ('Ks.Sug'%in%parNames){
        SugarUsers['halfSat','Sugar']=as.numeric(parset['Ks.Sug'])
    }
    
    
    if ('Ks.H2'%in%parNames){
        MethanogensH2['halfSat','H2']=as.numeric(parset['Ks.H2'])
    }
    
    if ('Ks.AA'%in%parNames){
        AAUsers['halfSat','AminoAcid']=as.numeric(parset['Ks.AA'])
    }
    
    
    if ('GmaxMH2'%in%parNames){
        MethanogensH2['maxGrowthRate','H2']=as.numeric(parset['GmaxMH2'])
    }
    
    if ('GmaxAA'%in%parNames){
        AAUsers['maxGrowthRate','AminoAcid']=as.numeric(parset['GmaxAA'])
    }
    
    if ('GmaxSug'%in%parNames){
        SugarUsers['maxGrowthRate','Sugar']=as.numeric(parset['GmaxSug'])
    }
    
    
    if ('kLa'%in%parNames){
        kLa=as.numeric(parset['kLa'])
    }
    
    if ('washOut'%in%parNames){
        sys.res['washOut',]=as.numeric(parset['washOut'])
    }
    
    if ('Vg'%in%parNames){ #make sure this is after (not before) changing washOut values
        sys.res['washOut',c('H2.gas','CO2.gas','CH4.gas')]=rep(as.numeric(parset['Vg']),3)
    }

    if ('khyd'%in%parNames){
        khyd=as.numeric(parset['khyd'])*khyd.orig
    }
    
    
    
    if ('khyd.NSC'%in%parNames){
        khyd['NSC']=as.numeric(parset['khyd.NSC'])
    }
    if ('khyd.NDF'%in%parNames){
            khyd['NDF']=as.numeric(parset['khyd.NDF'])
        }
    if ('khyd.Protein'%in%parNames){
        khyd['Protein']=as.numeric(parset['khyd.Protein'])
    }
    
    if ('AbsorptionScale'%in%parNames){
        vfa.absorption=as.numeric(parset['AbsorptionScale'])*vfa.absorption.orig
        }
    

                                        #Note stoichiom is affected by yield for sug users and AA users
        if ('yield.Sug'%in%parNames){
            source('calc.stoichiom.yields.R')
            SugarUsers['yield','Sugar']=as.numeric(parset['yield.Sug'])
            new.stoichioms=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$stoi
            SugarUsers['Rtype','H2O']=calc.stoichiom.yields(as.numeric(SugarUsers['yield','Sugar']),'Sugar')$waterStatus
            p.names=c('Sugar','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')
            for (pn in p.names){
                SugarUsers['stoichiom',pn]=as.numeric(new.stoichioms[pn])
            }
        }

        if ('yield.AA'%in%parNames){
            source('calc.stoichiom.yields.R')
            AAUsers['yield','AminoAcid']=as.numeric(parset['yield.AA'])
            new.stoichioms=calc.stoichiom.yields(as.numeric(AAUsers['yield','AminoAcid']),'AminoAcid')$stoi
            AAUsers['Rtype','H2O']=calc.stoichiom.yields(as.numeric(AAUsers['yield','AminoAcid']),'AminoAcid')$waterStatus
            #print(new.stoichioms)
            p.names=c('AminoAcid','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')
            for (pn in p.names){
                AAUsers['stoichiom',pn]=as.numeric(new.stoichioms[pn])
            }
        }

        if ('yield.H2'%in%parNames){
            MethanogensH2['yield','H2']=as.numeric(parset['yield.H2'])
        }

#}
