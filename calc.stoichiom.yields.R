calc.stoichiom.yields=function(Yx.gg,substrate){

    #Yx.gg is the number of grams of biomass per gram of substrate (g X/ g S)
    #based on the mass balance of the stoichiometry
    
    #Depending on the value of the Yx.gg, water will end up being a substrate or a product

    if (Yx.gg>1){
        stop('Yx.gg can not be more than 1')
    }
    if (Yx.gg<0){
        stop('Yx.gg can not be negative')
    }
    molarMass=c('Sugar'=180.16,'Acetate'=60.05,'Butyrate'=88.1,'Propionate'=74.08,'H2'=2,
                    'NH3'=17,'SIC'=44,'H2O'=18,'Biomass'=113,'AminoAcid'=134)

    
    if (substrate=='Sugar'){
        
        lambda=c(0.43,0.29,0.28) #Munoz Tanayo 2016 (estimated)

                                        #Glucose utilization  

                                        #  C6H12O6 + 2H2O  -> 2CH3COOH +2CO2 + 4H2                   (R1)
                                        #  3C6H12O6        -> 2CH3COOH  + 4CH3CH2COOH + 2CO2 + 2H2O  (R2)
                                        #  C6H12O6         -> CH3CH2CH2COOH + 2CO2 + 2H2             (R3)
                                        #  5C6H12O6 + 6NH3 -> 6C5H7O2N + 18H2O                       (R4)
           
       #convert yield to mol/mol
        Yx.sug=Yx.gg*180/113 #mol X /mol S
        #Yx.sug=0.16 #mol X/mol Sug    
    
        fsu.x=(5/6)*Yx.sug # mol fraction of sug taken up that is used for growth (5 mol of sugar -> 6 mol biomass)
        fsu=max(0,1-fsu.x) #fraction of sug used for products
        #print(fsu)
        
        Y.stoi=NA*seq(1,9)
        names(Y.stoi)=c('Sugar','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')
        
        Y.stoi['Sugar']=1
        Y.stoi['Acetate']=fsu*(2*lambda[1]+(2/3)*lambda[2])
        Y.stoi['Butyrate']=fsu*lambda[3]
        Y.stoi['Propionate']=fsu*(4/3)*lambda[2]
        Y.stoi['H2']=fsu*(4*lambda[1]+2*lambda[3])
        Y.stoi['NH3']=Yx.sug
        Y.stoi['SIC']=fsu*(2*lambda[1]+(2/3)*lambda[2]+2*lambda[3])
        Y.stoi['Biomass']=Yx.sug

        #find the yield of water from mass balance
        #Y.stoi['H2O']=fsu*(-2*lambda[1]+(2/3)*lambda[2]+Yx.sug*(18/6)) #(from Rafael's code)


        uptake.mass=sum(Y.stoi[c('Sugar','NH3')]*molarMass[c('Sugar','NH3')])
        production.mass=sum(Y.stoi[c('Acetate','Butyrate','Propionate','H2','SIC','Biomass')]*molarMass[c('Acetate','Butyrate','Propionate','H2','SIC','Biomass')])
    }


    if (substrate=='AminoAcid'){
        
            #convert yield to mol/mol
        Yx.aa=Yx.gg*134/113 #mol X /mol AA
   
        #Yac_aa = (1-Yaa)*0.67;
        #Ypr_aa = (1-Yaa)*0.062;
        #Ybu_aa = (1-Yaa)*0.24;
        #Yh2_aa = (1-Yaa)*0.82;
        #YIC_aa = (1-Yaa)*0.88;

        #YIN_aa = N_aa -Yaa*N_mb;

        Y.stoi=NA*seq(1,9)
        names(Y.stoi)=c('AminoAcid','Acetate','Butyrate','Propionate','H2','NH3','SIC','H2O','Biomass')

        N.aa=1.5
        N.mb=1
        Y.stoi['AminoAcid']=1
        Y.stoi['Acetate']=(1-Yx.aa)*0.67
        Y.stoi['Butyrate']=(1-Yx.aa)*0.24
        Y.stoi['Propionate']=(1-Yx.aa)*0.062
        Y.stoi['H2']=(1-Yx.aa)*0.82
        Y.stoi['SIC']= (1-Yx.aa)*0.88
        Y.stoi['Biomass']=Yx.aa
        Y.stoi['NH3']=N.aa -Yx.aa*N.mb

        uptake.mass=sum(Y.stoi['AminoAcid']*molarMass['AminoAcid'])
        production.mass=sum(Y.stoi[c('Acetate','Butyrate','Propionate','H2','SIC','Biomass','NH3')]*molarMass[c('Acetate','Butyrate','Propionate','H2','SIC','Biomass','NH3')])
    }

    #print(paste('uptake.mass',uptake.mass))
    #print(paste('production.mass',production.mass))

    if (uptake.mass>production.mass){
        waterStatus='P'
        Y.stoi['H2O']=(uptake.mass-production.mass)/molarMass['H2O']
    }else{
        waterStatus='Sw'
        Y.stoi['H2O']=(-uptake.mass+production.mass)/molarMass['H2O']
    }
    
    return(list(stoi=Y.stoi ,waterStatus=waterStatus))

}
