paramList=list(

    polymer.names=c('NDF','NSC','Protein'),
    vfa.names=c('Acetate','Butyrate','Propionate'),
    gas.names=c('H2'='H2.gas','SIC'='CO2.gas','CH4'='CH4.gas'),


                                        #hydrolysis of polymers
    khyd =  c('NDF'=0.05,'NSC'=0.20,'Protein'=0.22), #/h #proper orig I think
#khyd.orig  =  c('NDF'=1.0,'NSC'=1.0,'Protein'=1.0) #/h #NOTE after GSA have changed khyd.NSC to 0.1

#Liquid-gas transfer coef (/h)
    kLa  = 8.33, #(Batstone 200 /d = 8.33 /h) #1.07 

#VFA absorption through rumen wall
#vfa.absorption.orig  =  0.5*c('Acetate'=0.33,'Propionate'=0.51,'Butyrate'=0.46) #/h from dijkstra et al 1993
    vfa.absorption = c('Acetate'=0.33,'Propionate'=0.51,'Butyrate'=0.46), #/h from dijkstra et al 1993


#Fixed Parameters----------------------------------------------------------------------------------------------------------------

    eps = 2.2e-16, #tiny number to avoid dividing by zero

    kd  =  8.3333e-4, #death rate of microbes

#mass fraction of carb and pro in biomass
    fch.x  = 0.20,
    fpro.x  = 0.55,

    Z0 =  0.14, #used to convert SIC to soluble CO2 (Sco2)

#max specific utilisation rate constant for sugar and amino acid (for reference only).
#Note, Gmax=km*yield (mol/mol). #km.su    = 0.99; Ks.su    = 9e-3;   mol/mol/h

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


