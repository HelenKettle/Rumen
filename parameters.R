#hydrolysis of polymers
khyd.orig <<- c('NDF'=0.05,'NSC'=0.20,'Protein'=0.22) #/h #proper orig I think
#khyd.orig <<- c('NDF'=1.0,'NSC'=1.0,'Protein'=1.0) #/h #NOTE after GSA have changed khyd.NSC to 0.1

#Liquid-gas transfer coef (/h)
kLa <<-8.33 #(Batstone 200 /d = 8.33 /h) #1.07 

#VFA absorption through rumen wall
#vfa.absorption.orig <<- 0.5*c('Acetate'=0.33,'Propionate'=0.51,'Butyrate'=0.46) #/h from dijkstra et al 1993
vfa.absorption.orig <<- c('Acetate'=0.33,'Propionate'=0.51,'Butyrate'=0.46) #/h from dijkstra et al 1993


#Fixed Parameters----------------------------------------------------------------------------------------------------------------

eps<<-2.2e-16 #tiny number to avoid dividing by zero

kd <<- 8.3333e-4 #death rate of microbes

#mass fraction of carb and pro in biomass
fch.x <<-0.20
fpro.x <<-0.55
f.X<<-c('NSC'=fch.x,'Protein'=fpro.x)

Z0<<- 0.14 #used to convert SIC to soluble CO2 (Sco2)

#max specific utilisation rate constant for sugar and amino acid (for reference only). Note, Gmax=km*yield (mol/mol). #km.su    = 0.99; Ks.su    = 9e-3;   mol/mol/h

# Physicochemical parameters --------------------------------
Vol.l<<-100 # Volume in the liquid phase, L (Berends cite Taweel et al 2006)
Vol.g<<-40  # Volume in the gas phase, L (Berends - no ref)

Ptot <<- 1.01325;                                 # System pressure, bar. (approx same at atm)
R <<- 8.314*1e-2;                                 # bar*L/(mol*K)
T.rumen <<- 39+273.15; 	# Temperature, K

# Henrys constant  M/bar, at T = 25C (298.15K)
KH.co2.s <<-0.035;
KH.ch4.s <<- 0.0014;
KH.h2.s <<- 7.8e-4; 
#convert to rumen temperature
KH.co2 <<- KH.co2.s*exp(-(19410/(R*100))*(1/298.15-1/T.rumen));
KH.ch4 <<- KH.ch4.s*exp(-(14240/(R*100))*(1/298.15-1/T.rumen));
KH.h2  <<- KH.h2.s*exp(-(4180/(R*100))*(1/298.15-1/T.rumen));
KH<<-c('H2.gas'=KH.h2,'CO2.gas'=KH.co2,'CH4.gas'=KH.ch4) #mol/L/bar

# Equilibrium constants 
deltaH0.Ka.w <<- 55900;   deltaH0.Ka.co2 <<- 7646;  deltaH0.Ka.nh4 <<- 51965;

# Acid-base constants mol/L
K.w <<- exp(deltaH0.Ka.w/(R*100)*(1/298.15-1/T.rumen))*1e-14;
K.a.co2 <<- 10^(-6.35)*exp(deltaH0.Ka.co2/(R*100)*(1/298.15-1/T.rumen));
K.a.nh4 <<-  10^(-9.25)*exp(deltaH0.Ka.nh4/(R*100)*(1/298.15-1/T.rumen));
K.a.ac <<- 10^(-4.76)
K.a.bu <<- 10^(-4.82)
K.a.pr <<- 10^(-4.88)
K.a.vfa <<- K.a.ac


