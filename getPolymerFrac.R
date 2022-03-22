getPolymerFrac=function(diet,treatment,dietComposition){
    
    #Dividing feed into polymer fractions (expressed as grams per Kg of DM)
    
    #As a first estimate NSC can be calculated as 1000 - NDF - Protein - Ash - Oil (Richard Dewhurst)
    #polymer.frac.gPkg<<-c("NDF"=289,"NSC"=486.8,"Protein"=144) #Troy et al 2015 (mixed, control)

#    comp.mat=rbind(c(144,289,27.2,53),
#        c(133,227,27.7,35))
#    comp.mat=cbind(c('Control','Control'),c('Mixed','Concentrate'),comp.mat)
#    colnames(comp.mat)=c('treatment','basal.diet','CP','NDF','Oil','Ash')

#    comp.mat=read.csv(paste(dataFolder,'DietComposition.csv',sep=''),header=TRUE)

    comp.mat=dietComposition

    composition=comp.mat[comp.mat[,'treatment']==treatment & comp.mat[,'basal.diet']==diet,]
    
    if (!exists('composition')){stop('composition is not defined for this diet/treatment')}

    NSC=1000-sum(as.numeric(composition[c('NDF','CP','Ash','Oil')]))

    polymer.frac.gPkg=as.numeric(c(composition['NDF'],NSC,composition['CP']))

    names(polymer.frac.gPkg)=c('NDF','NSC','Protein')
    
    return(polymer.frac.gPkg)
}

