#' getPolymerFrac
#'
#' computes the fractions of NDF, NSC and Protein
#' Dividing feed into polymer fractions (expressed as grams per Kg of DM)
#' As a first estimate NSC is calculated as 1000 - NDF - Protein - Ash - Oil
#' 
#' @param basal.diet either 'Mixed' or 'Concentrate'
#' @param treatment either 'Control','Nitrate' or 'Rapeseed-cake'
#' @return named vector containg fractions: ['NDF','NSC','Protein']
#' @export
#' 
getPolymerFrac=function(basal.diet,treatment,dietCompositionMat){
    

    comp.mat=dietCompositionMat

    print(comp.mat)
    
    composition=comp.mat[comp.mat[,'Additive']==treatment & comp.mat[,'Diet']==basal.diet,]
    
    if (!exists('composition')){stop('composition is not defined for this diet/treatment')}

    NSC=1000-sum(as.numeric(composition[c('NDF','CP','Ash','Oil')]))

    polymer.frac.gPkg=as.numeric(c(composition['NDF'],NSC,composition['CP']))

    names(polymer.frac.gPkg)=c('NDF','NSC','Protein')
    
    return(polymer.frac.gPkg)
}

