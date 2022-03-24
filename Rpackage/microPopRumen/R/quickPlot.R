#' quickPlot
#'
#' plot up model variables
#'
#' @param out The output from rumenModel()

quickPlot=function(out){
    
    microbeNames=out$parms$microbeNames
    resourceNames=out$parms$resourceNames
    
    Bac=out$solution[,microbeNames]
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    polymerNames=c('NDF','NSC','Protein')
    vfaNames=c('Acetate','Butyrate','Propionate')
    gasNames=c('H2'='H2.gas','SIC'='CO2.gas','CH4'='CH4.gas')

    model.time=out$solution[,'time']-min(out$solution[,'time'])

    grDevices::dev.new()
    graphics::par(mfrow=c(2,2))
    yl=c('microbes','polymers','VFA','Gases (g)')
    pl=list(microbeNames,polymerNames,vfaNames,gasNames)
    for (n in 1:length(pl)){
        graphics::matplot(model.time,out$solution[,pl[[n]]],type='l',
                xlab='time (h)',ylab=yl[n],
                lwd=2,cex.lab=1.2,cex.axis=1.2)
        graphics::legend('topright',legend=pl[[n]],col=1:3,lty=1:3,lwd=2,bty='n',cex=0.8)
    }

}
