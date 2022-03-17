quickPlot=function(out){
    
    microbeNames=out$parms$microbeNames
    resourceNames=out$parms$resourceNames
    time=out$solution[,1]
    Bac=out$solution[,microbeNames]
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    polymerNames=c('NDF','NSC','Protein')
    vfaNames=c('Acetate','Butyrate','Propionate')
    gasNames=c('H2'='H2.gas','SIC'='CO2.gas','CH4'='CH4.gas')
 
    dev.new()
    par(mfrow=c(2,3))
    yl=c('microbes','polymers','VFA','Gases')
    pl=list(microbeNames,polymerNames,vfaNames,gasNames)
    for (n in 1:length(pl)){
        matplot(time/24,out$solution[,pl[[n]]],type='l',
                xlab='time (d)',ylab=yl[n],
                lwd=2,cex.lab=1.2,cex.axis=1.2)
        legend('right',legend=pl[[n]],col=1:3,lty=1:3,lwd=2)
    }

    plot(time/24,MPR,type='l',ylab='MPR (g/d)',xlab='time (d)')
}
