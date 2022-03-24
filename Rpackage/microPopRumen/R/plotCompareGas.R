plotCompareGas=function(out,spinUpTime.hours,gasMat=NULL){

    #TSmat with 'Time'(hours),'DMIR' (dry matter intake rate in kg/d),'MPRmPh' (methane production rate in moles per hour))

    #gas data is in myPars[['TSmat']]
    
    model.time=out$solution[,'time']-min(out$solution[,'time'])
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    if (!is.null(gasMat)){
        gasTime=gasMat[,1]+spinUpTime.hours
        mpr=c(MPR,gasMat[,2])
    }else{
        mpr=MPR
    }

    
    dev.new()
    plot(range(model.time),range(mpr),type='n',xlab='Time (h)',ylab='MPR (moles/h)')
    lines(model.time,MPR,col='black')
    if (!is.null(gasMat)){
        lines(gasTime,gasMat[,2],col='red')
    }
    abline(v=spinUpTime.hours,lty=2,col='blue')

    legend('topleft',c('model','data'),col=c('black','red'),lty=c(1,1),bty='n')

    text(spinUpTime.hours+0.01*diff(range(model.time)),min(mpr),'end of spin-up',adj=c(0,0),col='blue')
}
