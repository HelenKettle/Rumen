plotCompareGas=function(out,spinUpTime.hours,useGasData=TRUE){

    #TSmat with 'Time'(hours),'DMIR' (dry matter intake rate in kg/d),'MPRmPh' (methane production rate in moles per hour))

    #gas data is in myPars[['TSmat']]
    
    model.time=out$solution[,'time']-min(out$solution[,'time'])
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    if (useGasData){
        gasMat=out$myPars[['TSmat']]
        gasTime=gasMat[,'Time']-min(gasMat[,'Time'])
        mpr=c(MPR,gasMat[,'MPRmPh'])
    }else{
        mpr=MPR
    }

    
    dev.new()
    plot(range(model.time),range(mpr),type='n',xlab='Time (h)',ylab='MPR (moles/h)')
    lines(model.time,MPR,col='black')
    if (useGasData){
        lines(gasTime,gasMat[,'MPRmPh'],col='red')
    }
    abline(v=spinUpTime.hours,lty=2,col='blue')

    legend('topleft',c('model','data','spin up time'),col=c('black','red','blue'),lty=c(1,1,2),bty='n')
    
}
