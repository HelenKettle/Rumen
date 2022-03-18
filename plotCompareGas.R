plotCompareGas=function(out,spinUpTime.hours){

    #TSmat with 'Time'(hours),'DMIR' (dry matter intake rate in kg/d),'MPRmPh' (methane production rate in moles per hour))

    #gas data is in myPars[['TSmat']]
    
    model.time=out$solution[,'time']-min(out$solution[,'time'])
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    gasMat=out$myPars[['TSmat']]
    gasTime=gasMat[,'Time']-min(gasMat[,'Time'])
    
    dev.new()
    plot(range(model.time),range(c(MPR,gasMat[,'MPRmPh'])),type='n',xlab='Time (h)',ylab='MPR (moles/h)')
    lines(model.time,MPR,col='black')
    lines(gasTime,gasMat[,'MPRmPh'],col='red')
    abline(v=spinUpTime.hours,lty=2,col='blue')

    legend('topleft',c('model','data','spin up time'),col=c('black','red','blue'),lty=c(1,1,2),bty='n')
    
}
