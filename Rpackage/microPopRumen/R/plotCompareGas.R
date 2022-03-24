#' plotCompareGas
#'
#' plots modelled methane production rate against measured data
#'
#' @param out Output from rumenModel()
#' @param spinUpTime.hours time period to run model before using feed data
#' @param gasMat gas measurements data set (time (h), 'MPRmPh' (methane production rate in moles per hour))
#' 

plotCompareGas=function(out,spinUpTime.hours,gasMat=NULL){


    model.time=out$solution[,'time']-min(out$solution[,'time'])
    methane=out$solution[,'CH4.gas']
    MPR=methane*out$parms$Smats$washOut['CH4.gas']

    if (!is.null(gasMat)){
        gasTime=gasMat[,1]+spinUpTime.hours
        mpr=c(MPR,gasMat[,2])
    }else{
        mpr=MPR
    }

    
    grDevices::dev.new()
    graphics::plot(range(model.time),range(mpr),type='n',xlab='Time (h)',ylab='MPR (moles/h)')
    graphics::lines(model.time,MPR,col='black')
    if (!is.null(gasMat)){
        graphics::lines(gasTime,gasMat[,2],col='red')
    }
    graphics::abline(v=spinUpTime.hours,lty=2,col='blue')

    graphics::legend('topleft',c('model','data'),col=c('black','red'),lty=c(1,1),bty='n')

    graphics::text(spinUpTime.hours+0.01*diff(range(model.time)),min(mpr),'end of spin-up',adj=c(0,0),col='blue')
}
