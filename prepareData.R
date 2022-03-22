prepareData=function(feedMat,spinUpTime.hours,gasMat){

#returns a matrix: TSmat with 'Time'(hours),'DMIR' (dry matter intake rate in kg/d),'MPRmPh' (methane production rate in moles per hour))

    
#get input data file ('feedRate')
    if (useFeedData){
        JD.feed=feedMat[,1]*24 #change to hours
    }

    if (useGasData){
#get output data file (Gas data eg. MPR or HPR or CO2, in 'animalGasData')
        JD.gas=gasMat[,1]*24 #change to hours
    }
    
#check times (JDs) are the same and if not find overlap and interpolate to 1 min intervals
    one.min=1/60#in hours

    if (useGasData & useFeedData){
        req.times=seq(max(min(JD.gas),min(JD.feed)),min(max(JD.gas),max(JD.feed)),by=one.min)
    }else{
        if (useFeedData){
            req.times=
        
        req.gas.name='CH4.g.min.' 
        req.feed.name='DMIR(kg/d)'

        MPRgPm = approx(JD.gas,animalGasData[,req.gas.name],req.times,ties=mean)$y #g/m 

    #interpolate feed
        feedTS = approx(JD.feed,feedRate[,req.feed.name],req.times,ties=mean)$y/24 #now in kg/h

        molar.mass.CH4=16.04
        TSmat1=cbind('Time'=req.times,'DMIR'=feedTS,'MPRmPh'=60*MPRgPm/molar.mass.CH4) #mol/h
    #add in some spin up time

        req.times.start=req.times[1]
        spin.up.start=req.times.start-spinUpTime.hours
        spin.up.fin=spin.up.start+spinUpTime.hours
        spin.up.vec=seq(spin.up.start,spin.up.fin,by=one.min)
        DMIRbefore=rep(mean(TSmat1[,2]),length(spin.up.vec))
        MPRbefore=rep(mean(TSmat1[,3]),length(spin.up.vec))
        TSmat=rbind(cbind('Time'=spin.up.vec,'DMIR'=DMIRbefore,'MPRmPh'=MPRbefore),TSmat1)
    }
        
    return(TSmat)

}
