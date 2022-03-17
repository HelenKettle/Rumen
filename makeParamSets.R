#make parameter sets to run parallel simulations
makeParamSets=function(paramNames,rangesList,numDivs,defaultValuesList=NULL,type='OAT'){


    library(FME)
    
    numParams=length(paramNames)
    
    if (type=='OAT'){
        if (numDivs==1){
            NR=1
        }else{
            NR=numParams*numDivs
        }
        print(paste('number of parameter sets is: ',NR))
    }


    if (type=='OAT'){


        if (numDivs==1){
            mat=matrix(NA,nrow=NR,ncol=numParams)
            colnames(mat)=c(paramNames)
            
            for (par in paramNames){
                mat[1,par]=defaultValuesList[[par]]
            }
            
        }else{

            mat=matrix(NA,nrow=NR,ncol=(numParams+1))
            colnames(mat)=c(paramNames,'parVarName')
            ct=1

            for (par in paramNames){
                
                vec=seq(rangesList[[par]]['min'],rangesList[[par]]['max'],length=numDivs)
                
                                        #take mean/default values for other parameters
                other.pars=paramNames[paramNames!=par]
                defaults=NA*seq(1,(numParams-1))
                names(defaults)=other.pars
                for (op in other.pars){
                    defaults[op]=defaultValuesList[[op]]#mean(rangesList[[op]])
                }
                
                for (val in vec){
                    
                    mat[ct,par]=val
                    mat[ct,other.pars]=defaults
                    mat[ct,'parVarName']=par
                    ct=ct+1
                    
                }
            }
        }

    }else{ #global

        paramNames=rownames(rangesList)
        numPar=length(paramNames)

        if (numPar>4){stop('this function is only written for 4 parameters or less')}

        for (p in paramNames){
            vec=seq(rangesList[p,'min'],rangesList[p,'max'],length=numDivs[p])
            assign(paste(p,'vec',sep='.'),vec)
        }

        mat=NULL
        for (v1 in get(paste(paramNames[1],'vec',sep='.'))){
            for (v2 in get(paste(paramNames[2],'vec',sep='.'))){
                if (numPar>2){
                    for (v3 in get(paste(paramNames[3],'vec',sep='.'))){
                        if (numPar>3){
                            for (v4 in get(paste(paramNames[4],'vec',sep='.'))){
                                mat=rbind(mat,c(v1,v2,v3,v4))
                            }
                        }else{
                            mat=rbind(mat,c(v1,v2,v3))
                        }
                    }
                }else{
                    mat=rbind(mat,c(v1,v2))
                }
            }             
        }

        colnames(mat)=paramNames

    }
    
     
    print(paste('number of parameter sets is: ',nrow(mat)))

    return(mat)
}


    
