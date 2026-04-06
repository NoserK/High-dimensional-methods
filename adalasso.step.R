adalasso.step <-
function(formula, data, family, weights, adaptive.object, s, verbose, selection.criterion) {
    result = list()
    
    #Pull out the relevant data
    response.name = rownames(attr(terms(formula, data=data), 'factors'))[1]
    response.col = which(colnames(data)==response.name)
    predictor.names = attr(terms(formula, data=data), 'term.labels')
    
    #extract the response and the covariates
    y = as.matrix(data[,response.col])
    x = as.matrix(data[,-response.col])
    n <- nrow(x)

    # center and scale the covariates
    result[['meanx']] = adaptive.object[['meanx']]
    result[['scale']] = adaptive.object[['adaweight']]
    x.centered = sweep(x, 2, result[['meanx']], '-')
    x.scaled = sweep(x.centered, 2, result[['scale']], '*')
    
    #family == 'binomial' requires p, 1-p to both be specified:
    if (family == 'binomial') {
        result[['model']] = model = glmnet(x=x.scaled, y=as.matrix(cbind(1-y, y), nrow(x), 2), family=family, weights=weights, lambda=s, standardize=FALSE, intercept=TRUE)
        cv.obj = cv.glmnet(y=as.matrix(cbind(1-y, y), nrow(x), 2), x=x.scaled, nfolds=n, family=family, weights=weights, lambda=s, standardize=FALSE, intercept=TRUE)}
    } else {
        result[['model']] = model = glmnet(x=x.scaled, y=y, family=family, weights=weights, lambda=s, standardize=FALSE, intercept=TRUE)
        cv.obj = cv.glmnet(y=y, x=x.scaled, nfolds=n, family=family, weights=weights, lambda=s, standardize=FALSE, intercept=TRUE)
    }

    #Get the deviance and degrees of freedom, in order to compute the decision criteria
    dev = deviance(model)
    df = sapply(predict(model, type='nonzero'), length) + 1
    
    #Compute the selection criteria
    result[['AIC']] = dev + 2*df
    result[['BIC']] = dev + log(n)*df
    result[['AICc']] = dev + 2*df + 2*df*(df+1)/(n-df-1)
    result[['CV']] = cv.obj
    
    result[['lambda.index']] = lambda.index = which.min(result[[selection.criterion]])
    result[['lambda']] = model[['lambda']][lambda.index]    
    result[['coef']] = coef(result[['model']])

    return(result)
}
