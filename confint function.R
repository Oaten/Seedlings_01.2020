easyPredCI <- function(model, newdata=NULL, alpha=0.05) {

  # Marc Girondot - 2016-01-09

  if (is.null(newdata)) {

    if (any(class(model)=="glmerMod")) newdata <- model@frame

    if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) newdata <- model$data

    if (any(class(model)=="glmmadmb")) newdata <- model$frame

  }

  

  ## baseline prediction, on the linear predictor scale:

  pred0 <- predict(model, re.form=NA, newdata=newdata)

  ## fixed-effects model matrix for new data

    if (any(class(model)=="glmmadmb")) {

  X <- model.matrix(delete.response(model$terms), newdata)

    } else {

  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],

                    newdata)

    }



  if (any(class(model)=="glm")) {

    # Marc Girondot - 2016-01-09

    # Note that beta is not used

    beta <- model$coefficients

  } else {

    beta <- fixef(model) ## fixed-effects coefficients

  }



  V <- vcov(model)     ## variance-covariance matrix of beta

  

  # Marc Girondot - 2016-01-09

  if (any(!(colnames(V) %in% colnames(X)))) {

    dfi <- matrix(data = rep(0, dim(X)[1]*sum(!(colnames(V) %in% colnames(X)))), nrow = dim(X)[1])

    colnames(dfi) <- colnames(V)[!(colnames(V) %in% colnames(X))]

    X <- cbind(X, dfi)

  }

  

  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions



    ## inverse-link function

  # Marc Girondot - 2016-01-09

  if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) linkinv <- model$family$linkinv

  if (any(class(model)=="glmerMod")) linkinv <- model@resp$family$linkinv

  if (any(class(model)=="glmmadmb")) linkinv <- model$ilinkfun



  ## construct 95% Normal CIs on the link scale and

  ##  transform back to the response (probability) scale:

  crit <- -qnorm(alpha/2)

  linkinv(cbind(lwr=pred0-crit*pred.se,

                upr=pred0+crit*pred.se))

}

