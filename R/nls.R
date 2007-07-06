nls <-
  function (formula, data = parent.frame(), start, control = nls.control(),
            algorithm = c("default", "plinear", "port"), trace = FALSE,
            subset, weights, na.action, model = FALSE,
            lower = -Inf, upper = Inf, ...)
{
    ## canonicalize the arguments
  formula <- as.formula(formula)
  algorithm <- match.arg(algorithm)
  
  if(!is.list(data) && !is.environment(data))
    stop("'data' must be a list or an environment")

  mf <- match.call()                  # for creating the model frame
  varNames <- all.vars(formula) # parameter and variable names from formula
  ## for prediction we will need to know those which are in RHS
  form2 <- formula; form2[[2]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)

    ## adjust a one-sided model formula by using 0 as the response
  if (length(formula) == 2) {
    formula[[3]] <- formula[[2]]
    formula[[2]] <- 0
  }

    ## get names of the parameters from the starting values or selfStart model
  pnames <- {
    if (missing(start)) {
      if(!is.null(attr(data, "parameters"))) {
        names(attr(data, "parameters"))
      } else { ## try selfStart - like object
        cll <- formula[[length(formula)]]
        func <- get(as.character(cll[[1]]))
        if(!is.null(pn <- attr(func, "pnames")))
          as.character(as.list(match.call(func, call = cll))[-1][pn])
      }
    } else
    names(start)
  }

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

    ## Heuristics for determining which names in formula represent actual
    ## variables :

    ## If it is a parameter it is not a variable (nothing to guess here :-)
  if(length(pnames))
        varNames <- varNames[is.na(match(varNames, pnames))]
# Can fit a model with pnames even if no varNames 
  if((!length(pnames)) && !length(varNames)) stop("no parameters to fit")
    ## This aux.function needs to be as complicated because
    ## exists(var, data) does not work (with lists or dataframes):
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)),
                                   error = function(e) -1)
  n <- sapply(varNames, lenVar)
  if(any(not.there <- n == -1)) {
    nnn <- names(n[not.there])
    if(missing(start)) {
      if(algorithm == "plinear")
        ## TODO: only specify values for the non-lin. parameters
        stop("No starting values specified")
      ## Provide some starting values instead of erroring out later;
      ## '1' seems slightly better than 0 (which is often invalid):
      warning("No starting values specified for some parameters.\n",
              "Intializing ", paste(sQuote(nnn), collapse=", "),
              " to '1.'.\n",
              "Consider specifying 'start' or using a selfStart model")
      start <- as.list(rep(1., length(nnn)))
      names(start) <- nnn
      varNames <- varNames[i <- is.na(match(varNames, nnn))]
      n <- n[i]
    }
    else # has 'start' but forgot some
      stop("parameters without starting value in 'data': ",
           paste(nnn, collapse=", "))
  }

    ## If its length is a multiple of the response or LHS of the formula,
    ## then it is probably a variable.
    ## This may fail (e.g. when LHS contains parameters):
  respLength <- length(eval(formula[[2]], data, env))
  {
    if(length(n)<1){
#     Some problems might have no offician varNames
      varIndex <- logical(0)
      mf <- list(0)
      wts <- numeric(0)
    }
    else{
      varIndex <- n %% respLength == 0
      if(diff(range(n))>0){
#       'data' is a list that can not be coerced to a data.frame
        mf <- data
        if(missing(start))start <- getInitial(formula, mf)
        startEnv <- new.env(parent = environment(formula))
        for (i in names(start)) assign(i, start[[i]], envir = startEnv)
        rhs <- eval(formula[[3]], data, startEnv)
        n <- length(rhs)
      }
      else{

        mf$formula <-                # replace by one-sided linear model formula
          as.formula(paste("~", paste(varNames[varIndex], collapse = "+")),
                     env = environment(formula))
        mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
        mf$lower <- mf$upper <- NULL
        mf[[1]] <- as.name("model.frame")
        mf <- eval.parent(mf)
        n <- nrow(mf)
        mf <- as.list(mf)
      }
      wts <- {
        if (!mWeights) model.weights(mf)
        else rep(1, n)
      }
      if (any(wts < 0 | is.na(wts))) 
        stop("missing or negative weights not allowed")
    }
  }
# set up iteration 
  if (missing(start)) start <- getInitial(formula, mf)
  for(var in varNames[!varIndex])
    mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[ varNamesRHS %in% varNames[varIndex] ]

  m <- switch(algorithm,
              plinear = stats:::nlsModel.plinear(formula, mf, start, wts),
              port = stats:::nlsModel(formula, mf, start, wts, upper),
              ## Default:
              stats:::nlsModel(formula, mf, start, wts))

  ctrl <- nls.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
# Iterate   
  if (algorithm != "port") {
    if (!missing(lower) || !missing(upper))
      warning('Upper or lower bounds ignored unless algorithm = "port"')
    convInfo <- .Call(stats:::R_nls_iter, m, ctrl, trace)
    nls.out <- list(m = m, convInfo = convInfo,
                    data = substitute(data), call = match.call())
  }
  else { ## "port" i.e., PORT algorithm
    iv <- stats:::nls_port_fit(m, start, lower, upper, control, trace)
    nls.out <- list(m = m, data = substitute(data), call = match.call())
    ## FIXME: this is really a logical for  *NON*convergence:
    nls.out$convergence <- as.integer(if (iv[1] %in% 3:6) 0 else 1)
    nls.out$message <-
      switch(as.character(iv[1]),
             "3" = "X-convergence (3)",
             "4" = "relative convergence (4)",
             "5" = "both X-convergence and relative convergence (5)",
             "6" = "absolute function convergence (6)",
             
             "7" = "singular convergence (7)",
             "8" = "false convergence (8)",
             "9" = "function evaluation limit reached without convergence (9)",
             "10" = "iteration limit reached without convergence (9)",
             "14" = "storage has been allocated (?) (14)",
             
             "15" = "LIV too small (15)",
             "16" = "LV too small (16)",
             "63" = "fn cannot be computed at initial par (63)",
             "65" = "gr cannot be computed at initial par (65)",
             "300" = "initial par violates constraints")
    if (is.null(nls.out$message))
      nls.out$message <-
        paste("See PORT documentation.	Code (", iv[1], ")", sep = "")
    if (nls.out$convergence) {
      msg <- paste("Convergence failure:", nls.out$message)
      if(ctrl$warnOnly) {
        warning(msg)
      } else stop(msg)
    }

    ## we need these (evaluated) for profiling
    nls.out$call$lower <- lower
    nls.out$call$upper <- upper
  }
# Done 
    ## we need these (evaluated) for profiling
  nls.out$call$algorithm <- algorithm
  nls.out$call$control <- ctrl
  nls.out$call$trace <- trace
  
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <-
    attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if(model)
    nls.out$model <- mf
  if(!mWeights)
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- "nls"
  nls.out
}
