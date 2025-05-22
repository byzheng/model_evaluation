# This script uses the hydroGOF package:
# Zambrano-Bigiarini, M. (2014). hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.
# R package version 0.3-8 (archived). https://cran.r-project.org/web/packages/hydroGOF/index.html



# Some functions to quickly calculate RMSE, PBias, NSE, d, r, and r2 when multiple datasets in one
# These functions were used in the manuscript:
# "Model evaluation: the misuse of statistical techniques when evaluating observations versus predictions"
# Malcolm McPhee, Jonathan Richetti, Barry Croke, and Brad Walmsley
# 2024
# doi: xxx


#' Root Mean Square Error (RMSE)
#'
#' Calculates the Root Mean Square Error between simulated and observed values.
#'
#' @param sim A numeric vector, time series (`ts`), or zoo object with simulated values.
#' @param obs A numeric vector, time series (`ts`), or zoo object with observed values.
#' @param na.rm Logical; should missing values (NA) be removed?
#'
#' @return A single numeric value representing the RMSE.
#'
#' @details
#' The RMSE is a commonly used metric to evaluate the difference between predicted
#' and observed values. It is calculated as:
#' \deqn{\sqrt{\frac{1}{n} \sum_{i=1}^{n} (obs_i - sim_i)^2}}
#'
#' @examples
#' sim <- c(1.1, 2.0, 3.3, 4.2)
#' obs <- c(1.0, 2.1, 3.5, 4.0)
#' rmse(sim, obs)
#'
#' @references
#' Zambrano-Bigiarini, M. (2014). \emph{hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.}
#' R package version 0.3-8 (archived). \url{https://cran.r-project.org/package=hydroGOF}
#'
#' @export
rmse <- function(sim, obs, na.rm = TRUE) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' and 'obs' must be numeric, ts, or zoo.")
    }

    if (length(sim) != length(obs)) {
        stop("Invalid argument: 'sim' and 'obs' must have the same length.")
    }

    sqrt(mean((obs - sim)^2, na.rm = na.rm))
}


#' Percent Bias (PBIAS)
#'
#' Calculates the Percent Bias (PBIAS) between simulated and observed values.
#'
#' @param sim A numeric vector, time series (`ts`), or zoo object with simulated values.
#' @param obs A numeric vector, time series (`ts`), or zoo object with observed values.
#' @param na.rm Logical; should missing values (NA) be removed?
#'
#' @return A numeric value representing the percent bias, rounded to one decimal place.
#'
#' @details
#' PBIAS measures the average tendency of the simulated values to be larger or smaller
#' than their observed counterparts. Ideal value is 0. Positive values indicate model
#' underestimation; negative values indicate overestimation.
#'
#' \deqn{\mathrm{PBIAS} = 100 \times \frac{\sum_{i=1}^{n} (sim_i - obs_i)}{\sum_{i=1}^{n} obs_i}}
#'
#' @examples
#' sim <- c(1.2, 2.4, 3.1, 4.0)
#' obs <- c(1.0, 2.0, 3.0, 4.0)
#' pbias(sim, obs)
#'
#' @references
#' Zambrano-Bigiarini, M. (2014). \emph{hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.}
#' R package version 0.3-8 (archived). \url{https://cran.r-project.org/package=hydroGOF}
#'
#' @export
pbias <- function(sim, obs, na.rm = TRUE) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' and 'obs' must be numeric, ts, or zoo.")
    }

    vi <- valindex(sim, obs)  # user-defined or borrowed from hydroGOF
    obs <- obs[vi]
    sim <- sim[vi]

    denominator <- sum(obs)
    if (denominator == 0) {
        warning("'sum(obs) = 0', it is not possible to compute 'pbias'")
        return(NA_real_)
    }

    pbias <- 100 * sum(sim - obs) / denominator
    round(pbias, 1)
}



valindex = function (sim, obs, ...)
{
    if (length(obs) != length(sim)) {
        stop("Invalid argument: 'length(sim) != length(obs)' !! (",
             length(sim), "!=", length(obs), ") !!")
    }
    else {
        index <- which(!is.na(sim) & !is.na(obs))
        if (length(index) == 0)
            warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
        return(index)
    }
}

#' Nash-Sutcliffe Efficiency (NSE)
#'
#' Computes the Nash-Sutcliffe Efficiency (NSE) to assess the predictive accuracy of a hydrological model.
#'
#' @param sim A numeric vector, time series (`ts`), `zoo`, or `xts` object with simulated values.
#' @param obs A numeric vector, time series (`ts`), `zoo`, or `xts` object with observed values.
#' @param na.rm Logical; should missing values be removed?
#' @param FUN Optional preprocessing function to apply to both `sim` and `obs`.
#' @param epsilon Character or numeric; method to handle zero denominator. Options: `"Pushpalatha2012"`, `"other"`, or `0` (default).
#' @param epsilon.value Numeric value to use if `epsilon = "other"`.
#' @param ... Additional arguments passed to the preprocessing function `FUN`.
#'
#' @return A numeric value representing the Nash-Sutcliffe Efficiency.
#'
#' @details
#' The NSE compares the residual variance to the observed data variance.
#' Values close to 1 indicate good predictive performance.
#'
#' - NSE = 1: perfect model fit
#' - NSE = 0: model is as good as the mean of observations
#' - NSE < 0: model is worse than the mean of observations
#'
#' \deqn{\mathrm{NSE} = 1 - \frac{\sum_{i=1}^{n}(obs_i - sim_i)^2}{\sum_{i=1}^{n}(obs_i - \bar{obs})^2}}
#'
#' @references
#' Nash, J.E. and Sutcliffe, J.V. (1970). River flow forecasting through conceptual models part I — A discussion of principles. *Journal of Hydrology*, 10(3), 282–290.
#'
#' Zambrano-Bigiarini, M. (2014). \emph{hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.}
#' R package version 0.3-8 (archived). \url{https://cran.r-project.org/package=hydroGOF}
#'
#' @export
nse <- function(sim, obs, na.rm = TRUE, FUN = NULL,
                epsilon = c(0, "Pushpalatha2012", "other"),
                epsilon.value = NA, ...) {

    # Check input types
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo", "xts")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo", "xts"))) {
        stop("Invalid input: 'sim' and 'obs' must be numeric, ts, zoo, or xts.")
    }

    # Remove NAs if requested
    vi <- valindex(sim, obs)  # Ensure valindex is defined elsewhere
    sim <- sim[vi]
    obs <- obs[vi]

    # Preprocessing (optional)
    if (!is.null(FUN)) {
        new <- preproc(sim = sim, obs = obs, FUN = FUN,
                       epsilon = epsilon, epsilon.value = epsilon.value, ...)
        sim <- new[["sim"]]
        obs <- new[["obs"]]
    }

    # NSE calculation
    denominator <- sum((obs - mean(obs))^2)
    if (denominator == 0) {
        warning("Denominator is zero: cannot compute NSE.")
        return(NA_real_)
    }

    NS <- 1 - sum((obs - sim)^2) / denominator
    return(NS)
}


#' Willmott's Index of Agreement (d)
#'
#' Computes the Index of Agreement, a standardized measure of the degree of model prediction error.
#' The index varies between 0 and 1, where 1 indicates perfect agreement and 0 indicates no agreement.
#'
#' @param sim Numeric vector, `ts`, or `zoo` object of simulated values.
#' @param obs Numeric vector, `ts`, or `zoo` object of observed values.
#' @param na.rm Logical; should missing values be removed? (currently unused, but included for consistency)
#' @param ... Additional arguments (currently ignored).
#'
#' @return Numeric value between 0 and 1 representing the Index of Agreement.
#'
#' @details
#' The Index of Agreement (d) was developed by Willmott (1981) as a standardized measure of model prediction error.
#' \deqn{
#' d = 1 - \frac{\sum (obs_i - sim_i)^2}{\sum (|sim_i - \bar{obs}| + |obs_i - \bar{obs}|)^2}
#' }
#'
#' @references
#' Willmott, C.J. (1981). On the validation of models. *Physical Geography*, 2(2), 184–194.
#'
#' Zambrano-Bigiarini, M. (2014). \emph{hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.}
#' R package version 0.3-8 (archived). \url{https://cran.r-project.org/package=hydroGOF}
#'
#' @export
d <- function(sim, obs, na.rm = TRUE, ...) {
    # Input class check
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }

    # Remove missing values and align
    vi <- valindex(sim, obs)  # Assumes valindex() returns valid indices
    sim <- sim[vi]
    obs <- obs[vi]

    # Convert ts or zoo to numeric
    if (inherits(sim, c("ts", "zoo"))) sim <- as.numeric(sim)
    if (inherits(obs, c("ts", "zoo"))) obs <- as.numeric(obs)

    Om <- mean(obs)
    denominator <- sum((abs(sim - Om) + abs(obs - Om))^2)

    if (denominator == 0) {
        warning("'sum((abs(sim - Om) + abs(obs - Om))^2) = 0', cannot compute Index of Agreement")
        return(NA_real_)
    }

    d_val <- 1 - (sum((obs - sim)^2) / denominator)
    return(d_val)
}

#' Pearson's Correlation Coefficient (r)
#'
#' Computes the Pearson correlation coefficient between simulated and observed values.
#'
#' @param sim Numeric vector, `ts`, or `zoo` object of simulated values.
#' @param obs Numeric vector, `ts`, or `zoo` object of observed values.
#' @param ... Additional arguments passed to \code{\link[stats]{cor}} (currently ignored).
#'
#' @return Numeric value of Pearson's correlation coefficient.
#'
#' @seealso \code{\link[stats]{cor}}
#'
#' @export
rPearson <- function(sim, obs, ...) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }
    cor(sim, obs, method = "pearson", use = "pairwise.complete.obs")
}

#' Modelling Efficiency (MEF)
#'
#' Computes the Modelling Efficiency to evaluate the goodness-of-fit of simulated values against observed values.
#'
#' @param sim Numeric vector, `ts`, or `zoo` object of simulated values.
#' @param obs Numeric vector, `ts`, or `zoo` object of observed values.
#' @param na.rm Logical; should missing values be removed? (currently unused, but included for consistency)
#' @param ... Additional arguments (currently ignored).
#'
#' @return Numeric value representing the Modelling Efficiency.
#'
#' @details
#' MEF is defined as:
#' \deqn{
#' \mathrm{MEF} = 1 - \frac{\sum (obs_i - sim_i)^2}{\sum (obs_i - \bar{obs})^2}
#' }
#' where \eqn{obs_i} are observed values, \eqn{sim_i} are simulated values, and \eqn{\bar{obs}} is the mean of observed values.
#'
#' @references
#' Zambrano-Bigiarini, M. (2014). \emph{hydroGOF: Goodness-of-Fit Functions for Comparison of Simulated and Observed Hydrological Time Series.}
#' R package version 0.3-8 (archived). \url{https://cran.r-project.org/package=hydroGOF}
#'
#' @export
mef <- function(sim, obs, na.rm = TRUE, ...) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }
    if (length(sim) != length(obs)) {
        stop("'sim' and 'obs' must have the same length.")
    }

    vi <- valindex(sim, obs)  # Assumes valindex returns valid indices, excluding NAs
    sim <- sim[vi]
    obs <- obs[vi]

    # Convert ts or zoo to numeric
    if (inherits(sim, c("ts", "zoo"))) sim <- as.numeric(sim)
    if (inherits(obs, c("ts", "zoo"))) obs <- as.numeric(obs)

    1 - sum((obs - sim)^2) / sum((obs - mean(obs))^2)
}

#' Mean Squared Error of Prediction (MSEP)
#'
#' Computes the mean squared error between simulated and observed values.
#'
#' @param sim Numeric vector, ts, or zoo of simulated values.
#' @param obs Numeric vector, ts, or zoo of observed values.
#' @param na.rm Logical; currently unused, but reserved for NA handling.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Numeric MSEP value.
#'
#' @export
msep <- function(sim, obs, na.rm = TRUE, ...) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }
    if (length(sim) != length(obs)) {
        stop("'sim' and 'obs' must have the same length.")
    }

    vi <- valindex(sim, obs)
    sim <- sim[vi]
    obs <- obs[vi]

    if (inherits(sim, c("ts", "zoo"))) sim <- as.numeric(sim)
    if (inherits(obs, c("ts", "zoo"))) obs <- as.numeric(obs)

    mean((obs - sim)^2)
}

#' MSEP Bias Component
#'
#' Computes the bias component of the Mean Squared Error of Prediction.
#'
#' @inheritParams msep
#'
#' @return Numeric value of the bias component of MSEP.
#'
#' @export
msep_bias <- function(sim, obs, na.rm = TRUE, ...) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }
    if (length(sim) != length(obs)) {
        stop("'sim' and 'obs' must have the same length.")
    }

    vi <- valindex(sim, obs)
    sim <- sim[vi]
    obs <- obs[vi]

    if (inherits(sim, c("ts", "zoo"))) sim <- as.numeric(sim)
    if (inherits(obs, c("ts", "zoo"))) obs <- as.numeric(obs)

    (mean(sim) - mean(obs))^2
}

#' MSEP Slope Component
#'
#' Computes the slope component of the Mean Squared Error of Prediction.
#'
#' @inheritParams msep
#'
#' @return Numeric value of the slope component of MSEP.
#'
#' @export
msep_slope <- function(sim, obs, na.rm = TRUE, ...) {
    if (!inherits(sim, c("numeric", "integer", "ts", "zoo")) ||
        !inherits(obs, c("numeric", "integer", "ts", "zoo"))) {
        stop("Invalid argument type: 'sim' & 'obs' must be numeric, ts, or zoo objects.")
    }
    if (length(sim) != length(obs)) {
        stop("'sim' and 'obs' must have the same length.")
    }

    vi <- valindex(sim, obs)
    sim <- sim[vi]
    obs <- obs[vi]

    if (inherits(sim, c("ts", "zoo"))) sim <- as.numeric(sim)
    if (inherits(obs, c("ts", "zoo"))) obs <- as.numeric(obs)

    fit <- lm(obs ~ sim)
    slope <- coef(fit)[2]

    ((sum((sim - mean(sim))^2) / length(obs)) * (1 - slope)^2)
}


#' Calculate the Random (Deviance) Component of MSEP
#'
#' Computes the portion of the Mean Squared Error of Prediction (MSEP)
#' attributed to random error (unexplained variance), based on the squared correlation
#' between observed and simulated values.
#'
#' This function assumes `sim` and `obs` are numeric vectors or time series objects
#' of equal length. It first filters paired valid (non-NA) indices using `valindex()`.
#'
#' @param sim Numeric vector or time series of simulated values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently not used).
#'
#' @return Numeric scalar representing the random (deviance) component of MSEP.
#' Returns `NA` if variance is zero or correlation is undefined.
#'
#' @details
#' The random component of MSEP is calculated as:
#' \deqn{
#' \text{MSEP}_{random} = (1 - R^2) \times \frac{1}{n} \sum_{i=1}^n (obs_i - \bar{obs})^2
#' }
#' where \eqn{R^2} is the squared Pearson correlation coefficient between `obs` and `sim`.
#'
#' @examples
#' obs <- c(1, 2, 3, 4, 5)
#' sim <- c(1.1, 1.9, 3.1, 3.9, 5.2)
#' msep_random(sim, obs)
#'
#' @export
msep_random <- function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    if(length(obs) == 0) stop("No valid observations after filtering")

    if(var(obs) == 0 || var(sim) == 0) {
        warning("Variance of obs or sim is zero, correlation undefined. Returning NA.")
        return(NA_real_)
    }

    cor_val <- cor(obs, sim)
    if(is.na(cor_val)) {
        warning("Correlation returned NA. Returning NA.")
        return(NA_real_)
    }

    msep_random <- (1 - cor_val^2) * (sum((obs - mean(obs))^2) / length(obs))
    return(msep_random)
}


#' Sum of MSEP Components: Bias, Slope, and Random
#'
#' Computes the sum of the three components of the Mean Squared Error of Prediction (MSEP):
#' the bias component, the slope component, and the random (deviance) component.
#'
#' @param sim Numeric vector or time series of simulated values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric scalar representing the sum of the bias, slope, and random MSEP components.
#'
#' @details
#' This function uses the helper functions `msep_bias()`, `msep_slope()`, and `msep_random()`
#' to compute the individual components of MSEP and returns their sum.
#' It performs input validation and filters paired valid observations before calculation.
#'
#' @examples
#' obs <- c(1, 2, 3, 4, 5)
#' sim <- c(1.1, 1.9, 3.1, 3.9, 5.2)
#' SumMSEP123(sim, obs)
#'
#' @export
SumMSEP123 <- function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    SumMSEP123 <- msep_bias(obs, sim) + msep_slope(obs, sim) + msep_random(obs, sim)
    return(SumMSEP123)
}

#' Proportion of MSEP Due to Bias
#'
#' Calculates the percentage proportion of the Mean Squared Error of Prediction (MSEP)
#' that is attributable to bias in the predictions.
#'
#' @param sim Numeric vector or time series of simulated (predicted) values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric scalar representing the percentage proportion of MSEP due to bias.
#'
#' @details
#' The function first validates the inputs and synchronizes the `sim` and `obs` vectors
#' to ensure they have the same length and no missing paired values.
#' It uses the helper functions `msep_bias()` to calculate the bias component of MSEP, and `msep()`
#' to calculate the total MSEP, then returns the ratio as a percentage.
#'
#' @examples
#' obs <- c(10, 20, 30, 40, 50)
#' sim <- c(11, 19, 29, 42, 48)
#' prop_bias(sim, obs)
#'
#' @export
prop_bias = function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    prop_bias <- (msep_bias(sim, obs) / msep(sim, obs)) * 100
    return(prop_bias)
}



#' Proportion of MSEP Due to Slope
#'
#' Calculates the percentage proportion of the Mean Squared Error of Prediction (MSEP)
#' that is attributable to the slope component of the prediction error.
#'
#' @param sim Numeric vector or time series of simulated (predicted) values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric scalar representing the percentage proportion of MSEP due to slope.
#'
#' @details
#' The function validates the input vectors and synchronizes them by removing pairs with missing values.
#' It computes the slope component of MSEP using `msep_slope()` and the total MSEP using `msep()`,
#' then returns their ratio as a percentage.
#'
#' @examples
#' obs <- c(10, 20, 30, 40, 50)
#' sim <- c(11, 19, 29, 42, 48)
#' prop_slope(sim, obs)
#'
#' @export
prop_slope = function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    prop_slope <- (msep_slope(sim, obs) / msep(sim, obs)) * 100
    return(prop_slope)
}

#' Proportion of MSEP Due to Deviance
#'
#' Calculates the percentage proportion of the Mean Squared Error of Prediction (MSEP)
#' that is attributable to the deviance (random error) component of the prediction error.
#'
#' @param sim Numeric vector or time series of simulated (predicted) values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric scalar representing the percentage proportion of MSEP due to deviance (random error).
#'
#' @details
#' This function first validates and synchronizes input vectors by removing missing value pairs.
#' It calculates the deviance component of MSEP using `msep_random()` and divides it by the total
#' MSEP computed by `msep()`, then converts the ratio to a percentage.
#'
#' @examples
#' obs <- c(10, 20, 30, 40, 50)
#' sim <- c(11, 19, 29, 42, 48)
#' prop_deviance(sim, obs)
#'
#' @export
prop_deviance = function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    prop_deviance <- (msep_random(sim, obs) / msep(sim, obs)) * 100
    return(prop_deviance)
}

#' Paired T-Test for Mean Bias
#'
#' Performs a paired t-test between simulated and observed values to assess
#' if there is a significant mean bias. Returns the p-value of the test.
#'
#' @param sim Numeric vector or time series of simulated (predicted) values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments passed to \code{t.test} (currently unused).
#'
#' @return Numeric scalar representing the p-value of the paired t-test.
#'
#' @details
#' The function first validates the input classes and lengths, then synchronizes
#' the data by removing NA pairs using \code{valindex}. It performs a paired
#' t-test between \code{sim} and \code{obs} to test the null hypothesis that
#' their means are equal. A p-value less than 0.05 typically indicates
#' a significant mean bias.
#'
#' @examples
#' obs <- c(10, 20, 30, 40, 50)
#' sim <- c(11, 19, 29, 42, 48)
#' P_mean_ttest(sim, obs)
#'
#' @export
P_mean_ttest <- function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    P_slope_ttest <- t.test(sim, obs, paired = TRUE)
    return(P_slope_ttest$p.value)
}

#' Beta Coefficient from Linear Regression of Observed on Simulated
#'
#' Computes the slope (beta coefficient) from a linear regression of observed values on simulated values.
#' This coefficient represents the relationship strength and direction between the two variables.
#'
#' @param sim Numeric vector or time series of simulated (predicted) values.
#' @param obs Numeric vector or time series of observed values.
#' @param na.rm Logical; if TRUE, missing values are removed before computation. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric scalar representing the beta coefficient (slope) from the regression \code{obs ~ sim}.
#'
#' @details
#' The function validates input types and length, synchronizes data by removing missing values,
#' and fits a linear model with observed values as the response and simulated values as the predictor.
#' The returned beta coefficient corresponds to the slope of this linear relationship.
#'
#' @examples
#' obs <- c(10, 20, 30, 40, 50)
#' sim <- c(11, 19, 29, 42, 48)
#' beta(sim, obs)
#'
#' @export
beta = function (sim, obs, na.rm = TRUE, ...)
{
    if (is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
        is.na(match(class(obs), c("integer", "numeric", "ts", "zoo"))))
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim))
        stop("Invalid argument: 'sim' & 'obs' don't have the same length!")

    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.na(match(class(sim), c("ts", "zoo"))))
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo"))))
        obs <- as.numeric(obs)

    fitted_lm <- summary(lm(obs ~ sim))
    beta <- fitted_lm$coefficients[2, 1]
    return(beta)
}


# wrapper functions ----
#' Calculate Multiple Goodness-of-Fit Metrics Grouped by a Factor
#'
#' This function takes a data frame and computes a suite of goodness-of-fit metrics
#' between observed and predicted values, grouped by a specified factor.
#'
#' @param df A data frame containing observed and predicted values and grouping variable.
#' @param obs Unquoted column name of the observed values in `df`.
#' @param pred Unquoted column name of the predicted (simulated) values in `df`.
#' @param group Unquoted column name of the grouping variable in `df`.
#'
#' @return A transposed data frame where rows correspond to different metrics and columns correspond to groups.
#' The first row with group names is removed for a cleaner output.
#'
#' @details
#' The function calculates multiple metrics such as MSEP and its decomposition,
#' RMSE, Mean Bias, Nash-Sutcliffe Efficiency (NSE), Pearson's correlation, t-test p-value,
#' beta coefficient from regression, and others, for each group.
#'
#' Requires the `tidyverse` package and all the metric functions
#' (e.g., `msep`, `rmse`, `msep_bias`, etc.) to be loaded in the environment.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(group = rep(c("A", "B"), each = 5),
#'                  obs = c(1,2,3,4,5, 2,3,4,5,6),
#'                  pred = c(1.1,1.9,3.2,3.8,4.9, 2.1,2.9,3.8,5.1,6.2))
#' all_metrics_calc(df, obs, pred, group)
#' }
#'
#' @importFrom tidyverse enquo
#' @importFrom dplyr group_by reframe
#' @export
all_metrics_calc = function(df, obs, pred, group){
    library(tidyverse)

    group = enquo(group)
    obs = enquo(obs)
    pred = enquo(pred)

    metrics = df %>%
        group_by(!!group) %>%
        reframe(
            n = length(!!pred),
            Obs_NA = sum(is.na(!!obs)),
            Pred_NA = sum(is.na(!!pred)),
            obs_mean = mean(!!obs, na.rm=TRUE),
            pred_mean = mean(!!pred, na.rm=TRUE),
            MeanBias = mean(!!obs - !!pred, na.rm=TRUE),
            MSEP = msep(!!obs, !!pred),
            RMSE = rmse(!!obs, !!pred),
            MSEP_bias = msep_bias(!!obs, !!pred),
            MSEP_slope = msep_slope(!!obs, !!pred),
            MSEP_random = msep_random(!!obs, !!pred),
            SumMSEP = SumMSEP123(!!obs, !!pred),
            Prop_BIAS = prop_bias(!!obs, !!pred),
            Prop_SLOPE = prop_slope(!!obs, !!pred),
            Prop_DEVIANCE = prop_deviance(!!obs, !!pred),
            MEF = mef(!!obs, !!pred),
            PBIAS = pbias(!!obs, !!pred),
            d = d(!!obs, !!pred),
            NSE = nse(!!obs, !!pred),
            r = rPearson(!!obs, !!pred),
            r2 = r^2,
            P_mean = P_mean_ttest(!!obs, !!pred),
            beta_coef = beta(!!obs, !!pred)
        )

    transposed_metrics = t(metrics)
    transposed_metrics = transposed_metrics[-1,]  # Remove group names row

    return(transposed_metrics)
}

#' Calculate Recommended Goodness-of-Fit Metrics Grouped by a Factor
#'
#' This function groups the dataframe by a specified grouping variable and calculates a subset of
#' recommended goodness-of-fit metrics between observed and predicted values.
#'
#' @param df Data frame containing observed, predicted, and grouping columns.
#' @param obs Unquoted name of the observed values column.
#' @param pred Unquoted name of the predicted values column.
#' @param group Unquoted name of the grouping variable.
#'
#' @return A transposed data frame of metrics (rows) by groups (columns).
#'
#' @export
recommended_metrics_calc = function(df, obs, pred, group){
    library(tidyverse)

    group <- enquo(group)
    obs <- enquo(obs)
    pred <- enquo(pred)

    metrics <- df %>%
        group_by(!!group) %>%
        reframe(
            n = length(!!pred),
            obs_mean = mean(!!obs, na.rm = TRUE),
            pred_mean = mean(!!pred, na.rm = TRUE),
            MeanBias = mean(!!obs - !!pred, na.rm = TRUE),
            RMSE = rmse(!!obs, !!pred),
            Prop_BIAS = prop_bias(!!obs, !!pred),
            Prop_SLOPE = prop_slope(!!obs, !!pred),
            Prop_DEVIANCE = prop_deviance(!!obs, !!pred),
            NSE_MEF = mef(!!obs, !!pred)
        )

    transposed_metrics <- t(metrics)
    transposed_metrics <- transposed_metrics[-1,]  # remove group names row

    return(transposed_metrics)
}

metrics_calc = function(df, obs, pred, group){
    # Calculates selected goodness-of-fit metrics grouped by 'group'
    # Args:
    #   df: Data frame with data
    #   obs: unquoted name of observed values column
    #   pred: unquoted name of predicted values column
    #   group: unquoted name of grouping column
    #
    # Returns:
    #   Transposed data frame of metrics by group

    library(tidyverse)

    group <- enquo(group)
    obs <- enquo(obs)
    pred <- enquo(pred)

    metrics <- df %>%
        group_by(!!group) %>%
        reframe(
            n = length(!!pred),
            RMSE = rmse(!!obs, !!pred),         # assuming rmse handles NAs or add na.rm if supported
            PBIAS = pbias(!!obs, !!pred),
            d = d(!!obs, !!pred),
            NSE = nse(!!obs, !!pred),
            r2 = rPearson(!!obs, !!pred)^2
        )

    transposed_metrics <- t(metrics)
    transposed_metrics <- transposed_metrics[-1,]  # drop the group labels row

    return(transposed_metrics)
}

#########################################################
## little example on how the function expects the data ##

#datasets = rep(1:4, each=5)
#obs = runif(n=20, min=1, max=20)
#pred = runif(n=20, min=1, max=20)
#df = data.frame(datasets,obs,pred)
#names(df) = c('Datasets','Observations','Simulations')

#all_metrics_calc(df, Observations, Simulations, Datasets)
