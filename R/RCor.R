#' Kendall, Spearman and Pearson correlation and their generalizations for non-continuous data
#'
#' `RCor()` computes the specified correlation with corresponding confidence intervals and P-values for the associated independence test either in the iid or in the time series case.
#'
#' @param X a n x 1 numeric vector, matrix or data frame.
#' @param Y a n x 1 numeric vector, matrix or data frame.
#' @param alpha a numeric value specifying the significance level. The confidence level will be 1 - alpha.
#' @param method a character string specifying the correlation coefficient to be used for the independence test. Possible values are "tau", "tau_b", "tau_b_mod", "gamma", "rho", "rho_b" and "r". The recommendation for data with ties is "gamma". Specifying "tau_b_mod" only yields the independence test for IID data.
#' @param IID logical indicator determining whether the inference shall be conducted under iid (default) or time series assumptions (see CITATION for a precise description of the assumptions)
#' @param Fisher logical indicator determining whether the CIs shall be computed by using the Fisher transformation.

#' @return The value of the chosen correlation coefficient along with its confidence interval and an independence test as well as an uncorrelatedness test.
#' @export
#'
#' @import stats
#' @import DescTools
#' @import Matrix
#' @import dplyr
#' @examples
#' X <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
#' Y <- c(1, 1, 1, 1, 2, 1, 2, 2, 2)
#' RCor(X, Y)
RCor <- function(X, Y, alpha = 0.1, method = "gamma", IID = TRUE, Fisher = TRUE){
  if (!(is.numeric(X) && is.numeric(Y) && length(X) == length(Y))){
    stop("`X` and `Y` must be numeric vectors of the same length", call. = FALSE)
  }
  X <- X
  Y <- Y
  n <- length(X)
  if (method == "tau"){
    tau_info <- DescTools:::.DoCount(X, Y)
    tau <- (tau_info$C - tau_info$D) / choose(n, 2)
    tau_fis <- atanh(tau)
    X_TieProb3 <- sum((table(X)/length(X))^3)
    Y_TieProb3 <- sum((table(Y)/length(Y))^3)
    if (isTRUE(IID)){
      # Calculate Marc's variance estimator
      var_hat <- 4 * mean((4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - tau)^2)
      # Variance under independence assumption
      var_hat_ind <- 4/9 * (1 - X_TieProb3) * (1 - Y_TieProb3)
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * tau / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * tau / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- Tau_LRV(X, Y, tau)
      var_hat_ind <- Tau_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * tau / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * tau / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(tau + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), tau + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(tau_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - tau^2)), tanh(tau_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - tau^2)))
    }
    res <- dplyr::tribble(~Tau, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          tau, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
  if (method == "tau_b"){
    tau_info <- DescTools:::.DoCount(X, Y)
    tau <- (tau_info$C - tau_info$D) / choose(n, 2)
    tau_b <- stats::cor(X, Y, method = "kendall")
    tau_b_fis <- atanh(tau_b)
    X_TieProb <- sum((table(X)/length(X))^2)
    Y_TieProb <- sum((table(Y)/length(Y))^2)
    X_TieProb3 <- sum((table(X)/length(X))^3)
    Y_TieProb3 <- sum((table(Y)/length(Y))^3)
    if (isTRUE(IID)){
      # Calculate Marc's variance estimator
      G_XYXY <- G_XY(X, Y)
      G_XX <- G_X(X)
      G_YY <- G_Y(Y)
      x_eqX <- x_eq(X)
      y_eqY <- y_eq(Y)
      var_tau <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau)^2)
      var_xix <- 4 * mean((x_eqX - X_TieProb)^2)
      var_xiy <- 4 * mean((y_eqY - Y_TieProb)^2)
      var_tauxix <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX - X_TieProb))
      var_tauxiy <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (y_eqY - Y_TieProb))
      var_xixxiy <- 4 * mean((x_eqX - X_TieProb) * (y_eqY - Y_TieProb))
      var_hat <- (var_tau + tau * (var_tauxix / (1 - X_TieProb) + var_tauxiy / (1 - Y_TieProb)) + tau^2 / 4 * (var_xix / (1 - X_TieProb)^2 + var_xiy / (1 - Y_TieProb)^2 + (2 * var_xixxiy) / (1 - Y_TieProb) / (1 - X_TieProb))) / ((1 - X_TieProb) * (1 - Y_TieProb))
      # Variance under independence assumption
      var_hat_ind <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3) / (1 - X_TieProb) / (1 - Y_TieProb)
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * tau_b / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * tau_b / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- TauB_LRV(X, Y, tau, X_TieProb, Y_TieProb, bandwidth = "Dehling")
      var_hat_ind <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3) / (1 - X_TieProb) / (1 - Y_TieProb) * Rhob_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * tau_b / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * tau_b / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(tau_b + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), tau_b + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(tau_b_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - tau_b^2)), tanh(tau_b_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - tau_b^2)))
    }
    res <- dplyr::tribble(~TauB, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          tau_b, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
  if (method == "tau_b_mod"){
    ties_x <- 0
    ties_y <- 0
    for (i in 3:n) {
      for (j in 2:(i - 1)) {
        for (k in 1:(j - 1)) {
          ties_x <- ties_x + ifelse(X[i] == X[j] & X[j] == X[k], 1, 0)
          ties_y <- ties_y + ifelse(Y[i] == Y[j] & Y[j] == Y[k], 1, 0)
        }
      }
    }
    ties_x <- ties_x / choose(n, 3)
    ties_y <- ties_y / choose(n, 3)
    tau_info <- DescTools:::.DoCount(X, Y)
    tau <- (tau_info$C - tau_info$D) / choose(n, 2)
    tau_b_mod <- tau / sqrt((1 - ties_x) * (1 - ties_y))
    if (isTRUE(IID)){
      var_hat_ind <- 4/9
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * tau_b_mod / sqrt(var_hat_ind)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    res <- dplyr::tribble(~TauB_Mod, ~PValueIND,
                          tau_b_mod, p_val_ind)
    return(res)
  }
  if (method == "gamma"){
    gamma_info <- DescTools:::.DoCount(X, Y)
    gamma <- (gamma_info$C - gamma_info$D) / (gamma_info$C + gamma_info$D)
    gamma_fis <- atanh(gamma)
    tau <- (gamma_info$C - gamma_info$D) / choose(n, 2)
    X_TieProb <- sum((table(X)/length(X))^2)
    Y_TieProb <- sum((table(Y)/length(Y))^2)
    X_TieProb3 <- sum((table(X)/length(X))^3)
    Y_TieProb3 <- sum((table(Y)/length(Y))^3)
    XY_TieProb <- sum((table(X, Y)/length(X))^2)
    tie_prob <- X_TieProb + Y_TieProb - XY_TieProb
    if (isTRUE(IID)){
      # Calculate Marc's variance estimator
      G_XYXY <- G_XY(X, Y)
      G_XX <- G_X(X)
      G_YY <- G_Y(Y)
      x_eqX <- x_eq(X)
      y_eqY <- y_eq(Y)
      x_eq_y_eqXY <- x_eq_y_eq(X, Y)
      var_tau <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau)^2)
      var_nu <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob)^2)
      var_taunu <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX + y_eqY - x_eq_y_eqXY - tie_prob))
      var_hat <- (var_tau + gamma^2 * var_nu + 2 * gamma * var_taunu) / (1 - tie_prob)^2
      # Variance under independence assumption
      var_hat_ind <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3) / (1 - X_TieProb)^2 / (1 - Y_TieProb)^2
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * gamma / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * gamma / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- Gamma_LRV(X, Y, tau, tie_prob, bandwidth = "Dehling")
      var_hat_ind <- 4 / 9 * (1 - X_TieProb3) * (1 - Y_TieProb3) / (1 - X_TieProb)^2 / (1 - Y_TieProb)^2 * Rhob_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * gamma / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * gamma / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(gamma + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), gamma + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(gamma_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - gamma^2)), tanh(gamma_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - gamma^2)))
    }
    res <- dplyr::tribble(~Gamma, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          gamma, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
  if (method == "rho"){
    rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
    rho_fis <- atanh(rho)
    X_TieProb3 <- sum((table(X)/length(X))^3)
    Y_TieProb3 <- sum((table(Y)/length(Y))^3)
    if (isTRUE(IID)){
      # Calculate Marc's variance estimator
      G_XX <- G_X(X)
      G_YY <- G_Y(Y)
      var_hat <- 9 * mean((4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY)  + 1 - rho)^2)
      # Variance under independence assumption
      var_hat_ind <- (1 - X_TieProb3) * (1 - Y_TieProb3)
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * rho / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * rho / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- SRho_LRV(X, Y, rho, bandwidth = "Dehling")
      var_hat_ind <- 9 / 4 * Tau_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * rho / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * rho / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(rho + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), rho + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(rho_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - rho^2)), tanh(rho_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - rho^2)))
    }
    res <- dplyr::tribble(~Rho, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          rho, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
  if (method == "rho_b"){
    rho_b <- stats::cor(X, Y, method = "spearman")
    rho_b_fis <- atanh(rho_b)
    rho <- 12 * (n - 1) / n^3 * stats::cov(X, Y, method = "spearman")
    rho_x <- 12 * (n - 1) / n^3 * stats::cov(X, X, method = "spearman")
    rho_y <- 12 * (n - 1) / n^3 * stats::cov(Y, Y, method = "spearman")
    X_TieProb <- sum((table(X)/length(X))^2)
    Y_TieProb <- sum((table(Y)/length(Y))^2)
    if (isTRUE(IID)){
      # Calculate Marc's variance estimator
      g_xX <- g_x(X)
      g_yY <- g_y(Y)
      G_XX <- G_X(X)
      G_YY <- G_Y(Y)
      F_XX <- F_X(X)
      F_X_X <- F_X_(X)
      F_YY <- F_Y(Y)
      F_Y_Y <- F_Y_(Y)
      # Calculate min functions
      min1XX <- min1X(X)
      min2XX <- min2X(X)
      min3XX <- min3X(X)
      min4XX <- min4X(X)
      min1YY <- min1Y(Y)
      min2YY <- min2Y(Y)
      min3YY <- min3Y(Y)
      min4YY <- min4Y(Y)
      var_rho <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY)  + 1 - rho)^2)
      var_rhox <- 9 * mean((2 * (min1XX + min2XX + min3XX + min4XX + 2 * G_XX^2 - 4 * G_XX) + 1 - rho_x)^2)
      var_rhoy <- 9 * mean((2 * (min1YY + min2YY + min3YY + min4YY + 2 * G_YY^2 - 4 * G_YY) + 1 - rho_y)^2)
      var_rhorhox <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY)  + 1 - rho) * (2 * (min1XX + min2XX + min3XX + min4XX + 2 * G_XX^2 - 4 * G_XX) + 1 - rho_x))
      var_rhorhoy <- 9 * mean((4 * (g_xX + g_yY + G_XX * G_YY - G_XX - G_YY)  + 1 - rho) * (2 * (min1YY + min2YY + min3YY + min4YY + 2 * G_YY^2 - 4 * G_YY) + 1 - rho_y))
      var_rhoxrhoy <- 9 * mean((2 * (min1XX + min2XX + min3XX + min4XX + 2 * G_XX^2 - 4 * G_XX) + 1 - rho_x) * (2 * (min1YY + min2YY + min3YY + min4YY + 2 * G_YY^2 - 4 * G_YY) + 1 - rho_y))
      var_hat <- (var_rho - rho * (var_rhorhox / rho_x + var_rhorhoy / rho_y) + 0.25 * rho^2 * (var_rhox / rho_x^2 + var_rhoy / rho_y^2 + 2 * var_rhoxrhoy / (rho_x * rho_y))) / (rho_x * rho_y)
      # Variance under independence assumption
      var_hat_ind <- 1
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * rho_b / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * rho_b / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- Rhob_LRV(X, Y, rho, rho_x, rho_y, bandwidth = "Dehling")
      var_hat_ind <- Rhob_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * rho_b / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * rho_b / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(rho_b + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), rho_b + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(rho_b_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - rho_b^2)), tanh(rho_b_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - rho_b^2)))
    }
    res <- dplyr::tribble(~RhoB, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          rho_b, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
  if (method == "r"){
    r <- stats::cor(X, Y)
    r_fis <- atanh(r)
    mean_x <- mean(X)
    mean_y <- mean(Y)
    var_x <- stats::var(X)
    var_y <- stats::var(Y)
    if (isTRUE(IID)){
      C_XX_XX <- stats::var((X - mean_x)^2)
      C_YY_YY <- stats::var((Y - mean_y)^2)
      C_XX_YY <- stats::cov((X - mean_x)^2, (Y - mean_y)^2)
      C_XX_XY <- stats::cov((X - mean_x)^2, (Y - mean_y) * (X - mean_x))
      C_YY_XY <- stats::cov((Y - mean_y)^2, (Y - mean_y) * (X - mean_x))
      # Variance
      var_hat <- C_XX_YY / var_x / var_y + r^2 / 4 * (C_XX_XX / var_x^2 + 2 * C_XX_YY / var_x / var_y + C_YY_YY / var_y) - r * (C_XX_XY / sqrt(var_x)^3 / sqrt(var_y) + C_YY_XY / sqrt(var_y)^3 / sqrt(var_x))
      # Variance under independence assumption
      var_hat_ind <- 1
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * r / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * r / sqrt(var_hat)))
    } else if (isFALSE(IID)){
      var_hat <- Rho_LRV(X, Y)
      var_hat_ind <- Rho_ind_LRV(X, Y, bandwidth = "Dehling")
      p_val_ind <- stats::pnorm(-abs(sqrt(n) * r / sqrt(var_hat_ind)))
      p_val <- stats::pnorm(-abs(sqrt(n) * r / sqrt(var_hat)))
    } else stop("Please insert a valid option for the variable `IID`!", call. = FALSE)
    if (isFALSE(Fisher)){
      CI <- c(r + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n), r + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n))
    } else if (isTRUE(Fisher)){
      CI <- c(tanh(r_fis + stats::qnorm(alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - r^2)), tanh(r_fis + stats::qnorm(1 - alpha/2)*sqrt(var_hat)/sqrt(n)/(1 - r^2)))
    }
    res <- dplyr::tribble(~R, ~CI_lower, ~CI_upper, ~PValue, ~PValueIND,
                          r, CI[1], CI[2], p_val, p_val_ind)
    return(res)
  }
}

#' @keywords internal
Tau_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
if (length(X) != length(Y)){stop("X and Y must have equal length")}
n <- length(X)

# Determine bandwidth
if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
else stop("Please insert a valid bandwith calculation method")

# Calculate weights
h <- 1:(n-1)
w <- pmax(1 - abs(h) / (b + 1), 0)

# Calculate autocovariances in a vector with row = lag
x_autoc <- (n - 1) / n * stats::acf((rank(X) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # This is the acf of the demeaned grade. Therefore, demean = FALSE
y_autoc <- (n - 1) / n * stats::acf((rank(Y) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # This is the acf of the demeaned grade. Therefore, demean = FALSE

# Calculate estimator of LRV for tau under independence
Tau_ind_LRV <- 64 * sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))

return(Tau_ind_LRV)
}

#' @keywords internal
Rho_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
if (length(X) != length(Y)){stop("X and Y must have equal length")}
n <- length(X)

# Determine bandwidth
if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
else stop("Please insert a valid bandwith calculation method")

# Calculate weights
h <- 1:(n-1)
w <- pmax(1 - abs(h) / (b + 1), 0)

# Calculate autocovariances in a vector with row = lag
x_autoc <- stats::acf(X, plot = FALSE, type = "correlation", demean = TRUE, lag.max = n - 1)$acf
y_autoc <- stats::acf(Y, plot = FALSE, type = "correlation", demean = TRUE, lag.max = n - 1)$acf

# Calculate estimator of LRV for rho under independence
Rho_ind_LRV <- sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))

return(Rho_ind_LRV)
}

#' @keywords internal
Rhob_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
if (length(X) != length(Y)){stop("X and Y must have equal length")}
n <- length(X)

# Determine bandwidth
if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
else stop("Please insert a valid bandwith calculation method")

# Calculate weights
h <- 1:(n-1)
w <- pmax(1 - abs(h) / (b + 1), 0)

# Calculate autocovariances in a vector with row = lag
x_autoc <- stats::acf((rank(X) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf / stats::cov(X, X, method = "spearman") # This is the acf of the demeaned grade. Therefore, demean = FALSE
y_autoc <- stats::acf((rank(Y) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf / stats::cov(Y, Y, method = "spearman") # This is the acf of the demeaned grade. Therefore, demean = FALSE

# Calculate estimator of LRV for rho_b under independence
Rhob_ind_LRV <- sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))

return(Rhob_ind_LRV)
}

#' @keywords internal
Tau_LRV <- function(X, Y, kendall, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Define kernel realizations
  k_XY <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall

  # Calculate autocovariances in a vector with row = lag
  k_XY_autoc <- (n - 1) / n * stats::acf(k_XY, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE

  # Calculate estimator of LRV for tau
  Tau_LRV <- 4 * sum(k_XY_autoc[1], 2 * (w * k_XY_autoc[-1]))

  return(Tau_LRV)
}

#' @keywords internal
Rho_LRV <- function(X, Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Estimate values
  mean_x <- mean(X)
  mean_y <- mean(Y)
  sigma_xy <- (n - 1) / n * stats::cov(X, Y)
  var_x <- (n - 1) / n * stats::var(X)
  var_y <- (n - 1) / n * stats::var(Y)
  x_autoc <- (n - 1) / n * stats::acf(X, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y_autoc <- (n - 1) / n * stats::acf(Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2_autoc <- (n - 1) / n * stats::acf(X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y2_autoc <- (n - 1) / n * stats::acf(Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy_autoc <- (n - 1) / n * stats::acf(X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy_crossc <- (n - 1) / n * stats::ccf(X, Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xx2_crossc <- (n - 1) / n * stats::ccf(X, X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yx2_crossc <- (n - 1) / n * stats::ccf(Y, X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy2_crossc <- (n - 1) / n * stats::ccf(X, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yy2_crossc <- (n - 1) / n * stats::ccf(Y, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2y2_crossc <- (n - 1) / n * stats::ccf(X^2, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xxy_crossc <- (n - 1) / n * stats::ccf(X, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yxy_crossc <- (n - 1) / n * stats::ccf(Y, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2xy_crossc <- (n - 1) / n * stats::ccf(X^2, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y2xy_crossc <- (n - 1) / n * stats::ccf(Y^2, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf

  # Estimate Long-Run Variances
  x_LRVa <- sum(x_autoc[1], 2 * (w * x_autoc[-1]))
  y_LRVa <- sum(y_autoc[1], 2 * (w * y_autoc[-1]))
  x2_LRVa <- sum(x2_autoc[1], 2 * (w * x2_autoc[-1]))
  y2_LRVa <- sum(y2_autoc[1], 2 * (w * y2_autoc[-1]))
  xy_LRVa <- sum(xy_autoc[1], 2 * (w * xy_autoc[-1]))
  xy_LRVc <- sum(xy_crossc[1], 2 * (w * xy_crossc[-1]))
  xx2_LRVc <- sum(xx2_crossc[1], 2 * (w * xx2_crossc[-1]))
  yx2_LRVc <- sum(yx2_crossc[1], 2 * (w * yx2_crossc[-1]))
  xy2_LRVc <- sum(xy2_crossc[1], 2 * (w * xy2_crossc[-1]))
  yy2_LRVc <- sum(yy2_crossc[1], 2 * (w * yy2_crossc[-1]))
  x2y2_LRVc <- sum(x2y2_crossc[1], 2 * (w * x2y2_crossc[-1]))
  xxy_LRVc <- sum(xxy_crossc[1], 2 * (w * xxy_crossc[-1]))
  yxy_LRVc <- sum(yxy_crossc[1], 2 * (w * yxy_crossc[-1]))
  x2xy_LRVc <- sum(x2xy_crossc[1], 2 * (w * x2xy_crossc[-1]))
  y2xy_LRVc <- sum(y2xy_crossc[1], 2 * (w * y2xy_crossc[-1]))

  # Fill the matrix with long-Run variances
  Sigma <- diag(c(x_LRVa, y_LRVa, x2_LRVa, y2_LRVa, xy_LRVa))
  Sigma[upper.tri(Sigma)] <- c(xy_LRVc, xx2_LRVc, yx2_LRVc, xy2_LRVc, yy2_LRVc, x2y2_LRVc, xxy_LRVc, yxy_LRVc, x2xy_LRVc, y2xy_LRVc)
  Sigma <- as.matrix(Matrix::forceSymmetric(Sigma, uplo = "U"))

  # Create Jacobian matrices for Delta method
  A <- matrix(c(-2*mean_x, -mean_y, 0, 0, -mean_x, -2*mean_y, 1, 0, 0, 0, 0, 1, 0, 1, 0), nrow = 3)
  B <- c(-0.5 * sigma_xy / sqrt(var_y) * sqrt(var_x) ^ (-3), 1 / sqrt(var_y * var_x), -0.5 * sigma_xy / sqrt(var_x) * sqrt(var_y) ^ (-3))

  Rho_LRV <- B %*% A %*% Sigma %*% t(A) %*% B

  return(Rho_LRV)
}

#' @keywords internal
SRho_LRV <- function(X, Y, spearman, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Define kernel realizations
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  k_XY <- 4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - spearman

  # Calculate autocovariances in a vector with row = lag
  k_XY_autoc <- (n - 1) / n * stats::acf(k_XY, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE

  # Calculate estimator of LRV for srho
  SRho_LRV <- 9 * sum(k_XY_autoc[1], 2 * (w * k_XY_autoc[-1]))

  return(SRho_LRV)
}

#' @keywords internal
Gamma_LRV <- function(X, Y, kendall, tie_prob, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Define kernel realizations
  k_XY_tau <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall
  k_XY_tie <- x_eq(X) + y_eq(Y) - x_eq_y_eq(X, Y) - tie_prob

  # Calculate autocovariances in a vector with row = lag
  k_XY_tau_autoc <- (n - 1) / n * stats::acf(k_XY_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tau has mean 0. Therefore, demean = FALSE
  k_XY_tie_autoc <- (n - 1) / n * stats::acf(k_XY_tie, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tie has mean 0. Therefore, demean = FALSE
  k_XY_tautie_crossc <- (n - 1) / n * stats::ccf(k_XY_tau, k_XY_tie, plot = FALSE, type = "covariance", lag.max = n - 1)$acf

  # Calculate estimator of LRV for gamma
  sigma_tau_sq <- 4 * sum(k_XY_tau_autoc[1], 2 * (w * k_XY_tau_autoc[-1]))
  sigma_nu_sq <- 4 * sum(k_XY_tie_autoc[1], 2 * (w * k_XY_tie_autoc[-1]))
  sigma_taunu <- 4 * sum(c(sort(w), 1, w) * k_XY_tautie_crossc)
  gamma <- kendall/(1-tie_prob)

  Gamma_LRV <- (sigma_tau_sq + gamma^2 * sigma_nu_sq + 2 * gamma * sigma_taunu) / (1 - tie_prob)^2

  return(Gamma_LRV)
}

#' @keywords internal
TauB_LRV <- function(X, Y, kendall, kendall_X, kendall_Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Define kernel realizations
  k_XY_tau <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall
  k_X_tau <- 1 - x_eq(X) - kendall_X
  k_Y_tau <- 1 - y_eq(Y) - kendall_Y

  # Calculate autocovariances in a vector with row = lag
  k_XY_tau_autoc <- (n - 1) / n * stats::acf(k_XY_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tau has mean 0. Therefore, demean = FALSE
  k_X_tau_autoc <- (n - 1) / n * stats::acf(k_X_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_X_tau has mean 0. Therefore, demean = FALSE
  k_Y_tau_autoc <- (n - 1) / n * stats::acf(k_Y_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_Y_tau has mean 0. Therefore, demean = FALSE
  k_X_tau_crossc <- (n - 1) / n * stats::ccf(k_XY_tau, k_X_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_Y_tau_crossc <- (n - 1) / n * stats::ccf(k_XY_tau, k_Y_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_XY_tautau_crossc <- (n - 1) / n * stats::ccf(k_X_tau, k_Y_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf

  # Calculate estimator of LRV for taub
  sigma_tau_sq <- 4 * sum(k_XY_tau_autoc[1], 2 * (w * k_XY_tau_autoc[-1]))
  sigma_tauX_sq <- 4 * sum(k_X_tau_autoc[1], 2 * (w * k_X_tau_autoc[-1]))
  sigma_tauY_sq <- 4 * sum(k_Y_tau_autoc[1], 2 * (w * k_Y_tau_autoc[-1]))
  sigma_tautauX <- 4 * sum(c(sort(w), 1, w) * k_X_tau_crossc)
  sigma_tautauY <- 4 * sum(c(sort(w), 1, w) * k_Y_tau_crossc)
  sigma_tauXtauY <- 4 * sum(c(sort(w), 1, w) * k_XY_tautau_crossc)

  TauB_LRV <- (sigma_tau_sq - kendall * (sigma_tautauX / kendall_X - sigma_tautauY / kendall_Y) + kendall^2 / 4 * (sigma_tauX_sq / kendall_X^2 + sigma_tauY_sq / kendall_Y^2 + (2 * sigma_tauXtauY) / kendall_X / kendall_Y)) / (kendall_X * kendall_Y)
  return(TauB_LRV)
}

#' @keywords internal
Rhob_LRV <- function(X, Y, spearman, spearman_X, spearman_Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)

  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")

  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)

  # Define kernel realizations
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  k_XY_rho <- 4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - spearman
  k_X_rho <- 2 * (f_x1(X) + f_x2(X) + f_x3(X) + f_x4(X) + 2 * G_X(X)^2 - 4 * G_X(X)) + 1 - spearman_X
  k_Y_rho <- 2 * (f_y1(Y) + f_y2(Y) + f_y3(Y) + f_y4(Y) + 2 * G_Y(Y)^2 - 4 * G_Y(Y)) + 1 - spearman_Y

  # Calculate autocovariances in a vector with row = lag
  k_XY_rho_autoc <- (n - 1) / n * stats::acf(k_XY_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE
  k_X_rho_autoc <- (n - 1) / n * stats::acf(k_X_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_X_tie has mean 0. Therefore, demean = FALSE
  k_Y_rho_autoc <- (n - 1) / n * stats::acf(k_Y_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_Y_tie has mean 0. Therefore, demean = FALSE
  k_X_rho_crossc <- (n - 1) / n * stats::ccf(k_XY_rho, k_X_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_Y_rho_crossc <- (n - 1) / n * stats::ccf(k_XY_rho, k_Y_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_XY_rho_crossc <- (n - 1) / n * stats::ccf(k_X_rho, k_Y_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf

  # Calculate estimator of LRV for srho
  sigma_rho_sq <- 9 * sum(k_XY_rho_autoc[1], 2 * (w * k_XY_rho_autoc[-1]))
  sigma_rhoX_sq <- 9 * sum(k_X_rho_autoc[1], 2 * (w * k_X_rho_autoc[-1]))
  sigma_rhoY_sq <- 9 * sum(k_Y_rho_autoc[1], 2 * (w * k_Y_rho_autoc[-1]))
  sigma_rhorhoX <- 9 * sum(c(sort(w), 1, w) * k_X_rho_crossc)
  sigma_rhorhoY <- 9 * sum(c(sort(w), 1, w) * k_Y_rho_crossc)
  sigma_rhoXrhoY <- 9 * sum(c(sort(w), 1, w) * k_XY_rho_crossc)

  Rhob_LRV <- (sigma_rho_sq - spearman * (sigma_rhorhoX / spearman_X - sigma_rhorhoY / spearman_Y) + spearman^2 / 4 * (sigma_rhoX_sq / spearman_X^2 + sigma_rhoY_sq / spearman_Y^2 + (2 * sigma_rhoXrhoY) / spearman_Y / spearman_X)) / (spearman_X * spearman_Y)

  return(Rhob_LRV)
}

#' @keywords internal
G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
#' @keywords internal
G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
#' @keywords internal
F_X <- Vectorize(function(x_val) mean(X <= x_val))
#' @keywords internal
F_X_ <- Vectorize(function(x_val) mean(X < x_val))
#' @keywords internal
G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
#' @keywords internal
F_Y <- Vectorize(function(y_val) mean(Y <= y_val))
#' @keywords internal
F_Y_ <- Vectorize(function(y_val) mean(Y < y_val))
#' @keywords internal
g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
#' @keywords internal
g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
#' @keywords internal
f_x1 <- Vectorize(function(x_val) mean(min(F_X(x_val), F_X(X))))
#' @keywords internal
f_x2 <- Vectorize(function(x_val) mean(min(F_X_(x_val), F_X(X))))
#' @keywords internal
f_x3 <- Vectorize(function(x_val) mean(min(F_X(x_val), F_X_(X))))
#' @keywords internal
f_x4 <- Vectorize(function(x_val) mean(min(F_X_(x_val), F_X_(X))))
#' @keywords internal
f_y1 <- Vectorize(function(y_val) mean(min(F_Y(y_val), F_Y(Y))))
#' @keywords internal
f_y2 <- Vectorize(function(y_val) mean(min(F_Y_(y_val), F_Y(Y))))
#' @keywords internal
f_y3 <- Vectorize(function(y_val) mean(min(F_Y(y_val), F_Y_(Y))))
#' @keywords internal
f_y4 <- Vectorize(function(y_val) mean(min(F_Y_(y_val), F_Y_(Y))))
#' @keywords internal
x_eq <- Vectorize(function(x_val) mean(X == x_val))
#' @keywords internal
y_eq <- Vectorize(function(y_val) mean(Y == y_val))
#' @keywords internal
x_eq_y_eq <- Vectorize(function(x_val, y_val) mean(X == x_val & Y == y_val))
#' @keywords internal
min1X <- Vectorize(function(x_val) mean(pmin(F_XX, F_X(x_val))))
#' @keywords internal
min2X <- Vectorize(function(x_val) mean(pmin(F_XX, F_X_(x_val))))
#' @keywords internal
min3X <- Vectorize(function(x_val) mean(pmin(F_X_X, F_X(x_val))))
#' @keywords internal
min4X <- Vectorize(function(x_val) mean(pmin(F_X_X, F_X_(x_val))))
#' @keywords internal
min1Y <- Vectorize(function(y_val) mean(pmin(F_YY, F_Y(y_val))))
#' @keywords internal
min2Y <- Vectorize(function(y_val) mean(pmin(F_YY, F_Y_(y_val))))
#' @keywords internal
min3Y <- Vectorize(function(y_val) mean(pmin(F_Y_Y, F_Y(y_val))))
#' @keywords internal
min4Y <- Vectorize(function(y_val) mean(pmin(F_Y_Y, F_Y_(y_val))))

