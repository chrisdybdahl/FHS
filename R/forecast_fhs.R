#' Rolling Forecasts for Filtered Historical Simulation (FHS)
#'
#' @useDynLib FHS
#' @import rugarch
#' @import progress
#' @importFrom xts xts
#' @export
RollFHS <- function(
    data,
    c,
    n,
    m,
    p = 1,
    q = 1,
    r = 1,
    model = "sGARCH",
    dist = "norm",
    type = c("EWMA", "GARCH"),
    lambda = 0.94,
    nboot = 10000,
    verbose = 0,
    solver = "hybrid",
    solver.control = list(tol = 1e-7),
    ...
) {
  type <- match.arg(type)
  if (type == "EWMA" && (length(lambda) != 1 || !is.numeric(lambda) || lambda <= 0 || lambda >= 1)) {
    stop("The argument 'lambda' must be a single numeric value between 0 and 1 for EWMA.")
  }
  if (is.null(nboot) || !is.numeric(nboot) || nboot <= 0) {
    stop("The argument 'nboot' must be a positive integer.")
  }

  if (verbose > 0) {
    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsed || Estimated time remaining: :eta]",
     total = n,
     complete = "=",
     incomplete = "-",
     current = ">",
     show_after = 0.2,
     clear = TRUE
     )
  }

  df <- tail(data, n + m)
  data_xts <- xts::xts(df$Return, order.by = df$Date)

  spec_unrestricted <- rugarch::ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(garchOrder = c(p, q), model = model),
    distribution.model = dist
  )

  # Initialize storage for results
  var <- matrix(NA, nrow = n, ncol = length(c), dimnames = list(NULL, paste0("VaR_", c)))
  es <- matrix(NA, nrow = n, ncol = length(c), dimnames = list(NULL, paste0("ES_", c)))

  sigma_t <- NA
  params <- NA
  for (i in 1:n) {
    if (verbose > 0) pb$tick()
    window_start <- i
    window_end <- m + i - 1
    window <- data_xts[window_start:window_end]

    # Reparametrize for every r period
    if (i %% r == 0 || i == 1) {
      # Estimate conditional volatility
      if (type == "GARCH") {
        fit <- rugarch::ugarchfit(
          spec = spec_unrestricted,
          data = df,
          solver = solver,
          solver.control = solver.control
        )

        params <- rugarch::coef(fit)
      }
    }

    if (type == "GARCH") {
      spec_restricted <- rugarch::ugarchspec(
        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
        variance.model = list(garchOrder = c(p, q), model = model),
        distribution.model = dist,
        fixed.pars = list(params)
      )

      fit <- rugarch::ugarchfit(
        spec = spec_restricted,
        data = df,
        solver = solver,
        solver.control = solver.control
      )

      cond_sigma <- as.matrix(rugarch::sigma(fit), ncol = length(c))
      cond_sigma_ahead <- as.matrix(rugarch::sigma(rugarch::ugarchforecast(fit, n.ahead = 1)), ncol = length(c))

    } else if (type == "EWMA") {

      vol <- stats::var(window)
      cond_var <- ewma(window, vol, lambda = lambda)
      cond_var_ahead <-  lambda * cond_var[m] + (1 - lambda) * window[m]^2
      cond_sigma <- sqrt(cond_var)
      cond_sigma_ahead <- sqrt(cond_var_ahead)
    }

    # Scale returns
    scaled_returns <- as.numeric(window) / rep(cond_sigma, length.out = length(window))

    # Bootstrap sampling
    boot_samples <- sample(scaled_returns, size = nboot, replace = TRUE)


    # Rescale bootstrapped returns
    boot_losses <- boot_samples * rep(cond_sigma_ahead, length.out = length(boot_samples))

    # Compute VaR and ES
    var[i,] <- stats::quantile(boot_losses, probs = c)
    es[i,] <- mean(boot_losses[boot_losses < var[i,]])
  }

  if (verbose > 0) {
    pb$terminate()
  }

  # Create xts objects for VaR and ES
  dates <- tail(data$Date, n)
  VaR <- xts::xts(-var, order.by = dates)
  ES <- xts::xts(-es, order.by = dates)
  results_xts <- xts::merge.xts(VaR, ES)

  if (any(is.na(results_xts))) {
    warning("There were NA values, carry forward last non-NA")
    results_xts <- zoo::na.locf(results_xts, na.rm = FALSE)
  }

  # Return results
  return(results_xts)
}

