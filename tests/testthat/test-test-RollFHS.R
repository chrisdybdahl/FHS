test_that("RollFHS runs without errors and returns expected structure", {
  # Generate sample data
  set.seed(12)
  #dates <- seq.Date(from = as.Date("2023-01-01"), by = "day", length.out = 1000)
  #df <- data.frame(Date = dates, Return = rnorm(1000))
  df <- read.csv("C:/Users/chris/RStudioProjects/FHS/data/clean_returns.csv", row.names = 1)
  asset_df <- df["DEPMc1"]
  asset_df$Date <- as.Date(rownames(asset_df))
  colnames(asset_df) <- c("Return", "Date")
  rownames(asset_df) <- NULL


  # Parameters
  c <- 0.05
  n <- 989
  m <- 250
  r <- 10
  model <- "gjrGARCH"
  dist <- "sstd"
  type <- "EWMA"
  verbose <- 1
  solver <- "hybrid"
  solver.control <- list(trace = 0)

  # Call the function
  result <- RollFHS(asset_df, c = c, n = n, m = m, r = r,
                    model = model, dist = dist, type = type,
                    solver = solver, solver.control = solver.control,
                    verbose = verbose)

  # Check the result is of the correct type
  expect_s3_class(result, "xts")

  # Check the result has the correct dimensions
  expect_equal(nrow(result), n)  # Ensure it has 'n' rows
  expect_equal(ncol(result), 2)  # Should have columns for VaR and ES

  # Validate the data range
  expect_true(all(result$VaR >= 0))  # VaR should be negative
  expect_true(all(result$ES >= 0))   # ES should be negative
})
