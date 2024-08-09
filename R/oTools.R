#' Convert long format data to matrix (wide) format
#'
#' This function takes a long format dataset (each row is a different
#' participant and a different trial) and converts it to a matrix format (each
#' row is a participant and each column is a trial). By default, additional
#' columns and rows are added with summary statistics. An additional column is
#' added at the end for the sum score of each participant. An additional three
#' rows are added for mean, variance, and item-to-rest correlation of each
#' trial.
#'
#' @param data A dataframe in the long format
#' @param idCol The name of the column that contains the participant id
#' @param corrCol The name of the column that contains the score of the trial
#' @param trialNCol The name of the column that contains the trial number
#' @param summaryStats A boolean indicating whether to add summary statistics
#' @return The data in the wide format
#' @export
long2Matrix <- function(data, idCol, corrCol, trialNCol, summaryStats = TRUE) {
  # Get subjects and trials
  sbjs <- data[idCol] |>
    unique() |>
    pull() |>
    sort()
  trials <- unique(data[trialNCol]) |>
    pull() |>
    as.character()

  # Create a tibble with subject id row
  taskMatrix <- tibble(!!idCol := sbjs) # nolint

  for (trial in trials) {
    trialData <- data[data[trialNCol] == trial, c(idCol, corrCol)]
    # Rename corr to the trial
    names(trialData)[2] <- trial

    taskMatrix <- taskMatrix |>
      left_join(
        trialData,
        by = idCol,
        relationship = "one-to-one"
      )
  }

  # Set the SbjID column as row names
  taskMatrix <- column_to_rownames(taskMatrix, var = idCol)

  if (summaryStats) {
    # Save raw matrix for later
    rawInfo <- taskMatrix

    # Add mean row
    taskMatrix <- rbind(taskMatrix, mean = apply(taskMatrix, 2, mean))

    # Add participant mean row
    taskMatrix <- cbind(taskMatrix, Sum = apply(taskMatrix, 1, sum))

    # Add correlation row
    totals <- taskMatrix[seq_len(nrow(taskMatrix) - 1), "Sum"]
    tmp <- c()
    for (col in seq_len(ncol(rawInfo))) {
      tmp <- c(tmp, cor(rawInfo[col], totals - rawInfo[col]))
    }
    taskMatrix["cor", ] <- c(tmp, 1)

    # Add variance row
    sbjs <- as.character(t(unique(data[idCol])))
    taskMatrix <- rbind(taskMatrix, var = apply(taskMatrix[sbjs, ], 2, var))
  }

  return(taskMatrix)
}


#' Create parcels given a long format dataset
#'
#' This function takes a long format dataset (each row is a different
#' participant and a different trial) and creates parcels of the data. If the
#' randomize parameter is set to TRUE, the trials will be split between parcels
#' randomly. Otherwise, the trials will be split between parcels in order (e.g.
#' 2 parcels will use even/odd trial numbers). The parcelsr are returned as a
#' list of tibbles.
#'
#' @param data A dataframe in the long format
#' @param trialNCol The name of the column that contains the trial number
#' @param nParcels The number of parcels to create, default 2
#' @param randomize A boolean indicating whether to randomize the trials
#' @return A list of tibbles with the parcels
#' @export
makeParcels <- function(data, trialNCol, nParcels = 2, randomize = FALSE) {
  # Randomize trialNs if needed
  if (randomize) {
    trialNumbers <- data |>
      select(!!rlang::sym(trialNCol)) |>
      unique() |>
      pull()

    # Create a mapping from original and randomized
    trialMap <- tibble(
      !!trialNCol := trialNumbers, # nolint
      randomized = sample(trialNumbers)
    )

    # Remap all the trial numbers
    data <- data |>
      left_join(trialMap, by = trialNCol) |>
      mutate(!!rlang::sym(trialNCol) := .data$randomized) |> # nolint
      select(-.data$randomized)
  }

  parcels <- list()

  for (i in 1:nParcels) {
    parcels[[i]] <- data |>
      filter(!!rlang::sym(trialNCol) %% nParcels == (i - 1))
  }

  return(parcels)
}

#' Calculate reliability given a wide format dataset (matrix).
#'
#' This function takes a a wide format dataset (matrix form) from this package,
#' which includes extra rows and columns for summary statistics and calculates
#' reliability. If the matrix does not have the extra summary statistics, you
#' should just use splitHalf from the psych package.
#'
#' @param mat A matrix of data
#' @return splitHalf object
#' @export
matrix2Rel <- function(mat) {
  tmp <- mat[, !is.na(mat["cor", ])]
  tmp <- tmp[seq_len(nrow(tmp) - 3), seq_len(ncol(tmp) - 1)]

  return(splitHalf(tmp, check.keys = FALSE))
}

#' Calculate the disattenuated correlation
#'
#' @param rxy The correlation between two variables
#' @param rxx The reliability of the first variable
#' @param ryy The reliability of the second variable
#' @return The disattenuated correlation
#' @export
disattCor <- function(rxy, rxx, ryy) {
  return(rxy / sqrt(rxx * ryy))
}

#' Calculate aggregate equal weighted reliability
#'
#' Given a dataframe of data where each column is a measure and each row is a
#' participant alongside an array of reliabilities (one for each measures),
#' calculate the aggregate reliability with equal weighting for each measure.
#'
#' @param data A dataframe of data
#' @param rel An array of reliabilities
#' @return The aggregate reliability
#' @export
aggRel <- function(data, rel) {
  # Calculate correlation matrix
  corMat <- cor(data, use = "complete.obs")
  diag(corMat) <- NA

  numer <- sum(rel) + sum(corMat, na.rm = TRUE)
  denom <- length(rel) + sum(corMat, na.rm = TRUE)

  return(numer / denom)
}

#' Create a correlation matrix pairs plot.
#'
#' Create a pairs plot for a set of measures. The diagonal typically includes
#' the reliabilities and a KDE of the distribution of scores across
#' participants. The bottom triangle plots scatter plots with a fit line. The
#' upper triangle indicates the correlation between the two measures with BF
#' and a pie chart representing the relative evidence between a positive
#' correlation and a point null.
#'
#' @param data A dataframe of data (each column is a measure, each row is a participant)
#' @param reliability An array of reliabilities for each measure
#' @param nullInterval The null interval for the correlation BF
#' @param rscale The scale for the correlation BF
#' @param relSymbol The symbol to use for reliability, can be an array of symbols for each measure
#' @param draw_dist A boolean indicating whether to draw the distribution on the diagonal
#' @param showN A boolean indicating whether to show the number of observations
#' @param completedOnly A boolean indicating whether to only use complete cases
#' @return A ggplot object
#' @export
#'
corMatrixPlot <- function(
    data,
    reliability,
    nullInterval = NULL,
    rscale = NULL,
    relSymbol = "r",
    draw_dist = FALSE,
    showN = FALSE,
    completedOnly = TRUE) {
  if (completedOnly) {
    data <- data[complete.cases(data), ]
  }

  # Setup test labels and reliability symbols
  testLabels <- names(data)
  nTests <- length(testLabels)
  if (length(relSymbol) == 1) {
    relSymbol <- rep(relSymbol, times = nTests)
  }

  # Create base plot matrix
  plots <- list()
  for (i in 1:(nTests * nTests)) {
    plots[[i]] <- GGally::ggally_text(paste("Plot #", i, sep = ""))
  }
  corMatPlot <- GGally::ggmatrix(
    plots = plots,
    nrow = nTests,
    ncol = nTests,
    xAxisLabels = testLabels,
    yAxisLabels = testLabels,
  ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(color = "black"),
      strip.background = element_rect(
        fill = "white"
      )
    )

  # Edit each cell
  for (row in 1:corMatPlot$nrow) {
    for (col in 1:corMatPlot$ncol) {
      if (row == col) { # Diagonal reliability
        rel <- round(reliability[row], digits = 2)
        rel <- str_replace(rel, "0.", ".")
        if (str_length(rel) == 2) {
          rel <- paste0(rel, "0")
        }

        corMatPlot[row, row] <- ggplot(
          data,
          aes_string(x = testLabels[row])
        ) +
          geom_density(
            color = ifelse(draw_dist, "gray", NA)
          ) +
          ggpp::annotate(
            "text_npc",
            npcx = "center",
            npcy = "middle",
            label = paste0(
              relSymbol[row], " = ", rel,
              "\n\n", testLabels[row]
            )
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
          )
      } else if (row < col) { # Off diagonal correlations
        x <- dplyr::pull(data[, col])
        y <- dplyr::pull(data[, row])

        # Calculate BF
        corBF <- BayesFactor::correlationBF(
          x, y,
          nullInterval = nullInterval,
          rscale = rscale
        )

        # Correlation with BFs
        if (length(corBF) == 2) {
          bfSamples <- BayesFactor::posterior(
            corBF[2],
            iterations = 10000,
            progress = FALSE
          )
          rho <- summary(bfSamples)$statistics[1, 1]
          bf <- BayesFactor::extractBF(corBF[2])$bf
        } else {
          bfSamples <- BayesFactor::posterior(
            corBF,
            iterations = 10000,
            progress = FALSE
          )
          rho <- summary(bfSamples)$statistics[1, 1]
          bf <- BayesFactor::extractBF(corBF)$bf
        }

        # Make the text
        rhoText <- str_replace(
          as.character(round(rho, digits = 2)), "0.", "."
        )
        BFText <- as.character(round(bf, digits = 2))

        if (showN) {
          # Count number of full observations
          nObs <- sum(complete.cases(data[, c(row, col)]))
          corText <- paste0(
            "atop(", paste0("r(", nObs, ') == "', rhoText), '",',
            paste0("\nBF[+0] == ", BFText), ")"
          )
        } else {
          corText <- paste0(
            "atop(", paste0('r == "', rhoText), '",',
            paste0("\nBF[+0] == ", BFText), ")"
          )
        }

        # Place text in the cell
        pieData <- data.frame(
          Value = c(bf, 1),
          Hypo = c("H1", "H0")
        )
        pieData$Fraction <- pieData$Value / sum(pieData$Value)
        pieData$YMax <- cumsum(pieData$Fraction)
        pieData$YMin <- c(0, head(pieData$YMax, n = -1))

        corMatPlot[row, col] <- ggplot(pieData) +
          aes(
            ymax = .data$YMax,
            ymin = .data$YMin, xmax = 4, xmin = 1, fill = .data$Hypo # nolint
          ) +
          geom_rect(color = "white", alpha = .5) +
          coord_polar("y", start = (1 / (bf + 1)) * pi + pi) +
          scale_fill_manual(values = c("coral", "aquamarine3")) +
          xlim(c(1, 4)) +
          annotate(
            "text",
            x = 1,
            y = 1,
            label = corText,
            parse = TRUE
          ) +
          theme(
            panel.background = element_rect(color = "white"),
            axis.text = element_blank(),
            axis.ticks = element_blank()
          )


        # Create scatter plot for the lower diagonal cell
        corMatPlot[col, row] <- tryCatch(
          {
            # Calculate line of best fit
            plotLM <- rstanarm::stan_glm(
              as.formula(paste0(testLabels[col], "~", testLabels[row])),
              data = data,
              refresh = 0,
            )
            drawsLines <- as.data.frame(as.matrix(plotLM))
            colnames(drawsLines) <- c("intercept", "slope", "sigma")
            ggplot(data, aes_string(
              x = testLabels[row],
              y = testLabels[col]
            )) +
              geom_abline(
                data = drawsLines,
                aes(intercept = .data$intercept, slope = .data$slope), # nolint
                color = "gray",
                size = 0.1,
                alpha = 0.1
              ) +
              geom_abline(
                intercept = coef(plotLM)[1],
                slope = coef(plotLM)[2]
              ) +
              geom_point(alpha = 0.3)
            # Aspect ratio issue https://github.com/ggobi/ggally/issues/415
          },
          error = function(e) {
            ggplot(data, aes_string(
              x = testLabels[row],
              y = testLabels[col]
            )) +
              geom_point(alpha = 0.3)
          }
        )
      }
    }
  }

  return(corMatPlot)
}
