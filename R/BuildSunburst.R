#' Create a Sunburst Plot for Treatment Pathways
#'
#' This function creates an interactive sunburst plot to visualize treatment pathways
#' from CohortPathway analysis results. The plot shows the hierarchical structure of
#' treatment sequences, with each ring representing a step in the pathway.
#'
#' @param cpResults A list containing the results from CohortPathway analysis.
#'        Must include 'pathwaysAnalysisPathsData' and 'isCombo' data frames.
#' @param generationSet A data frame containing information about cohorts
#'        that will be used to generate descriptive names for the events in the diagram and target names of the plot
#' @param nPaths Integer specifying the maximum number of steps to include in the plot.
#' @param minCount Integer specifying the minimum count value for a path to be included.
#' @return An HTML widget object containing the interactive sunburst plot.
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' library(CohortPathway)
#' sunburstPlot <- CohortPathways::createPathwaySunburst(cpResults, generationSet)
#' }
createPathwaySunburst <- function(
    cpResults,
    generationSet,
    nPaths = 3,
    minCount = 5,
    plotWidth = "80%",
    plotHeight = 600) {
  
  rlang::check_installed("sunburstR")
  rlang::check_installed("htmlwidgets")
  rlang::check_installed("d3r")
  
  # Input validation
  checkmate::assertList(
    cpResults,
    min.len = 7,
    types = "data.frame"
  )
  
  # Extract "isCombo" data frame from "cpResults" list and restrict paths to maximum value (nPaths)
  isCombo <- purrr::pluck(cpResults, "isCombo") |>
    dplyr::filter(numberOfEvents <= nPaths)
  
  checkmate::assertDataFrame(x = isCombo, min.rows = 1, min.cols = 1)
  
  # Get event names
  eventNames <- .splitEvenToPowers(
    df = isCombo, 
    generationSet = generationSet,
    cpResults = cpResults
  )
  
  # Filter rows to those with count above minCount value (default=5)
  pathsData <- cpResults$pathwaysAnalysisPathsData|>
    dplyr::select(-c(pathwayAnalysisGenerationId, targetCohortId)) 
  # |>
  #   dplyr::filter(countValue > minCount)
  
  # Remove columns with all NA values
  pathsData <- pathsData[, colSums(!is.na(pathsData)) > 0]
  
  # Stop function if the selected number of paths exceeds the number of paths in the analysis data (cpResults$pathwaysAnalysisPathsData)
  if (nPaths > ncol(pathsData)-1) {
    
    stop(paste0
         (
        "Error: Number of selected paths (", 
        nPaths, 
        ") exceeds the number of paths in the analysis data (", 
        ncol(pathsData)-1, ").",
        "Please select a value of ",
        ncol(pathsData)-1,
        " or greater."
     )
    )
  }
  
  # Generate "step" column names
  colStepNames <- c(
    paste0("step", 1:nPaths),
    "countValue"
  )
  
  # Subset the data frame selecting the columns from above
  pathsData <- pathsData[, colStepNames]
  
  # Select all "step" column names
  stepCols <- grep("^step", names(pathsData), value = TRUE)
  
  # Loop over "step" columns and join
  stepNames <- purrr::map_dfc(
    stepCols,
    function(colname) {

      joinColSuffix <- gsub("step", "", colname)

      pathsData |>
        dplyr::select(tidyselect::all_of(colname)) |>
        dplyr::left_join(eventNames, by = setNames("comboId", colname)) |>
        dplyr::transmute(
          !!paste0("pathName", joinColSuffix) := pathName
        )
    }
  )
  
  # # Combine original data with new step name columns
  pathsDataFinal <- dplyr::bind_cols(pathsData, stepNames) |>
    dplyr::select(!dplyr::contains("step")) |>
    dplyr::filter(countValue > minCount)
  
  # Convert tabular data to JSON
  pathsDataJson <- d3r::d3_nest(
    pathsDataFinal, 
    value_cols = "countValue"
  )
  
  # Create sunburst plot
  sunburstPlot <- sunburstR::sunburst(
    data = pathsDataJson, 
    width = plotWidth, 
    height = plotHeight, 
    valueField = "countValue",
    legend = list(w = 490, h = 50, r = 100, s = 5),
    count = TRUE
  )
  
  return(sunburstPlot)
}


.prepareEventNames <- function(generationSet, cpResults) {
  
  event_names <- purrr::pluck(
    cpResults, "pathwayAnalysisCodesLong"
  ) |>
    dplyr::select(.data$code, cohortId = .data$eventCohortId) |>
    dplyr::distinct() |>
    dplyr::inner_join(
      generationSet |>
        dplyr::select(.data$cohortId, .data$cohortName),
      by = dplyr::join_by(cohortId)
    ) |>
    dplyr::group_by(.data$code) |>
    dplyr::reframe(
      combination = paste(.data$cohortName, collapse = " & ")
    )
  
  return(event_names)
}


# Function to split an even number into a chosen number of power-of-two summands.
.splitEvenToPowers <- function(df, generationSet, cpResults) {
  
  # Set variables
  comboId <- df$comboId
  isCombo <- df$isCombo
  numberOfEvents <- df$numberOfEvents
  
  # Create empty data frame
  data_frame <- data.frame(
    comboId = integer(), 
    isCombo = integer(),
    numberOfEvents = integer(),
    splitNumbers = character()
  )
  
  # Loop through rows of the data frame
  for (i in 1:nrow(df)) {
    
    # If isCombo is 0, simply return the input number without splitting
    if (isCombo[i] == 0) {
      
      # Create data frame
      data_frame_no_combos <- data.frame(
        comboId = comboId[i],
        numberOfEvents = numberOfEvents[i],
        isCombo = isCombo[i],
        splitNumbers = as.character(comboId[i])
      )
      
      # Append rows to data frame
      data_frame <- rbind(data_frame, data_frame_no_combos)
      
      
    } else {
      
      # Check if the number is even
      if (comboId[i] %% 2 != 0) {
        stop("Please input an even number!")
      }
      
      # Get the binary representation as a vector (least-significant bit first)
      bits <- as.integer(intToBits(comboId[i]))
      
      # Identify positions where the bit is 1 (subtract one for exponent)
      exponents <- which(bits == 1) - 1
      
      # For even numbers, ignore the 2^0 component (which equals 1)
      exponents <- exponents[exponents != 0]
      
      # Compute the corresponding powers of two from the exponents
      powers <- 2^exponents
      
      # Sort the summands in descending order (largest first)
      currentParts <- sort(powers, decreasing = TRUE)
      
      # Define minimal and maximum possible parts for a valid split
      minParts <- length(currentParts)
      maxParts <- comboId[i] / 2  # since the smallest summand allowed is 2
      
      if (numberOfEvents[i] < minParts) {
        stop(paste("The minimal splitting has", minParts, "numberOfEvents. Cannot merge components further."))
      }
      
      if (numberOfEvents[i] > maxParts) {
        stop(paste("The maximal splitting into powers >1 is", maxParts, "numberOfEvents."))
      }
      
      # Iteratively split the summands until the desired number of parts is reached
      while (length(currentParts) < numberOfEvents[i]) {
        
        # We cannot split further if every summand is 2
        if (all(currentParts == 2)) {
          stop("Cannot further split without producing ones.")
        }
        
        # Choose the largest summand that is greater than 2
        candidates <- currentParts[currentParts > 2]
        idx <- which(currentParts == max(candidates))[1]
        valueToSplit <- currentParts[idx]
        
        # Replace the chosen summand with two equal halves
        currentParts <- currentParts[-idx]   # Remove the selected summand
        currentParts <- c(currentParts, valueToSplit / 2, valueToSplit / 2)
        
        # Resort in descending order for consistency.
        currentParts <- sort(currentParts, decreasing = TRUE)
      }
      
      # Convert the numeric vector to a single string, with elements separated by commas
      resultString <- paste(currentParts, collapse = ",")
      
      # Create data frame
      data_frame_with_combos <- data.frame(
        comboId = comboId[i],
        numberOfEvents = numberOfEvents[i],
        isCombo = isCombo[i],
        splitNumbers = resultString
      )
      
      # Append values to data frame
      data_frame <- rbind(data_frame, data_frame_with_combos)
      
    }
  }
  
  # Maximum number of columns to create
  maxCols <- max(unique(data_frame$numberOfEvents))
  
  # Create a vector of column names
  newColnames <- paste0("eventCohortCode_", 1:maxCols)
  
  # Convert the vector of column names to a data frame
  df <- setNames(data.frame(matrix(ncol = length(newColnames), nrow = 0)), newColnames)
  
  # Add NA in all rows (no. of rows is equal to the length of the original data frame)
  df[nrow(data_frame),] <- NA
  
  # Bind empty data frame columns back to the original data frame
  data_frame <- cbind(data_frame, df)
  
  # Split "splitNumbers" column to multiple columns
  data_frame <- tidyr::separate(
    data_frame, 
    col = splitNumbers,
    into = newColnames, 
    sep = ",", 
    convert = TRUE, 
    remove = FALSE
  )
  
  # Get event cohort ids of and code (comboId)
  eventCohortIdAndCode <- cpResults[["pathwayAnalysisCodesLong"]] |> 
    dplyr::filter(isCombo == 0) |> 
    dplyr::select(c(eventCohortId, code))
  
  # Extract "eventCohortCode_" column names
  eventCohortCodeCols <- grep("^eventCohortCode_", names(data_frame), value = TRUE)
  
  # Join columns dynamically
  eventCohortIds <- purrr::map_dfc(
    eventCohortCodeCols, 
    function(colname) {
      
      joinColSuffix <- gsub("eventCohortCode", "", colname)  # Extract the number, e.g. '1'
      
      data_frame |>
        dplyr::select(tidyselect::all_of(colname)) |>
        dplyr::left_join(eventCohortIdAndCode, by = setNames("code", colname)) |>
        dplyr::mutate(
          !!paste0("eventCohortId", joinColSuffix) := eventCohortId, .keep = "unused"
        )
   }
  )
  
  # Extract event cohort names
  cohortDefinitionSet <- generationSet |>
    dplyr::select(c(cohortId, cohortName))
  
  # Extract "eventCohortId_" column names
  eventCohortIdCols <- grep("^eventCohortId_", names(eventCohortIds), value = TRUE)
  
  # Join columns dynamically
  eventCohortNames <- purrr::map_dfc(
    eventCohortIdCols,
    function(colname) {
      
      joinColSuffix <- gsub("eventCohortId_", "", colname)
      
      eventCohortIds |>
        dplyr::select(tidyselect::all_of(colname)) |>
        dplyr::left_join(cohortDefinitionSet, by = setNames("cohortId", colname)) |>
        dplyr::transmute(
          !!paste0("eventCohortName_", joinColSuffix) := cohortName
        )
    }
  )
  
  # Combine the original data with the newly joined names
  data_frame <- dplyr::bind_cols(
    data_frame |> dplyr::select(comboId), 
    eventCohortNames
  )
  
  # Columns to concatenate
  colsToConcat <- grep("^eventCohortName_", names(data_frame), value = TRUE)
  
  # Unite columns
  data_frame <- data_frame |>
    tidyr::unite(
      "pathName",                       # New column name
      tidyselect::all_of(colsToConcat), # Columns to unite
      sep = " | ",                      # Separator
      na.rm = TRUE                      # Skip NA values
    )
  
  return(data_frame)
}

