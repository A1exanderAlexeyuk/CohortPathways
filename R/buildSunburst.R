#' Create a Sunburst Plot for Treatment Pathways
#'
#' This function creates an interactive sunburst plot to visualize treatment pathways
#' from CohortPathway analysis results using the plotly package. The plot shows the
#' hierarchical structure of treatment sequences, with each ring representing a step
#' in the pathway.
#'
#' @param cpResults A list containing the results from CohortPathway analysis.
#'        Must include 'pathwaysAnalysisPathsData' and 'isCombo' data frames.
#' @param generationSet A data frame containing information about cohorts
#'        that will be used to generate descriptive names for the events in the
#'        diagram and the plot titles.
#' @param maxPaths Integer specifying the maximum number of steps (depth) to include in the plot.
#' @param minCount Integer specifying the minimum count value for a path to be included.
#'
#' @return A named list of plotly objects, with each element corresponding to a target cohort.
#' @export
#'
#' @examples
#' \dontrun{
#' library(CohortPathways)
#' # Assuming cpResults and generationSet are loaded
#' sunburstPlots <- createPathwaySunburst(
#'   cpResults = cpResults,
#'   generationSet = generationSet,
#'   maxPaths = 5,
#'   minCount = 10
#' )
#' # To display the plot for the first target cohort
#' sunburstPlots[[1]]
#' }
createPathwaySunburst <- function(cpResults,
                                  generationSet,
                                  maxPaths = 5,
                                  minCount = 5) {
  # Check for required packages
  rlang::check_installed(c("plotly", "RColorBrewer"))
  
  # Input validation
  checkmate::assertList(cpResults, min.len = 1, types = "data.frame")
  checkmate::assertDataFrame(generationSet, min.rows = 1)
  checkmate::assertInt(maxPaths, lower = 1)
  checkmate::assertInt(minCount, lower = 0)
  
  # Extract and prepare data
  pathwaysAnalysisPathsDatas <- purrr::pluck(cpResults, "pathwaysAnalysisPathsData") |>
    dplyr::filter(.data$countValue >= minCount) |>
    dplyr::group_by(.data$targetCohortId) |>
    dplyr::group_split()
  
  isCombo <- purrr::pluck(cpResults, "isCombo")
  eventNames <- .prepareEventNames(generationSet, cpResults)
  codeToName <- rlang::set_names(eventNames$combination, as.character(eventNames$code))
  comboMap <- rlang::set_names(isCombo$isCombo, as.character(isCombo$comboId))
  
  # Helper function to get descriptive name for a code
  getName <- function(code) {
    name <- codeToName[as.character(code)]
    isComboFlag <- comboMap[as.character(code)]
    if (!is.na(isComboFlag) && isComboFlag == 1) {
      return(paste0(name, " (combo)"))
    }
    return(name)
  }
  
  # Get names for target cohorts to use in plot titles
  targetNames <- purrr::map_chr(pathwaysAnalysisPathsDatas, function(df) {
    tId <- unique(df$targetCohortId)
    generationSet |>
      dplyr::filter(.data$cohortId %in% tId) |>
      dplyr::pull(.data$cohortName) |>
      unique()
  })
  
  # Generate a plot for each target cohort
  plots <- lapply(seq_along(pathwaysAnalysisPathsDatas), function(i) {
    pathwaysData <- pathwaysAnalysisPathsDatas[[i]]
    
    if (nrow(pathwaysData) == 0) {
      warning(paste("No pathways with count >=", minCount, "for target cohort", targetNames[i]))
      return(NULL)
    }
    
    # Initialize data frame for sunburst plot
    sunburstData <- dplyr::tibble(
      ids = "Start",
      labels = "Start",
      parents = "",
      values = sum(pathwaysData$countValue),
      color = "#1f77b4" # Blue for Start
    )
    
    # Process each pathway
    for (r in 1:nrow(pathwaysData)) {
      path <- c("Start")
      row <- pathwaysData[r, ]
      
      for (j in 1:maxPaths) {
        stepCol <- paste0("step", j)
        if (!stepCol %in% names(row) || is.na(row[[stepCol]])) {
          break
        }
        
        code <- as.character(row[[stepCol]])
        name <- getName(code)
        
        # Create a unique ID for each node in the hierarchy
        parentId <- paste(path, collapse = "/")
        path <- c(path, name)
        currentId <- paste(path, collapse = "/")
        
        # Add new node if it doesn't exist
        if (!currentId %in% sunburstData$ids) {
          isComboFlag <- comboMap[code]
          nodeColor <- if (!is.na(isComboFlag) && isComboFlag == 1) "#ff7f0e" else "#aec7e8"
          
          sunburstData <- sunburstData |>
            dplyr::add_row(
              ids = currentId,
              labels = name,
              parents = parentId,
              values = 0,
              color = nodeColor
            )
        }
      }
    }
    
    # Aggregate values for each path segment
    for (r in 1:nrow(pathwaysData)) {
      path_val <- pathwaysData$countValue[r]
      path_nodes <- c("Start")
      for (j in 1:maxPaths) {
        stepCol <- paste0("step", j)
        if (!stepCol %in% names(pathwaysData) || is.na(pathwaysData[[stepCol]][r])) {
          break
        }
        code <- as.character(pathwaysData[[stepCol]][r])
        name <- getName(code)
        path_nodes <- c(path_nodes, name)
        currentId <- paste(path_nodes, collapse="/")
        
        idx <- which(sunburstData$ids == currentId)
        sunburstData$values[idx] <- sunburstData$values[idx] + path_val
      }
    }
    
    
    # Create the plot
    fig <- plotly::plot_ly(
      data = sunburstData,
      ids = ~ids,
      labels = ~labels,
      parents = ~parents,
      values = ~values,
      type = 'sunburst',
      branchvalues = 'total',
      hovertemplate = '<b>%{label}</b><br>Patients: %{value}<br><extra></extra>',
      marker = list(colors = ~color)
    )
    
    fig <- fig |>
      plotly::layout(
        title = paste("Treatment Pathways for:", targetNames[i]),
        margin = list(l = 10, r = 10, b = 10, t = 50)
      )
    
    return(fig)
  })
  
  # Remove NULL plots and return named list
  plots <- purrr::compact(plots)
  return(plots |> rlang::set_names(targetNames[seq_along(plots)]))
}

#' Create a Sankey Plot for Treatment Pathways
#'
#' This function creates an interactive Sankey diagram to visualize the flow of
#' patients through different treatment pathways.
#'
#' @param cpResults A list containing the results from CohortPathway analysis.
#' @param generationSet A data frame containing information about event cohorts.
#' @param maxPaths The maximum number of steps (depth) to include in the plot.
#' @param minCount The minimum count of a pathway to be included in the visualization.
#' @return A Sankey network plot created with the networkD3 package.
#' @export
#'
#' @examples
#' \dontrun{
#' library(CohortPathways)
#' # Assuming cpResults and generationSet are loaded
#' sankeyPlot <- createPathwaySankey(
#'   cpResults = cpResults,
#'   generationSet = generationSet,
#'   maxPaths = 5,
#'   minCount = 10
#' )
#' sankeyPlot
#' }
createPathwaySankey <- function(cpResults,
                                generationSet,
                                maxPaths = 5,
                                minCount = 10) {
  # Check for required packages
  rlang::check_installed(c("networkD3", "htmlwidgets", "RColorBrewer"))
  
  # Input validation
  checkmate::assertList(cpResults, min.len = 1)
  checkmate::assertDataFrame(generationSet, min.rows = 1)
  checkmate::assertInt(maxPaths, lower = 1)
  checkmate::assertInt(minCount, lower = 0)
  
  # Prepare data
  pathwayData <- purrr::pluck(cpResults, "pathwaysAnalysisPathsData") |>
    dplyr::filter(.data$countValue >= minCount)
  
  if (nrow(pathwayData) == 0) {
    stop("No pathways meet the minimum count threshold.")
  }
  
  eventNames <- .prepareEventNames(generationSet, cpResults)
  codeToName <- rlang::set_names(eventNames$combination, as.character(eventNames$code))
  
  # Build links data frame
  links <- dplyr::tibble()
  stepCols <- paste0("step", 1:maxPaths)
  
  # Add links from Start to Step 1
  step1_links <- pathwayData |>
    dplyr::filter(!is.na(.data$step1)) |>
    dplyr::group_by(.data$step1) |>
    dplyr::summarise(value = sum(.data$countValue)) |>
    dplyr::transmute(
      source_name = "Start",
      target_name = paste(codeToName[as.character(.data$step1)], "1", sep = "_"),
      value = .data$value
    )
  links <- dplyr::bind_rows(links, step1_links)
  
  # Add links between consecutive steps
  for (i in 1:(maxPaths - 1)) {
    sourceCol <- stepCols[i]
    targetCol <- stepCols[i + 1]
    
    if (!sourceCol %in% names(pathwayData) || !targetCol %in% names(pathwayData)) {
      break
    }
    
    step_links <- pathwayData |>
      dplyr::filter(!is.na(.data[[sourceCol]]) & !is.na(.data[[targetCol]])) |>
      dplyr::group_by(.data[[sourceCol]], .data[[targetCol]]) |>
      dplyr::summarise(value = sum(.data$countValue), .groups = "drop") |>
      dplyr::transmute(
        source_name = paste(codeToName[as.character(.data[[sourceCol]])], i, sep = "_"),
        target_name = paste(codeToName[as.character(.data[[targetCol]])], i + 1, sep = "_"),
        value = .data$value
      )
    links <- dplyr::bind_rows(links, step_links)
  }
  
  # Create nodes data frame
  node_names <- unique(c(links$source_name, links$target_name))
  nodes <- dplyr::tibble(
    name = node_names,
    group = gsub("_\\d+$", "", .data$name) # Group by event name for coloring
  ) |>
    dplyr::arrange(.data$name)
  
  # Create a color scale
  unique_groups <- unique(nodes$group)
  colors <- RColorBrewer::brewer.pal(length(unique_groups), "Set3")
  color_scale <- paste0('d3.scaleOrdinal().domain(["', paste(unique_groups, collapse = '","'), '"]).range(["', paste(colors, collapse = '","'), '"])')
  
  
  # Convert node names to 0-based IDs for networkD3
  links$source <- match(links$source_name, nodes$name) - 1
  links$target <- match(links$target_name, nodes$name) - 1
  nodes$name <- gsub("_\\d+$", "", nodes$name) # Clean up names for display
  
  # Create the Sankey plot
  sankey <- networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "group",
    fontSize = 12,
    nodeWidth = 30,
    sinksRight = FALSE,
    colourScale = color_scale
  )
  
  return(sankey)
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
