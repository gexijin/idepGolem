#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_01_load_data_ui <- function(id){
  ns <- shiny::NS(id)
  tagList(
    selectizeInput(
      inputId = ns("select_org"),
      label = strong("1. Select or search for your species."),
      choices = NULL,
      multiple = TRUE,
      options = list(
        maxItems = 1,
        placeholder = "Best matching species",
        onInitialize = I("function() { this.setValue(\"\"); }")))
  )
}
    
#' 01_load_data Server Functions
#'
#' @noRd 
mod_01_load_data_server <- function(id, idep_data) {
  species_choice <- setNames(
    object = as.list(idep_data$org_info$id), # values
    nm = idep_data$org_info$name2 # name of variables
  )
  species_choice <- append(
    setNames(
      object = "NEW", # values
      nm = "**NEW SPECIES**" # name of variables $
    ),
    species_choice
  )
  species_choice <- append(
    setNames(
      object = "BestMatch",
      nm = "Best matching species"
    ),
    species_choice
  )

  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = species_choice,
        server = TRUE
      )
    }
  )
}
    
## To be copied in the UI
# mod_01_load_data_ui("01_load_data_ui_1")
    
## To be copied in the server
# mod_01_load_data_server("01_load_data_ui_1")
