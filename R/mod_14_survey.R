
#' 14_survey UI Function
#'
#' @description A shiny Module.
#' @param id,input,output,session Internal parameters for {shiny}.
#' @noRd
#' @importFrom shiny NS tagList
mod_14_survey_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # Allow full-width layout inside modal
    tags$style(HTML("
      .modal-body .form-group,
      .modal-body .shiny-input-container {
        width: 100% !important;
      }
    ")),

    # LocalStorage + Shiny bridge
    tags$head(
      tags$script(src = "www/survey_storage.js"),
        tags$script(HTML(paste0(
          "var surveyNamespace = '", ns(""), "';",
          "var surveyInputId = '", ns("survey_done_storage"), "';"
    ))))
  )
}

#' 14_survey Server Functions
#' @noRd
mod_14_survey_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Configuration (GDPR-related)
    # Retention: delete entries older than this many months
    retention_months <- 12

    # Create data directory in the app's working directory
    data_dir <- file.path(getwd(), "data")
    if (!dir.exists(data_dir)) {
      dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    }
    survey_file <- file.path(data_dir, "survey_responses.csv")

    if (!file.exists(survey_file)) {
      write.csv(
        data.frame(
          timestamp = character(),
          org = character(),
          org_other = character(),
          role = character(),
          role_other = character(),
          use = character(),
          use_other = character(),
          frequency = character(),
          improve = character(),
          improve_other = character(),
          notes = character(),
          stringsAsFactors = FALSE
        ),
        survey_file, row.names = FALSE
      )
    }

    # Helper: retention purge
    purge_old_entries <- function() {
      # Purge survey responses older than retention_months
      resp <- tryCatch(read.csv(survey_file, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(resp) && nrow(resp) > 0) {
        # robust timestamp parse
        ts_parsed <- suppressWarnings(as.POSIXct(resp$timestamp, tz = "UTC"))
        cutoff <- Sys.time() - (retention_months * 30.4375 * 24 * 3600) # average month length
        keep_idx <- !is.na(ts_parsed) & ts_parsed >= cutoff
        resp <- resp[keep_idx, , drop = FALSE]
        write.csv(resp, survey_file, row.names = FALSE)
      }
    }
    purge_old_entries()

    # Show modal once
    observe({
      if (!is.null(input$survey_done_storage)) {
        # Show modal only if survey is NOT done (FALSE)
        if (!isTRUE(input$survey_done_storage)) {

          # Notice text
          notice <- tags$div(
            style = "background:#f8f9fa;border:1px solid #e3e6ea;border-radius:6px;
                    padding:12px;padding-bottom:7px;margin-bottom:12px;font-size:1.00em;",
            tags$p("Thank you for using iDEP! Please complete a quick 6 question survey (~30 seconds).")
          )

          skip_button <- actionButton(ns("decline_survey"), "Skip for now", class = "btn-secondary")

          showModal(modalDialog(
            title = tagList(
              div(
                style = "display:flex;justify-content:space-between;align-items:center;gap:10px;width:100%;",
                tags$span("We'd love your feedback!"),
                skip_button
              )
            ),
            size = "m",
            easyClose = FALSE,
            footer = tagList(
              actionButton(ns("submit_survey"), "Submit", class = "btn-primary")
            ),

            # Privacy notice (top)
            notice,

            # Q1: Organization (single)
            div(style = "width:100%;", # ensure q&a displays full width of modal
              selectInput(
                inputId = ns("q1_org"),
                label = tagList(
                  "1. What best describes your organization?"
                ),
                choices = c(
                  "..." = "",
                  "University / academic institute",
                  "Hospital or affiliated research center",
                  "Pharma / biotech company",
                  "Government agency",
                  "Nonprofit",
                  "Core facility / service provider",
                  "Other (please specify)"
                ),
                selected = "",
                selectize = FALSE
              ),
              conditionalPanel(
                sprintf("input['%s'] == 'Other (please specify)'", ns("q1_org")),
                textInput(ns("q1_other"), "Please specify:", "")
              ),

              # Q2: Role (multi-select)
                selectizeInput(
                  inputId = ns("q2_role"),
                  label = tagList(
                    "2. What best describes your roles? (Select all that apply.)"
                  ),
                  choices = c(
                    "PI / group leader",
                    "Bioinformatician / data scientist",
                    "Bench scientist",
                    "Core facility / service staff",
                    "Educator",
                    "Student / trainee",
                    "Other (please specify)"
                  ),
                  selected = character(0),
                  multiple = TRUE,
                  options = list(placeholder = "Select all that apply")
                ),
              uiOutput(ns("q2_other_ui")),

              # Q3: How you use iDEP (single)
              selectInput(
                inputId = ns("q3_use"),
                label = tagList(
                  "3. How do you primarily use iDEP in your work?"
                ),
                choices = c(
                  "..." = "",
                  "Explore datasets / generate hypotheses",
                  "Routine workflow analysis",
                  "QC or sanity checks",
                  "Teaching or training",
                  "Other (please specify)"
                ),
                selected = "",
                selectize = FALSE
              ),
              conditionalPanel(
                sprintf("input['%s'] == 'Other (please specify)'", ns("q3_use")),
                textInput(ns("q3_other"), "Please specify:", "")
              ),

              # Q4: Frequency (single)
              selectInput(
                inputId = ns("q4_freq"),
                label = tagList(
                  "4. How often do you use iDEP?"
                ),
                choices = c(
                  "..." = "",
                  "First time user",
                  "Occasional (1–5× per year)",
                  "Monthly (1–5× per month)",
                  "Weekly (1–5× per week)",
                  "Daily (5+ times per week)"
                ),
                selected = "",
                selectize = FALSE
              ),

              # Q5: Improvements (multi-select)
                selectizeInput(
                  inputId = ns("q5_improve"),
                  label = tagList(
                    "5. To improve iDEP, what would be most valuable? (Select as many as you want)"
                  ),
                  choices = c(
                    "More public datasets",
                    "Community events / webinars",
                    "Reproducibility features",
                    "New analysis modules (e.g., single-cell, multi-omics, time series)",
                    "Add AI assistant features to guide me through analysis",
                    "Nothing (good as is)",
                    "Other (please specify)"
                  ),
                  selected = character(0),
                  multiple = TRUE,
                  options = list(placeholder = "Select as many as apply")
                ),
              uiOutput(ns("q5_other_ui")),

              # Q6: Optional notes
              textAreaInput(
                ns("q6_notes"),
                "6. Any other suggestions or feedback (Optional)",
                "",
                rows = 3
              )
            )
          ))
        }
      }
    })

    # Conditional UI for Q2 "Other" text input
    output$q2_other_ui <- renderUI({
      if ("Other (please specify)" %in% input$q2_role) {
        textInput(ns("q2_other"), "Please specify:", "")
      }
    })

    # Conditional UI for Q5 "Other" text input
    output$q5_other_ui <- renderUI({
      if ("Other (please specify)" %in% input$q5_improve) {
        textInput(ns("q5_other"), "Please specify:", "")
      }
    })

    # Submission
    observeEvent(input$submit_survey, {

      # Ensure if 'other' is selected, the short response is answered
      if (identical(input$q1_org, "Other (please specify)") &&
          (is.null(input$q1_other) || input$q1_other == "")) {
        showNotification("Please specify your organization in Q1.", type = "error")
        return()
      }

      if ("Other (please specify)" %in% input$q2_role &&
          (is.null(input$q2_other) || input$q2_other == "")) {
        showNotification("Please specify your role in Q2.", type = "error")
        return()
      }

      if (identical(input$q3_use, "Other (please specify)") &&
          (is.null(input$q3_other) || input$q3_other == "")) {
        showNotification("Please specify your answer in Q3.", type = "error")
        return()
      }

      if ("Other (please specify)" %in% input$q5_improve &&
          (is.null(input$q5_other) || input$q5_other == "")) {
        showNotification("Please specify your answer in Q5.", type = "error")
        return()
      }

      # Store full survey responses
      # Helper to safely get input values (handle NULL and empty strings)
      safe_input <- function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        x <- trimws(as.character(x))
        if (x == "") return(NA_character_)
        return(x)
      }
      
      # Helper for multi-select inputs
      safe_collapse <- function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        return(paste(x, collapse = "; "))
      }
      
      new_entry <- data.frame(
        timestamp = as.character(Sys.time()),
        org = safe_input(input$q1_org),
        org_other = if (identical(input$q1_org, "Other (please specify)")) safe_input(input$q1_other) else NA_character_,
        role = safe_collapse(input$q2_role),
        role_other = if ("Other (please specify)" %in% input$q2_role) safe_input(input$q2_other) else NA_character_,
        use = safe_input(input$q3_use),
        use_other = if (identical(input$q3_use, "Other (please specify)")) safe_input(input$q3_other) else NA_character_,
        frequency = safe_input(input$q4_freq),
        improve = safe_collapse(input$q5_improve),
        improve_other = if ("Other (please specify)" %in% input$q5_improve) safe_input(input$q5_other) else NA_character_,
        notes = safe_input(input$q6_notes),
        stringsAsFactors = FALSE
      )

      # Save responses with na parameter to ensure consistent NA writing
      data.table::fwrite(new_entry, survey_file, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE, na = "NA")

      purge_old_entries()

      removeModal()
      showNotification("Thank you for your feedback!", type = "message")

      # mark localStorage
      session$sendCustomMessage("markSurveyComplete", list())
    })

    # Opt-out
    observeEvent(input$decline_survey, {
      removeModal()
      session$sendCustomMessage("markSurveyComplete", list())
    })
  })
}


## To be copied in the UI
# mod_14_survey_ui("14_survey_1")

## To be copied in the server
# mod_14_survey_server("14_survey_1")
