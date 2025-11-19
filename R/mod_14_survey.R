
#' 14_survey UI Function
#'
#' @description A shiny Module.
#' @param id,input,output,session Internal parameters for {shiny}.
#' @noRd
#' @importFrom shiny NS tagList
mod_14_survey_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Text full width of modal
    tags$style(HTML("
      .modal-body .form-group,
      .modal-body .shiny-input-container {
        width: 100% !important;
      }
    ")),
    # JS logic: trigger modal only if localStorage flag not set
    tags$script(HTML(sprintf("
      $(document).on('shiny:connected', function() {
        var ns_input = '%s';
        if (!localStorage.getItem('survey_done')) {
          // delay to let UI render
          setTimeout(function() {
            Shiny.setInputValue(ns_input, Math.random());
          }, 500);
        }
      });
      Shiny.addCustomMessageHandler('survey_complete', function(message) {
        localStorage.setItem('survey_done', true);
      });
    ", ns("show_survey"))))
  )
}

#' 14_survey Server Functions
#' @noRd
mod_14_survey_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # File paths
    survey_file <- "survey_responses.csv"
    survey_log_file <- "survey_log.csv"

    # Ensure files exist
    if (!file.exists(survey_file)) {
      write.csv(
        data.frame(
          timestamp = character(),
          ip_hash = character(),
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
    if (!file.exists(survey_log_file)) {
      write.csv(
        data.frame(ip_hash = character(), email = character(), stringsAsFactors = FALSE),
        survey_log_file, row.names = FALSE
      )
    }

    # Helper: get hashed client IP
    get_hashed_ip <- reactive({
      ip <- session$clientData$remoteAddress
      if (is.null(ip) || ip == "") ip <- "unknown"
      digest::digest(ip)
    })

    # Show modal once
    observeEvent(input$show_survey, {
      ip_hash <- get_hashed_ip()

      # Read log safely
      log_data <- tryCatch(read.csv(survey_log_file, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(log_data) && ip_hash %in% log_data$ip_hash) return()

      # Show modal
      showModal(modalDialog(
        title = "Thank you for using iDEP or ShinyGO. We'd love your feedback. Just 6 quick questions (~30 seconds).",
        size = "m",
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("submit_survey"), "Submit", class = "btn-primary")
        ),

        # Q1: Organization (single)
        div(style = "width:100%;", # ensure q&a displays full width of modal
          radioButtons(
            ns("q1_org"),
            label = tagList(
              "Q1. What best describes your organization? (Select one.)",
              tags$span("*", style = "color:red;", title = "Required")
            ),
            choices = c(
              "University or academic research institute",
              "Hospital or medical center",
              "Pharma or biotech company",
              "Government agency",
              "Nonprofit or foundation",
              "Core facility / service lab",
              "Other (please specify)"
            ),
            selected = character(0)
          ),
          conditionalPanel(
            sprintf("input['%s'] == 'Other (please specify)'", ns("q1_org")),
            textInput(ns("q1_other"), "Please specify:", "")
          ),

          # Q2: Role (multi-select)
          checkboxGroupInput(
            ns("q2_role"),
            label = tagList(
              "Q2. What best describes your role and involvement with data analysis/tools? (Select all that apply.)",
              tags$span("*", style = "color:red;", title = "Required")
            ),
            choices = c(
              "Principal investigator",
              "Bioinformatician",
              "Bench / experimental scientist",
              "Core facility / service lab staff",
              "Student / trainee",
              "Instructor / teacher",
              "Other (please specify)"
            )
          ),
          conditionalPanel(
            sprintf("input['%s'].includes('Other (please specify)')", ns("q2_role")),
            textInput(ns("q2_other"), "Please specify:", "")
          ),

          # Q3: How you use iDEP/ShinyGO (single)
          radioButtons(
            ns("q3_use"),
            label = tagList(
              "Q3. How do you primarily use iDEP or ShinyGO in your work? (Select one.)",
              tags$span("*", style = "color:red;", title = "Required")
            ),
            choices = c(
              "Exploring new datasets and generating hypotheses",
              "Routine analysis as part of a standard workflow",
              "Quick QC or sanity checks on results",
              "Teaching or training (e.g., courses, workshops)",
              "Method development or benchmarking",
              "Other (please specify)"
            ),
            selected = character(0)
          ),
          conditionalPanel(
            sprintf("input['%s'] == 'Other (please specify)'", ns("q3_use")),
            textInput(ns("q3_other"), "Please specify:", "")
          ),

          # Q4: Frequency (single)
          radioButtons(
            ns("q4_freq"),
            label = tagList(
              "Q4. How often do you use iDEP or ShinyGO? (Select one.)",
              tags$span("*", style = "color:red;", title = "Required")
            ),
            choices = c(
              "This is my first time!",
              "Occasionally (About 1–5 times per year)",
              "Monthly (About 1–5 times per month)",
              "Weekly (About 1–5 times per week)",
              "Daily"
            ),
            selected = character(0)
          ),

          # Q5: Improvements (multi-select)
          checkboxGroupInput(
            ns("q5_improve"),
            label = tagList(
              "Q5. If we could improve iDEP / ShinyGO for you, what would be most valuable? (Select as many as you want)",
              tags$span("*", style = "color:red;", title = "Required")
            ),
            choices = c(
              "Support for more public datasets",
              "User community (webinars, community forums, etc.)",
              "Features that make it easier to cite, document, or reproduce analyses",
              "Local deployment options",
              "Consultation or custom analysis help",
              "New analysis modules (e.g., single-cell, multi-omics, time series)",
              "More polished reports and export formats",
              "Nothing (good as is)",
              "Other (please specify)"
            )
          ),
          conditionalPanel(
            sprintf("input['%s'].includes('Other (please specify)')", ns("q5_improve")),
            textInput(ns("q5_other"), "Please specify:", "")
          ),

          # Q6: Optional notes
          textAreaInput(
            ns("q6_notes"),
            "Q6. Anything you want us to know? (Optional)",
            "",
            rows = 3
          ),

          # Q7: Optional email
          textAreaInput(
            ns("q7_email"),
            "Q7. Keep up with new features — leave your email (Optional)",
            "",
            rows = 1
          )
        )
      ))
    }, ignoreNULL = TRUE, once = TRUE) # <-- ensures it triggers only once per session

    # Submission
    observeEvent(input$submit_survey, {
      required <- c(input$survey_city, input$survey_role,
                    input$survey_like, input$survey_improve, input$survey_features)
      if (any(sapply(required, function(x) x == ""))) {
        showNotification("Please fill in all required fields before submitting.", type = "error")
        return()
      }

      ip_hash <- get_hashed_ip()
      ts <- Sys.time()

      # Save full survey
      new_entry <- data.frame(
        timestamp = as.character(Sys.time()),
        ip_hash = ip_hash,
        org = input$q1_org,
        org_other = ifelse(input$q1_org == "Other (please specify)", input$q1_other, ""),
        role = paste(input$q2_role, collapse = "; "),
        role_other = ifelse("Other (please specify)" %in% input$q2_role, input$q2_other, ""),
        use = input$q3_use,
        use_other = ifelse(input$q3_use == "Other (please specify)", input$q3_other, ""),
        frequency = input$q4_freq,
        improve = paste(input$q5_improve, collapse = "; "),
        improve_other = ifelse("Other (please specify)" %in% input$q5_improve, input$q5_other, ""),
        notes = input$q6_notes,
        stringsAsFactors = FALSE
      )
      write.table(new_entry, survey_file, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)

      # Save minimal log
      new_log <- data.frame(
        ip_hash = ip_hash,
        email = ifelse(input$q7_email == "", NA, input$q7_email),
        stringsAsFactors = FALSE
      )
      write.table(new_log, survey_log_file, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)

      removeModal()
      showNotification("Thank you for your feedback!", type = "message")
      session$sendCustomMessage("survey_complete", TRUE)
    })
  })
}


## To be copied in the UI
# mod_14_survey_ui("14_survey_1")

## To be copied in the server
# mod_14_survey_server("14_survey_1")