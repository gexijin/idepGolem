
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

    # Configuration (GDPR-related)
    # Retention: delete entries older than this many months
    retention_months <- 12

    # Privacy policy URL and contact for deletion/withdrawal requests
    privacy_policy_url <- "https://www.orditus.com/privacy-policy/"  # update with iDEP specific privacy policy
    privacy_contact <- "info@orditus.com"

    # Salt/secret for hashing (set in config.yml file for security)
    config <- yaml::read_yaml("config.yml")
    ip_salt <- config$survey_ip_salt
    if (is.na(ip_salt) || ip_salt == "") {
      # Fallback: generate a volatile in-memory salt (not ideal for cross-session consistency).
      # For production, REQUIRE a persistent secret via env var or config file.
      ip_salt <- digest::digest(paste0("fallback-", runif(1)), algo = "sha256")
      warning("SURVEY_IP_SALT not set. Using a volatile fallback salt. Set SURVEY_IP_SALT for consistent deduplication.")
    }

    # File paths
    survey_file <- "survey_responses.csv"  # full responses + salted IP hash + consent flag
    survey_log_file <- "survey_log.csv"    # emails only (if provided)

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
        data.frame(email = character(), stringsAsFactors = FALSE), # only email here, no ip
        survey_log_file, row.names = FALSE
      )
    }

    # Helper: salted/secret-keyed hashed client IP (pseudonymization)
    get_hashed_ip <- reactive({
      ip <- session$clientData$remoteAddress
      if (is.null(ip) || ip == "") ip <- "unknown"
      # Use HMAC with sha256 to bind hash to server-held secret (resistant to rainbow tables)
      digest::hmac(key = ip_salt, object = ip, algo = "sha256")
    })

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

    # Show modal once (client + server dedup)
    observeEvent(input$show_survey, {
      ip_hash <- get_hashed_ip()

      # Server-side dedup: if this salted IP already submitted, don't show modal
      resp <- tryCatch(read.csv(survey_file, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(resp) && nrow(resp) > 0 && ip_hash %in% resp$ip_hash_salted) return()

      # Privacy notice text (transparent, explicit, concise)
      privacy_notice <- tags$div(
        style = "background:#f8f9fa;border:1px solid #e3e6ea;border-radius:6px;padding:12px;margin-bottom:12px;font-size:0.95em;",
        tags$p(
          tags$b("Privacy notice: "),
          "Your responses will be used solely to improve this application. ",
          "Email is optional and responses are retained for up to ", retention_months, " months and then deleted."
        ),
        tags$p(
          "By submitting, you agree to the storage and processing of your data for feedback purposes. ",
          "You can withdraw consent and request access or deletion at any time by contacting ",
          tags$a(href = paste0("mailto:", privacy_contact), privacy_contact), ". ",
          "See our ",
          tags$a(href = privacy_policy_url, target = "_blank", "privacy policy"),
          " for details on processing and your rights."
        )
      )

      # Consent checkbox (must be checked)
      consent_ui <- checkboxInput(
        inputId = ns("survey_consent"),
        label = tags$span(
          HTML("I consent to my data being stored and processed for feedback purposes."),
          tags$span("*", style = "color:red;margin-left:6px;", title = "Required")
        ),
        value = FALSE
      )

      showModal(modalDialog(
        title = "Thank you for using iDEP or ShinyGO. We'd love your feedback. Just 6 quick questions (~30 seconds).",
        size = "m",
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("submit_survey"), "Submit", class = "btn-primary"),
          actionButton(ns("decline_survey"), "I do not want to participate", class = "btn-secondary")
        ),

        # Privacy notice (top)
        privacy_notice,

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
        ),

        # Consent (bottom, required)
        tags$hr(),
        consent_ui
      ))
    }, ignoreNULL = TRUE, once = TRUE) # <-- ensures it triggers only once per session

    # Submission
    observeEvent(input$submit_survey, {

      # Ensure required questions are answered
      if (is.null(input$q1_org) || input$q1_org == "") {
        showNotification("Please answer Question 1.", type = "error")
        return()
      }
      if (identical(input$q1_org, "Other (please specify)") &&
          (is.null(input$q1_other) || input$q1_other == "")) {
        showNotification("Please specify your organization in Q1.", type = "error")
        return()
      }

      if (is.null(input$q2_role) || length(input$q2_role) == 0) {
        showNotification("Please answer Question 2.", type = "error")
        return()
      }
      if ("Other (please specify)" %in% input$q2_role &&
          (is.null(input$q2_other) || input$q2_other == "")) {
        showNotification("Please specify your role in Q2.", type = "error")
        return()
      }

      if (is.null(input$q3_use) || input$q3_use == "") {
        showNotification("Please answer Question 3.", type = "error")
        return()
      }
      if (identical(input$q3_use, "Other (please specify)") &&
          (is.null(input$q3_other) || input$q3_other == "")) {
        showNotification("Please specify your answer in Q3.", type = "error")
        return()
      }

      if (is.null(input$q4_freq) || input$q4_freq == "") {
        showNotification("Please answer Question 4.", type = "error")
        return()
      }

      if (is.null(input$q5_improve) || length(input$q5_improve) == 0) {
        showNotification("Please answer Question 5.", type = "error")
        return()
      }
      if ("Other (please specify)" %in% input$q5_improve &&
          (is.null(input$q5_other) || input$q5_other == "")) {
        showNotification("Please specify your answer in Q5.", type = "error")
        return()
      }

      # Consent required
      if (!isTRUE(input$survey_consent)) {
        showNotification("Consent is required to submit this survey.", type = "error")
        return()
      }

      ip_hash <- get_hashed_ip()

      # Store full survey responses
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

      # Save responses
      data.table::fwrite(new_entry, survey_file, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)

      # Save email log if provided
      if (!is.null(input$q7_email) && nzchar(input$q7_email)) {
        data.table::fwrite(data.table::data.table(email = input$q7_email),
                           survey_log_file, sep = ",", append = TRUE,
                           col.names = FALSE, row.names = FALSE)
      }

      removeModal()
      showNotification("Thank you for your feedback!", type = "message")
      session$sendCustomMessage("survey_complete", TRUE)
      purge_old_entries()  # retention purge after retention_months
    })

    # Opt-out
    observeEvent(input$decline_survey, {
      removeModal()
      session$sendCustomMessage("survey_complete", TRUE)
    })
  })
}


## To be copied in the UI
# mod_14_survey_ui("14_survey_1")

## To be copied in the server
# mod_14_survey_server("14_survey_1")