#' 10_doc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_10_doc_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "About",
    fluidPage(
      h3(
        "If you find iDEP helpful, please ",
        a(
          "send us a brief email.",
          href = "mailto:gelabinfo@gmail.com?Subject=iDEP support letter"
        )
      ),
      p(" If you state your general research area and how iDEP
      makes you more productive, we can use it as a support letter when we
        apply for the next round
        of funding.
        Hundreds of strong, enthusiastic letters sent to us in 2019
        were essential when we applied for the current
        grant from NIH/NHGRI (R01HG010805),
        which expires in 20 months. Your letters will help sustain
        and improve this service."),
      p(
        "iDEP is developed and maintained by a small team at ",
        a(
          "South Dakota State University (SDSU). ",
          href = "https://www.sdstate.edu/"
        ),
        "Our team consists
      of Xijin Ge (PI), Jianli Qi (research staff), and
      two talented graduate students (Emma Spors and Ben Derenge).
      None of us are trained as software engineers. But
      we share the passion about  developing an
      user-friendly tool for all biologists,
      especially those who do not have access to bioinformaticians."
      ),
      p("Past contributors include Eun Wo Son, Runan Yao,
      Gavin Doering, Roberto Villegas-Diaz, and Eric Tulowetzke.
      Eric still helps us fix bugs after leaving the lab.
      Much of the new version of iDEP is rewritten by Gavin Doering.
      The iDEP logo was designed by Emma Spors.
      Technical support is kindly provided by the Office of Information
      Technology (OIT) at SDSU. Mirror site is enabled by a JetStream2
      allocation award (BIO210175), which is supported by NSF.
      "),
      h3("Please cite the BMC Genomics paper!"),
      p(
        "If you use iDEP, please cite: Ge, Son & Yao, iDEP: an integrated web application for differential
        expression and pathway analysis of RNA-Seq data, ",
        a("BMC Bioinformatics 19:1-24, 2018.",
          href = "https://doi.org/10.1186/s12859-018-2486-6",
          target = "_blank"
        ),
        "Merely mentioning iDEP with an URL is insufficient. It is difficult to track.",
        "Consider citing other tools that form the foundation of iDEP, such as ",
        a("ENSEMBL, ",
          href = "https:/doi.org/10.1093/nar/gkab1049",
          target = "_blank"
        ),
        a(" STRING-db,",
          href = "https://doi.org/10.1093/nar/gky1131",
          target = "_blank"
        ),
        a(" DESeq2, ",
          href = "https://doi.org/10.1186/s13059-014-0550-8",
          target = "_blank"
        ),
        a(" limma",
          href = "https://doi.org/10.1093/nar/gkv007",
          target = "_blank"
        ),
        " and many others.",
        " If you use the KEGG diagram, please also cite ",
        a("pathview, ",
          href = "https://doi.org/10.1093/bioinformatics/btt285",
          target = "_blank"
        ),
        "and ",
        a("KEGG.",
          href = "https://doi.org/10.1093/nar/gkaa970",
          target = "_blank"
        )
      ),
      p(
        "According to Google Scholar, more than",
        a("450 papers cited iDEP,",
          href = "https://scholar.google.com/scholar?oi=bibs&hl=en&cites=6502699637682046008,17999801138713500070,11001860275874506471",
          target = "_blank"
        ),
        " as of Sept 18, 2022.",
        "But our website has been accessed over 316,000 times by 83,000 users,
        spending 17.5 minutes each time. For every 1000 users, only 6 cited
        the iDEP paper, which is disappointingly low. Consider citing iDEP even
        if you used it just for preliminary analysis. Otherwise, this tool might
        just vanish."
      ),
      h3("Source code and database"),
      p(
        "Source code is available on ",
        a(
          "GitHub,",
          href = "https://github.com/espors/idepGolem/"
        ),
        " which also includes instructions to install
        iDEP on your local machine using our ",
        a(
          "database.",
          href = "http://bioinformatics.sdstate.edu/data/"
        )
      ),
      h3("Contact us"),
      p(
        "Please email",
        a(
          " Jenny",
          href = "mailto:gelabinfo@gmail.com?Subject=iDEP"
        ),
        " (recommended) or ",
        a("Dr. Ge",
          href = "mailto:xijin.ge@sdstate.edu?Subject=iDEP"
        ),
        " (unpredictable). Follow us on ",
        a("Twitter",
          href = "https://twitter.com/StevenXGe"
        ),
        " for recent updates. ",
        "File bug reports on ",
        a("GitHub.",
          href = "https://github.com/espors/idepGolem"
        )
      ),
      h3("NO GUARANTEE OF ACCURACY"),
      p("iDEP is developed by a small team with limited resources.
       We have not thoroughly tested it. So please verify all
       findings using other tools or R scripts. We tried our best
       to ensure our analysis is correct, but there is no guarantee.
      "),
      p(" By offering so many combinations of methods to analyze a data set,
      iDEP enables users to rationalize. It is human nature to focus on results
      that we like to see (confirmation bias).
      It is unfortunate that you can almost found further support
      for almost any theory
      from the massive but noisy literature.
      We encourage users to be critical of the results obtained using iDEP.
      Try to focus on robust results, rather than those that only should up
      with a certain parameter using a particular method.
      "),
      h3("Reproducibility"),
      p("Download the reports in the multiple
       tabs, which has the parameters and results. Also download and try
       the the R code on the DEG1 tab."),
      h4("Previous versions of iDEP are still useable:"),
      a("iDEP 0.96 with Ensembl Release 104, released July 30, 2022 ",
        href = "http://bioinformatics.sdstate.edu/idep96/"
      ),
      br(),
      a("iDEP 0.95 with Ensembl Release 104, released Feb. 8, 2022 ",
        href = "http://bioinformatics.sdstate.edu/idep95/"
      ),
      br(),
      a("iDEP 0.94 with Ensembl Release 104, released on Oct. 15, 2021 ",
        href = "http://bioinformatics.sdstate.edu/idep94/"
      ),
      br(),
      a("iDEP 0.93 with Ensembl Release 103, released May 20, 2021 ",
        href = "http://bioinformatics.sdstate.edu/idep93/"
      ),
      br(),
      a("iDEP 0.92 with Ensembl Release 100, released May 20, 2021 ",
        href = "http://bioinformatics.sdstate.edu/idep92/"
      ),
      br(),
      a("iDEP 0.90 with Ensembl Release 96, released May 19, 2021 ",
        href = "http://bioinformatics.sdstate.edu/idep90/"
      ),
      br(),
      a("iDEP 0.85 with Ensembl Release 95, released March 29, 2019 ",
        href = "http://bioinformatics.sdstate.edu/idep85/"
      ),
      br(),
      a("iDEP 0.82 with Ensembl  Release 92, released July 11, 2018 ",
        href = "http://bioinformatics.sdstate.edu/idep82/"
      ),
      br(),
      a("iDEP 0.73 with Ensembl  Release 91,  released December 2017 ",
        href = "http://bioinformatics.sdstate.edu/idep73/"
      ),
      htmlOutput(ns("session_info"))
    )
  )
}

#' 10_doc Server Functions
#'
#' @noRd
mod_10_doc_server <- function(id, pre_process, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$session_info <- renderUI({
      i <- c("<br><h4>R session info: </h4>")
      i <- c(i, capture.output(sessionInfo()))
      HTML(paste(i, collapse = "<br/>"))
    })
  })
}

## To be copied in the UI
# mod_09_doc_ui("10_doc_ui_1")

## To be copied in the server
# mod_09_doc_server("10_doc_ui_1")
