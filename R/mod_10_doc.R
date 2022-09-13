#' 10_doc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_10_doc_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    title = "About",
    fluidPage(
      h3("If you find iDEP helpful, please ",
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
        was essential when we applied for the current 
        grant from NIH/NHGRI (R01HG010805), 
        which expires in 20 months. Your letters will help sustain 
        and improve this service."
      ),
      p("iDEP is developed and maintained by a small team at ",
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
      h3("Citation"),
      p("Ge, Son & Yao, iDEP: an integrated web application for differential 
        expression and pathway analysis of RNA-Seq data, ",
        a("BMC Bioinformatics 19:1-24, 2018.", 
           href = "https://doi.org/10.1186/s12859-018-2486-6",
           target = "_blank"
        )
      ),
      p("Consider citing other tools that form the foundation of iDEP, such as ", 
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
      h3("Source code and database"),
      p("Source code is available on ",
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
      h3("NO GUARANTEE OF ACCURACY"),
      p("iDEP is developed by a small team with limited resources.
       We have not thoroughly tested it. So please verify all 
       findings using other tools or R scripts. We tried our best 
       to ensure our analysis is correct, but there is no guarantee.
      "),
      p(" By offering so many combinations of methods to analyze a data set,
      iDEP enables people to pick up the results that like to see (confirmation bias).
      It is unfortunate that you can almost found further support for almost any theory 
      from the massive but noisy literature. 
      We encourage users to be critical of the results obtained using iDEP.
      Try to focus on robust results, rather than those that only should up
      with a certain parameter using a particular method. 
      ")
    )
  )

}
    
#' 10_doc Server Functions
#'
#' @noRd 
mod_10_doc_server <- function(id, pre_process, idep_data, tab){
  moduleServer( id, function(input, output, session){
    ns <- session$ns 
  })
}
    
## To be copied in the UI
# mod_09_doc_ui("10_doc_ui_1")
    
## To be copied in the server
# mod_09_doc_server("10_doc_ui_1")
