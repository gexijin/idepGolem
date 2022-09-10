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
      h3("Letters of support always welcome"),
      p("iDEP is developed and maintained by a small team. 
      If you find iDEP helpful, please send us a brief  ",
        a(
          "email.",
          href = "mailto:gelabinfo@gmail.com?Subject=iDEP"
        ),
        " These emails will help us secure the next round 
        of funding to sustain and improve this tool. 
        Currently, we are supported by a grant from 
        NIH/NHGRI (R01HG010805)."
      ),
      p("Our small team consists 
      of Xijin Ge (PI), Jianli Qi(research staff), and several 
      graduate students. 
      Graduate students who are currently working on iDEP include
       Emma Spors and Ben Derenge.
      Past students include Eun Wo Son, Runan Yao,
      Gavin Doering, Roberto Villegas-Diaz, and Eric Tulowetzke. 
      Eric still helps us fix bugs after leaving the lab. 
      Much of the new version of iDEP is rewritten by Gavin Doering. 
      The iDEP logo was designed by Emma Spors.
      "),
      h3("Citation:"),
      p("Ge, Son & Yao, iDEP: an integrated web application for differential 
        expression and pathway analysis of RNA-Seq data, ", 
        a("BMC Bioinformatics 19:1-24, 2018.", 
           href="https://doi.org/10.1186/s12859-018-2486-6",
           target="_blank"
        )
      ),
      p("Consider citing other tools that form the foundation of iDEP, such as ", 
        a("ENSEMBL, ", 
           href="https:/doi.org/10.1093/nar/gkab1049",
           target="_blank"
        ),
        a("STRING-db,", 
           href="https://doi.org/10.1093/nar/gky1131",
           target="_blank"
        ),
        a("DESeq2", 
           href="https://doi.org/10.1186/s13059-014-0550-8",
           target="_blank"
        ),
        " and many others.",
        " If you use the KEGG diagram, please also cite ", 
        a("pathview, ", 
           href="https://doi.org/10.1093/bioinformatics/btt285",
           target="_blank"
        ),
        "and ",
        a("KEGG.", 
           href="https://doi.org/10.1093/nar/gkaa970",
           target="_blank"
        )
      ),
      h4(
        "Source code on ",
        a(
          "GitHub.",
          href = "https://github.com/espors/idepGolem/"
        )
      )
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
