# Author details
## Script name: app.R
## Purpose of the script: Run Transcriptomic Analysis for all species in BioMart.
## Author(s): Jyotirmoy Das, Ph.D.
## Date Created: 2023-03-01
## Date Last Modified: 2023-09-06
## Copyright statement: eula.txt
## Contact information: jyotirmoy.das@liu.se
## Please cite: 
## @article{},
## Notes: 
## github: https://github.com/JD2112/transcriptr/tree/main
## mannual: 
## Please read the manual/tutorial file for the sample preparation of the analysis.


## Load libraries
suppressMessages(suppressWarnings(library(shiny)))
suppressMessages(suppressWarnings(library(shinyWidgets)))
suppressMessages(suppressWarnings(library(shinydashboard)))
suppressMessages(suppressWarnings(library(shinydashboardPlus)))
suppressMessages(suppressWarnings(library(shinyjs)))
suppressMessages(suppressWarnings(library(Cairo)))
suppressMessages(suppressWarnings(library(viridisLite)))
suppressMessages(suppressWarnings(library(DT)))
suppressMessages(suppressWarnings(library(chromoMap)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(require(DO.db)))
suppressWarnings(suppressMessages(require(DOSE)))
suppressWarnings(suppressMessages(require(GO.db)))
suppressWarnings(suppressMessages(require(GOSemSim)))
suppressWarnings(suppressMessages(require(clusterProfiler)))
suppressWarnings(suppressMessages(require(org.Hs.eg.db)))
suppressWarnings(suppressMessages(require(org.Mm.eg.db)))
suppressWarnings(suppressMessages(require(org.Rn.eg.db)))
suppressWarnings(suppressMessages(require(org.Dm.eg.db)))
suppressWarnings(suppressMessages(require(org.At.tair.db)))
suppressWarnings(suppressMessages(require(org.Dr.eg.db)))
suppressWarnings(suppressMessages(require(org.Sc.sgd.db)))
suppressWarnings(suppressMessages(require(org.Ce.eg.db)))
suppressWarnings(suppressMessages(require(org.Bt.eg.db)))
suppressWarnings(suppressMessages(require(org.Gg.eg.db)))
suppressWarnings(suppressMessages(require(org.Ss.eg.db)))
suppressWarnings(suppressMessages(require(org.Mmu.eg.db)))
suppressWarnings(suppressMessages(require(org.Cf.eg.db)))
suppressWarnings(suppressMessages(require(org.EcK12.eg.db)))
suppressWarnings(suppressMessages(require(org.Xl.eg.db)))
suppressWarnings(suppressMessages(require(org.Pt.eg.db)))
suppressWarnings(suppressMessages(require(org.Ag.eg.db)))
suppressWarnings(suppressMessages(require(org.EcSakai.eg.db)))
suppressWarnings(suppressMessages(require(org.Mxanthus.db)))
suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(ReactomePA)))
suppressWarnings(suppressMessages(library(htmlwidgets)))
suppressWarnings(suppressMessages(library(shinyalert)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(shinyFiles))
suppressMessages(suppressMessages(library(FactoMineR)))
suppressMessages(suppressMessages(library(factoextra)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(UpSetR)))
suppressWarnings(suppressMessages(library(Vennerable)))
suppressWarnings(suppressMessages(require(AnnotationDbi)))
suppressWarnings(suppressMessages(require(Biobase)))
suppressWarnings(suppressMessages(require(BiocGenerics)))
suppressWarnings(suppressMessages(library(ggforce)))
suppressMessages(suppressMessages(library(ggfortify)))
suppressMessages(suppressMessages(library(explor)))
suppressMessages(suppressMessages(library(qpdf)))
suppressMessages(suppressMessages(library(ggpubr)))
suppressMessages(suppressMessages(library(shinylogs)))
suppressMessages(suppressMessages(library(sf)))
suppressMessages(suppressMessages(library(rgdal)))
suppressMessages(suppressWarnings(library(shinythemes)))
suppressMessages(suppressWarnings(library(shinyBS)))
suppressMessages(suppressWarnings(library(zip)))
suppressMessages(suppressWarnings(library(pathview)))
suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(spsComps)))
suppressMessages(suppressWarnings(library(systemPipeShiny)))
suppressMessages(suppressWarnings(library(spsUtil)))
suppressMessages(suppressWarnings(library(R.utils)))

# new libraries
suppressMessages(suppressWarnings(library(shinyhelper)))
suppressMessages(suppressWarnings(library(ggvenn)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(colourpicker)))
suppressMessages(suppressWarnings(library(crosstalk)))
suppressMessages(suppressWarnings(library(visNetwork)))
suppressMessages(suppressWarnings(library(gprofiler2)))
suppressMessages(suppressWarnings(library(mygene)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(igvShiny)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(later)))


#------- library ended ---------#


linebreaks <- function(n) {
  HTML(strrep(br(), n))
}

jscode <- "
shinyjs.init = function() {
    document.getElementById('heatmapplot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('distPlot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('gnattPlot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('volcanoplot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('pathplot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('cnvtumorplot').addEventListener('contextmenu', event => event.preventDefault());
    document.getElementById('piechart').addEventListener('contextmenu', event => event.preventDefault());
}"

js <- '
$(document).keyup(function(event) {
  if ($("#password").is(":focus") && (event.keyCode == 13)) {
      $("#ok").click();
  }
});
'
jsdis <- '
$(document).ready(function(){
  $("#myplot").on("shiny:recalculating", function(){
    $("#run").prop("disabled", true).text("Running...");
  }).on("shiny:recalculated", function(){
    $("#run").prop("disabled", false).text("Run");
  });
})
'

jscode1 <- HTML("
$(document).on('shiny:connected', function(event) {
  $('.sidebar-toggle').on('click', function() {
    if ($('body')[0].className != 'skin-blue sidebar-mini sidebar-collapse') {
      $('#sidebarCollapsed').css('display', 'none')
      $('nav.navbar-static-top').css('width', '1800px')
      $('nav.navbar-static-top').css('margin-left', '0px')
      $('header.main-header').css('width', '1800px')
      $('.sidebar-toggle').css('position', 'relative')
      $('span.logo').css('display', 'none')
    } else {
      $('#sidebarCollapsed').css('display', 'block')
      $('nav.navbar-static-top').css('margin-left', '230px')
      $('header.main-header').css('width', '884.44px')
      $('nav.navbar-static-top').css('width', '1800.44px')
      $('span.logo').css('display', 'block')
    }
  })
});
")

csscode <- HTML("
.sidebar-mini.sidebar-collapse .content-wrapper {
      margin-left: 0px !important;
}")

csscodebox <- HTML("
.box-collapse .content-wrapper {
      margin-left: 0px !important;
}")


modify_stop_propagation <- function(x) {
  x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}

sidebar <- dashboardSidebar(
  width = 230, collapsed = FALSE,
    tags$head(tags$script(jscode1)),
    tags$head(tags$style(csscode)),
  sidebarMenu(id="transcriptrmenu",
    menuItem("Dashboard", tabName = "main", icon=icon("database")),    
    modify_stop_propagation(
      menuItem("TranscriptR", startExpanded = TRUE, icon=icon(name="anchor", class= NULL, lib = "font-awesome"),
               menuSubItem("Quality Analysis", tabName = "qc", icon = icon("venus-mars")),
               menuSubItem("Transcriptome Analysis", tabName = "transana", icon = icon("venus-mars")),
               menuSubItem("IGViewer", tabName = "igviewer", icon = icon("chart-line")),
               menuSubItem("Dimensional Aanlysis", tabName = "pca", icon = icon("calculator")),
               menuSubItem("Volcano Plot", tabName = "volcano", icon = icon("jedi")),
               menuSubItem("HeatMap", tabName = "heatmap", icon = icon("ethernet")))
    ),
    modify_stop_propagation(
      menuItem("Enrichment Analysis", startExpanded = TRUE, icon=icon(name="anchor", class= NULL, lib = "font-awesome"),
               menuSubItem("Gene Ontology", tabName = "go", icon = icon("venus-mars")),
               menuSubItem("Pathway Analysis", tabName = "pathway", icon = icon("chart-line"))
               )
    ),
    modify_stop_propagation(
    menuItem("Set analysis", startExpanded = TRUE, icon=icon(name="anchor", class= NULL, lib = "font-awesome"),
      menuSubItem("Venn Analysis", tabName = "venn", icon = icon("venus-mars")),
      menuSubItem("UpSet Plot", tabName = "upset", icon = icon("chart-line")))
    ),
    menuItem("User Manual", tabName = "manual", icon = icon("file-pdf"))
  ),
  div(class="main-footer",
      style = " margin-left: 0;
              height: -30px;
              font-size: 15px;
              position: absolute;
              bottom: 0;
              left: 0;
              right: 0;
              text-align: center;
              font-family: Lucida Console;
              padding:0px;
              z-index: 1000;
              color: white"
  )
)


bodyMain <- tabItem(tabName= "main",
                    h2("Dashboard"),
                    box(
                        title = "Disclaimer",
                        status = "navy",
                        width = 12,
                        #height = "380px",
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        collapsed = FALSE,
                        closable = FALSE,
                        div(style = "font-size:16px; color:orange",
                            HTML("<b>For non-commercial Academic and Research purpose only!</b>")
                        ),
                        div(style = "font-size:16px; color:black",
                            HTML("Here we introduce TranscriptR, a complete pipeline for the transcriptome analysis which not only offers data visualization and normalization but also provide additional features such as the annotation of the genomic features resulting from the analysis, pairwise comparisons of DEGs with different graphical representation plus functional and pathway enrichment as downstream analysis directly from the sequencing data, all packed in a minimal, elegant and intuitive graphical user interface which brings the analysis of transcriptome data")
                        ),
                        div(style = "font-size:16px; font-weight: bold; color:orange",
                            HTML("NOTICE: YOUR INPUT DATA FILES AND RESULT WILL BE DELETED AS SOON AS YOU CLOSE THE BROWSER. WE DONOT READ/STORE YOUR DATA AND NOT RESPONSIBLE FOR IT.")
                        ),
                        div(style = "font-size:16px; font-weight: bold; color:red; font-family: monospace",
                        HTML("WARNING - Please note we included the possibility of the human transcriptome data analysis from the sequence data. None of the owner of the code or data center is responsible for handling the human data. It is upto the user if they want to use the online version for human data analysis."), br(),
                        tags$a(href = "https://hub.docker.com/jd21/transcriptr", "Otherwise the user can download the docker container and use it.", target= "_blank"), br(), 
                        tags$a(href = "https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiK1t-O35uBAxUCQ_EDHaXdCWoQFnoECCsQAQ&url=https%3A%2F%2Fwww.phgfoundation.org%2Fmedia%2F123%2Fdownload%2Fgdpr-and-genomic-data-report.pdf%3Fv%3D1&usg=AOvVaw0UnPevgJgQUtlOm3SUuOx7&opi=89978449", "Please read more here", target= "_blank")
                        )
                        ),
                    
                    tags$style(HTML("
                      .box-header {
                        padding: 0 10px 0 0;
                      }
                      .box-header h3 {
                        width: 100%;
                        padding: 10px;
                      }")),

                    fluidRow(
                      useShinyjs(),
                      extendShinyjs(text = jscode, functions = c()),

                      box(
                        title = "Notification",
                        status = "navy",
                        width = 6,
                        height = "380px",
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        collapsed = FALSE,
                        closable = FALSE,
                        boxToolSize = "lg",
                        htmlOutput("frame")
                        ),

                      box(
                        title = "Main Packages",
                        status = "navy",
                        width = 6,
                        height = "400px",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        collapsed = FALSE,
                        closable = FALSE,
                        boxToolSize = "xs",
                        userList(
                          userListItem(
                            image = "shiny.png", 
                            title = "Shiny"
                          ),                          
                          userListItem(
                            image = "tidyr.png", 
                            title = "tidyr"
                          ),
                          userListItem(
                            image = "tidyverse.png", 
                            title = "tidyverse"
                          ),
                          userListItem(
                            image = "ggplot2.png", 
                            title = "ggplot2"
                          ),
                          userListItem(
                            image = "dplyr.png", 
                            title = "dplyr"
                          ),
                          userListItem(
                            image = "purrr.png", 
                            title = "purrr"
                          ),
                          userListItem(
                            image = "readr.png", 
                            title = "readr"
                          ),
                          userListItem(
                            image = "readxl.png", 
                            title = "readxl"
                          )                          
                        )
                      ),
                      box(
                        title = "Supported by",
                        status = "navy",
                        width = 12,
                        height = "400px",
                        solidHeader = FALSE,
                        collapsible = FALSE,
                        collapsed = FALSE,
                        closable = FALSE,
                        boxToolSize = "xs",
                        column(4,
                         tags$a(img(src= "LiU-primary-logo.svg"), href="https://liu.se/en", target="_blank", align = "center")
                        ),
                        column(4,
                        br(),
                        br(),
                        br(),
                         tags$a(img(src= "core-facility-logo.svg", width = "80%"), href="https://liu.se/en/organisation/liu/medfak/coref", target="_blank", align = "center")
                         ),
                        column(4,
                         tags$a(img(src= "SciLifeLab_Logotype_Green_POS.svg"), href= "https://www.scilifelab.se/units/clinical-genomics-linkoping/", target="_blank", align = "center")
                        )
                      )
                    ))
                    

bodyTranscrip <- tabItem(tabName = "transana",                    
                    fluidRow(
                      useShinyjs(),
                      
                        extendShinyjs(text = jscode, functions = c()),
                        tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
                      h3("Transcriptomic Analysis") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptomicshelp"),
                      tags$head(
                        tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                           .inline .form-group{display: table-row;}")
                      ),
                      tags$head(
                        tags$style(
                          "#inTabset {
                            position: fixed;
                            width: 100%;
                            background-color: white;
                            top: 0;
                            }",
                          ".tab-content  {
                            margin-top: 10px;
                            margin-left: 10px;
                          }"
                        )
                      ),                      
                        tags$style(HTML("
                        .tabbable > .nav > li[class=active] > a {
                          background-color: #012D5A;
                          color: #FFF;
                          border-color: #012D5A;
                        }")),
                tabsetPanel(
                        tabPanel("Data Input and Parameter settings",
                        br(),                        
                        tabsetPanel(
                          id = "tabs",
                        tabPanel(value = "fqupload", 
                          "Fastq upload",
                        br(),
                          box(
                            title = "Input data files",
                            closable = FALSE,
                            width = 3,
                            status = "navy",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                            tags$tr(
                              width = "100%",
                              tags$td(width = "40%", tags$div(style = "font-size:14pX;", h4("Data Selection")) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "datauploadhelp")),
                              br(),
                              tags$td(
                                width = "60%",
                                awesomeRadio("datasel",
                                  label = NULL,
                                  choices = c(
                                    "Continue from Quality Analysis" = "qcdata",
                                    "Upload data" = "dataupload"
                                  ),
                                  selected = "qcdata",
                                  status = "success"
                                ),
                                conditionalPanel(
                                  condition = "input.datasel == 'dataupload' ",
                                  fileInput("file_transcriptome",
                                    label = "Upload FASTQ files (datafiles)",
                                    accept = c(".fastq.gz", ".fq.gz", ".gz"),
                                    multiple = TRUE
                                  )
                                )
                              )
                            )
                          ),
                          box(
                            title = "Data file viewer",
                            closable = FALSE,
                            width = 3,
                            status = "navy",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                            verbatimTextOutput("fastq_files") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "datafileviewerhelp")
                          ),
                          div(                            
                            style = "float: inherit;
                                  color: #fff;                                   
                                  border-color: #2e6da4",
                          actionBttn("next1", label = NULL, 
                            size = "md",
                            #color = "royal",
                            icon = icon("arrow-right")
                            #style = "gradient"
                            ))
                          ),
                          tabPanel("Sample info", value = "sinfo",
                        box(
                          title = "Sample Sheet/Meta Data file",
                          closable = FALSE,
                          width = 3,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          fileInput("file_transmeta",
                                    label = "Select Sample information file",
                                    accept = ".txt",
                                    multiple = FALSE
                                    ) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "sampleinfohelp"),
                         DT::dataTableOutput("transsamplefile") 
                        ),
                        box(
                          title = "Comparison file",
                          closable = FALSE,
                          width = 3,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          fileInput("file_transcomp",
                                    label = "Select Comparison file",
                                    accept = ".txt",
                                    multiple = FALSE
                                    ) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "samplemetahelp"),
                         DT::dataTableOutput("transcompfile") 
                        ),
                        div(
                          #style = "position:absolute;bottom:10em;left: 65em",
                          style = "float: inherit;
                                  color: #fff;                                   
                                  border-color: #2e6da4",
                        actionBttn("next2", label = NULL, 
                          size = "md",                          
                          #color = "royal",
                          icon = icon("arrow-right"),
                          #style = "gradient"
                          )),
                        div(
                          #style = "position:absolute;bottom:10em;left: 60em",
                          style = "float: inherit;
                                  color: #fff;                                   
                                  border-color: #2e6da4",
                          actionBttn("back1", label = NULL, 
                            size = "md",
                            #color = "royal",
                            icon = icon("arrow-left")
                            #style = "gradient"
                            #style=""
                            ))
                        ),
                        #linebreaks(15),
                        tabPanel("Params settings", value = "paramss",
                        box(
                          title = "Parameter Configurations setup",
                          closable = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          h4("Please select the parameters below to run the analysis pipeline.") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "parametershelp"),
                          
                        tags$table(
                          tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select sample type")),
                            tags$td(width = "60%", 
                            awesomeRadio("endtype",
                                      label = NULL,
                                      choices = c(
                                        "Single-end" = "send",
                                        "Paired-end" = "pend"
                                      ), 
                                      selected = "pend",
                                      status = "success"),
                            )),
                        tags$tr(width="100%",
                          tags$td(width = "40%", tags$div(style = "font-size:18pX; color: #354319;", "Species and database") 
                          )
                          ),      
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select Species")),
                            tags$td(width = "60%", selectInput("species",
                                      label = NULL,
                                      choices = c(
                                      "Human" = "homo_sapiens",
                                      "Mouse" = "mus_musculus",
                                      "Rat" = "rattus_norvegicus",
                                      "----------------------------",
                                      "acanthochromis_polyacanthus",
                                      "accipiter_nisus",
                                      "ailuropoda_melanoleuca",
                                      "amazona_collaria",
                                      "amphilophus_citrinellus",
                                      "amphiprion_ocellaris",
                                      "amphiprion_percula",
                                      "anabas_testudineus",
                                      "anas_platyrhynchos_platyrhynchos",
                                      "anas_platyrhynchos",
                                      "anas_zonorhyncha",
                                      "anolis_carolinensis",
                                      "anser_brachyrhynchus",
                                      "anser_cygnoides",
                                      "aotus_nancymaae",
                                      "apteryx_haastii",
                                      "apteryx_owenii",
                                      "apteryx_rowi",
                                      "aquila_chrysaetos_chrysaetos",
                                      "astatotilapia_calliptera",
                                      "astyanax_mexicanus_pachon",
                                      "astyanax_mexicanus",
                                      "athene_cunicularia",
                                      "balaenoptera_musculus",
                                      "betta_splendens",
                                      "bison_bison_bison",
                                      "bos_grunniens",
                                      "bos_indicus_hybrid",
                                      "bos_mutus",
                                      "bos_taurus_hybrid",
                                      "bos_taurus",
                                      "bubo_bubo",
                                      "buteo_japonicus",
                                      "caenorhabditis_elegans",
                                      "cairina_moschata_domestica",
                                      "calidris_pugnax",
                                      "calidris_pygmaea",
                                      "callithrix_jacchus",
                                      "callorhinchus_milii",
                                      "camarhynchus_parvulus",
                                      "camelus_dromedarius",
                                      "canis_lupus_dingo",
                                      "canis_lupus_familiaris",
                                      "canis_lupus_familiarisbasenji",
                                      "canis_lupus_familiarisboxer",
                                      "canis_lupus_familiarisgreatdane",
                                      "canis_lupus_familiarisgsd",
                                      "capra_hircus_blackbengal",
                                      "capra_hircus",
                                      "carassius_auratus",
                                      "carlito_syrichta",
                                      "castor_canadensis",
                                      "catagonus_wagneri",
                                      "catharus_ustulatus",
                                      "cavia_aperea",
                                      "cavia_porcellus",
                                      "cebus_imitator",
                                      "cercocebus_atys",
                                      "cervus_hanglu_yarkandensis",
                                      "chelonoidis_abingdonii",
                                      "chelydra_serpentina",
                                      "chinchilla_lanigera",
                                      "chlorocebus_sabaeus",
                                      "choloepus_hoffmanni",
                                      "chrysemys_picta_bellii",
                                      "chrysolophus_pictus",
                                      "ciona_intestinalis",
                                      "ciona_savignyi",
                                      "clupea_harengus",
                                      "colobus_angolensis_palliatus",
                                      "corvus_moneduloides",
                                      "cottoperca_gobio",
                                      "coturnix_japonica",
                                      "cricetulus_griseus_chok1gshd",
                                      "cricetulus_griseus_crigri",
                                      "cricetulus_griseus_picr",
                                      "crocodylus_porosus",
                                      "cyanistes_caeruleus",
                                      "cyclopterus_lumpus",
                                      "cynoglossus_semilaevis",
                                      "cyprinodon_variegatus",
                                      "cyprinus_carpio_carpio",
                                      "cyprinus_carpio_germanmirror",
                                      "cyprinus_carpio_hebaored",
                                      "cyprinus_carpio_huanghe",
                                      "cyprinus_carpio",
                                      "danio_rerio",
                                      "dasypus_novemcinctus",
                                      "delphinapterus_leucas",
                                      "denticeps_clupeoides",
                                      "dicentrarchus_labrax",
                                      "dipodomys_ordii",
                                      "dromaius_novaehollandiae",
                                      "drosophila_melanogaster",
                                      "echeneis_naucrates",
                                      "echinops_telfairi",
                                      "electrophorus_electricus",
                                      "eptatretus_burgeri",
                                      "equus_asinus",
                                      "equus_caballus",
                                      "erinaceus_europaeus",
                                      "erpetoichthys_calabaricus",
                                      "erythrura_gouldiae",
                                      "esox_lucius",
                                      "falco_tinnunculus",
                                      "felis_catus",
                                      "ficedula_albicollis",
                                      "fukomys_damarensis",
                                      "fundulus_heteroclitus",
                                      "gadus_morhua",
                                      "gallus_gallus_gca000002315v5",
                                      "gallus_gallus_gca016700215v2",
                                      "gallus_gallus",
                                      "gambusia_affinis",
                                      "gasterosteus_aculeatus",
                                      "geospiza_fortis",
                                      "gopherus_agassizii",
                                      "gopherus_evgoodei",
                                      "gorilla_gorilla",
                                      "gouania_willdenowi",
                                      "haplochromis_burtoni",
                                      "heterocephalus_glaber_female",
                                      "heterocephalus_glaber_male",
                                      "hippocampus_comes",
                                      "hucho_hucho",
                                      "ictalurus_punctatus",
                                      "ictidomys_tridecemlineatus",
                                      "jaculus_jaculus",
                                      "junco_hyemalis",
                                      "kryptolebias_marmoratus",
                                      "labrus_bergylta",
                                      "larimichthys_crocea",
                                      "lates_calcarifer",
                                      "laticauda_laticaudata",
                                      "latimeria_chalumnae",
                                      "lepidothrix_coronata",
                                      "lepisosteus_oculatus",
                                      "leptobrachium_leishanense",
                                      "lonchura_striata_domestica",
                                      "loxodonta_africana",
                                      "lynx_canadensis",
                                      "macaca_fascicularis",
                                      "macaca_mulatta",
                                      "macaca_nemestrina",
                                      "malurus_cyaneus_samueli",
                                      "manacus_vitellinus",
                                      "mandrillus_leucophaeus",
                                      "marmota_marmota_marmota",
                                      "mastacembelus_armatus",
                                      "maylandia_zebra",
                                      "meleagris_gallopavo",
                                      "melopsittacus_undulatus",
                                      "meriones_unguiculatus",
                                      "mesocricetus_auratus",
                                      "microcebus_murinus",
                                      "microtus_ochrogaster",
                                      "mola_mola",
                                      "monodelphis_domestica",
                                      "monodon_monoceros",
                                      "monopterus_albus",
                                      "moschus_moschiferus",
                                      "mus_caroli",
                                      "mus_musculus_129s1svimj",
                                      "mus_musculus_aj",
                                      "mus_musculus_akrj",
                                      "mus_musculus_balbcj",
                                      "mus_musculus_c3hhej",
                                      "mus_musculus_c57bl6nj",
                                      "mus_musculus_casteij",
                                      "mus_musculus_cbaj",
                                      "mus_musculus_dba2j",
                                      "mus_musculus_fvbnj",
                                      "mus_musculus_lpj",
                                      "mus_musculus_nodshiltj",
                                      "mus_musculus_nzohlltj",
                                      "mus_musculus_pwkphj",
                                      "mus_musculus_wsbeij",
                                      "mus_pahari",
                                      "mus_spicilegus",
                                      "mus_spretus",
                                      "mustela_putorius_furo",
                                      "myotis_lucifugus",
                                      "myripristis_murdjan",
                                      "naja_naja",
                                      "nannospalax_galili",
                                      "neogobius_melanostomus",
                                      "neolamprologus_brichardi",
                                      "neovison_vison",
                                      "nomascus_leucogenys",
                                      "notamacropus_eugenii",
                                      "notechis_scutatus",
                                      "nothobranchius_furzeri",
                                      "nothoprocta_perdicaria",
                                      "numida_meleagris",
                                      "ochotona_princeps",
                                      "octodon_degus",
                                      "oncorhynchus_kisutch",
                                      "oncorhynchus_mykiss",
                                      "oncorhynchus_tshawytscha",
                                      "oreochromis_aureus",
                                      "oreochromis_niloticus",
                                      "ornithorhynchus_anatinus",
                                      "oryctolagus_cuniculus",
                                      "oryzias_javanicus",
                                      "oryzias_latipes_hni",
                                      "oryzias_latipes_hsok",
                                      "oryzias_latipes",
                                      "oryzias_melastigma",
                                      "oryzias_sinensis",
                                      "otolemur_garnettii",
                                      "otus_sunia",
                                      "ovis_aries_rambouillet",
                                      "ovis_aries",
                                      "pan_paniscus",
                                      "pan_troglodytes",
                                      "panthera_leo",
                                      "panthera_pardus",
                                      "panthera_tigris_altaica",
                                      "papio_anubis",
                                      "parambassis_ranga",
                                      "paramormyrops_kingsleyae",
                                      "parus_major",
                                      "pavo_cristatus",
                                      "pelodiscus_sinensis",
                                      "pelusios_castaneus",
                                      "periophthalmus_magnuspinnatus",
                                      "peromyscus_maniculatus_bairdii",
                                      "petromyzon_marinus",
                                      "phascolarctos_cinereus",
                                      "phasianus_colchicus",
                                      "phocoena_sinus",
                                      "physeter_catodon",
                                      "piliocolobus_tephrosceles",
                                      "podarcis_muralis",
                                      "poecilia_formosa",
                                      "poecilia_latipinna",
                                      "poecilia_mexicana",
                                      "poecilia_reticulata",
                                      "pogona_vitticeps",
                                      "pongo_abelii",
                                      "procavia_capensis",
                                      "prolemur_simus",
                                      "propithecus_coquereli",
                                      "pseudonaja_textilis",
                                      "pteropus_vampyrus",
                                      "pundamilia_nyererei",
                                      "pygocentrus_nattereri",
                                      "rhinolophus_ferrumequinum",
                                      "rhinopithecus_bieti",
                                      "rhinopithecus_roxellana",
                                      "saccharomyces_cerevisiae",
                                      "saimiri_boliviensis_boliviensis",
                                      "salarias_fasciatus",
                                      "salmo_salar",
                                      "salmo_trutta",
                                      "salvator_merianae",
                                      "sander_lucioperca",
                                      "sarcophilus_harrisii",
                                      "sciurus_vulgaris",
                                      "scleropages_formosus",
                                      "scophthalmus_maximus",
                                      "serinus_canaria",
                                      "seriola_dumerili",
                                      "seriola_lalandi_dorsalis",
                                      "sinocyclocheilus_anshuiensis",
                                      "sinocyclocheilus_grahami",
                                      "sinocyclocheilus_rhinocerous",
                                      "sorex_araneus",
                                      "sparus_aurata",
                                      "spermophilus_dauricus",
                                      "sphaeramia_orbicularis",
                                      "sphenodon_punctatus",
                                      "stachyris_ruficeps",
                                      "stegastes_partitus",
                                      "strigops_habroptila",
                                      "strix_occidentalis_caurina",
                                      "struthio_camelus_australis",
                                      "suricata_suricatta",
                                      "sus_scrofa_bamei",
                                      "sus_scrofa_berkshire",
                                      "sus_scrofa_hampshire",
                                      "sus_scrofa_jinhua",
                                      "sus_scrofa_landrace",
                                      "sus_scrofa_largewhite",
                                      "sus_scrofa_meishan",
                                      "sus_scrofa_pietrain",
                                      "sus_scrofa_rongchang",
                                      "sus_scrofa_tibetan",
                                      "sus_scrofa_usmarc",
                                      "sus_scrofa_wuzhishan",
                                      "sus_scrofa",
                                      "taeniopygia_guttata",
                                      "takifugu_rubripes",
                                      "terrapene_carolina_triunguis",
                                      "tetraodon_nigroviridis",
                                      "theropithecus_gelada",
                                      "tupaia_belangeri",
                                      "tursiops_truncatus",
                                      "urocitellus_parryii",
                                      "ursus_americanus",
                                      "ursus_maritimus",
                                      "ursus_thibetanus_thibetanus",
                                      "varanus_komodoensis",
                                      "vicugna_pacos",
                                      "vombatus_ursinus",
                                      "vulpes_vulpes",
                                      "xenopus_tropicalis",
                                      "xiphophorus_couchianus",
                                      "xiphophorus_maculatus",
                                      "zalophus_californianus",
                                      "zonotrichia_albicollis",
                                      "zosterops_lateralis_melanops"
                                      ),
                                      selected = 'homo_sapiens',
                                      multiple = FALSE))),
                        br(),
                        br(),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select release version")),
                            tags$td(width = "60%", 
                            awesomeRadio("dbversion",
                                      label = NULL,
                                      choices = c(
                                        "Latest" = "latest",
                                        "Custom" = "custom"
                                      ), 
                                      selected = "latest",
                                      status = "success"),
                          conditionalPanel(
                          condition = "input.dbversion == 'custom' ",
                          textInput("customdbversion",
                          "Enter version number",
                          "")
                        )
                        )),
                        #),
                          
                        br(),
                        #tags$table( 
                          tags$tr(width="100%",
                          tags$td(width = "40%", tags$div(style = "font-size:18pX; color: #354319;", "Alignment parameters"))), 
                          tags$tr(width = "100%",
                              tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select sample suffix/es")),
                              tags$td(width = "60%", 
                              conditionalPanel(
                            condition = "input.endtype == 'send' ",
                            textInput("R1suffix_se", label = "Choose R1 suffix", value = "_R1_001.fastq.gz"),
                          ),
                          conditionalPanel(
                            condition = "input.endtype == 'pend' ",
                            textInput("R1suffix_pe", label = "Choose R1 suffix", value = "_R1_001.fastq.gz"),
                            textInput("R2suffix_pe", label = "Choose R2 suffix", value = "_R2_001.fastq.gz")
                          )
                          )),                          
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Overhang")),
                            tags$td(width = "60%", numericInput("overh", label = NULL, value = 100))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose Sequence Type")),
                            tags$td(width = "60%", selectInput("seqtype", label = NULL, choices = c(
                              "Unstranded" = 0,
                              "Stranded" = 1,
                              "Reverse Stranded" = 2
                            ),
                            selected = 0,
                            multiple = FALSE))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose ID")),
                            tags$td(width = "60%", selectInput("attribute", label = NULL, choices = c(
                              "Gene ID" = "gene_id",
                              "Transcript ID" = "transcript_id"
                            ),
                            selected = "gene_id",
                            multiple = FALSE))),
                        tags$tr(width="100%",
                          tags$td(width = "40%", tags$div(style = "font-size:18pX; color: #354319;", "DEG filtration"))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Enter CPM value")),
                            tags$td(width = "60%", numericInput("cpm", label = NULL, value = 1, min = 0))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Number of samples")),
                            tags$td(width = "60%", numericInput("nsamples", label = NULL, value = 3,min = 1))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose logFC Value")),
                            tags$td(width = "60%", numericInput("logFC", label = NULL, value = 1.5, min = 0))),
                        tags$tr(width = "100%",
                            tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose FDR Value")),
                            tags$td(width = "60%", numericInput("FDR", label = NULL, value = 0.05, step = 0.01, min = 0, max = 1)))
                        )
                        ),
                        box(
                          title = "Parameter Configurations (JSON)",
                          closable = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          verbatimTextOutput("jsonview"),
                          downloadButton("download", "Download JSON")
                        ),
                        box(
                          title = "Run Analysis",
                          closable = FALSE,
                          width = 6,
                          status = "danger",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          actionButton(
                            "testtranscript",
                            "Submit Job",
                              width = "400px",
                            icon("paper-plane"), 
                            style="color: #fff; background-color: #000004; border-color: #2e6da4"
                          ) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptrrunhelp")
                        ),
                        div(
                          #style = "position:absolute;bottom:10em;left: 60em",
                          style = "float: inherit;
                                  color: #fff;                                   
                                  border-color: #2e6da4",
                          actionBttn("back2", label = NULL, 
                            size = "md",
                            #color = "royal",
                            icon = icon("arrow-left")#,
                            #style = "gradient"
                            ))
                            )
                        )
                        ),
                        tabPanel("MultiQC",
                          box(
                          title = "MultiQC result",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 12,
                          status = "navy", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "multiqc",
                            status = "danger",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                             actionButton("showqc", "Show MultiQC") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "multiqchelp"),                        
                             htmlOutput("qcplot", height = "800px", width = "100%")
                        )
                        ),
                        tabPanel("Plots",
                          box(
                          title = "",
                          closable = FALSE,
                          width = 12,
                          #status = "navy",
                          solidHeader = FALSE,
                          collapsible = FALSE,                          
                          div(class="text-center",
                          awesomeRadio("degplotsel",
                            label = "Choose analysis type",
                            choices = c(
                              "EdgeR analysis" = "edger",
                              "DESeq analysis" = "deseq",
                              "limma analysis" = "limma"
                            )
                          ) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptrplothelp")
                          )
                          ),
                          box(
                          title = "BoxPlots",
                          closable = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          actionButton("showboxplot", "Show BoxPlot"),
                          uiOutput("boxplottrans", height = 700, width = 400)
                        ),
                          box(
                          title = "MDSPlots",
                          closable = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          actionButton("showmdsplot", "Show MDS Plot"),
                          uiOutput("mdsplottrans", height = 700, width = 400)
                        )
                        ),
                        tabPanel("DEG table",
                          box(
                          title = "DEG Table",
                          closable = FALSE,
                          width = 9,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          column(4,
                          awesomeRadio("degtablesel",
                            label = "Choose analysis type",
                            choices = c(
                              "EdgeR analysis" = "edger",
                              "DESeq analysis" = "deseq",
                              "limma analysis" = "limma"
                            )
                            ) 
                            ),
                          column(4, 
                          conditionalPanel(
                            condition = "input.degtablesel == 'edger' ",
                            selectInput("showdegtab1",
                                    width = "400px",
                                    "Select DEG Table (only TSV files)",
                                    choices = list.files("/transcriptr/workdir/results/edgeR_results", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/edgeR_results", pattern = ".tsv")) #For local run
                          ),
                          conditionalPanel(
                            condition = "input.degtablesel == 'deseq' ",
                            selectInput("showdegtab2",
                                    width = "400px",
                                    "Select DEG Table (only TSV files)",
                                    choices = list.files("/transcriptr/workdir/results/deseq2_results", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/deseq2_results", pattern = ".tsv"))
                          ),
                          conditionalPanel(
                            condition = "input.degtablesel == 'limma' ",
                            selectInput("showdegtab3",
                                    width = "400px",
                                    "Select DEG Table (only TSV files)",
                                    choices = list.files("/transcriptr/workdir/results/limma_results", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/limma_results", pattern = ".tsv"))
                          )
                          
                                    ),
                          column(4, 
                          div(style = "display:inline-block;", 
                              numericInput("page", "Go to page:", value = 1, min = 1)),
                          div(style = "display:inline-block;", 
                              actionButton("gotopage", "Go")) %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptrdeghelp")),
                          DT::dataTableOutput("degtable")
                        ),
                        conditionalPanel(
                          condition = "input.species == 'homo_sapiens'",
                          box(
                          title = "Gene Information",
                          closable = FALSE,
                          width = 3,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          htmlOutput("geneframe", height="400px", width = "100%")                           
                        )
                        ),                        
                        box(
                          title = "CPM BoxPlot" %>%
                                  shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "cpmboxplothelp"),
                          closable = FALSE,
                          width = 3,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          plotlyOutput("cpmplot",height="300px", width = "100%")
                        )
                        ),
                        tabPanel("InterSectDE",
                          box(
                          title = "Intersection among DE results",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 12,
                          status = "navy", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "is",
                            status = "danger",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          column(4,
                          verbatimTextOutput("edgerPath"),
                          selectInput("showedger",
                                    "Select DEG Table from EdgeR results",
                                    choices = list.files("/transcriptr/workdir/results/edgeR_results/", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/edgeR_results", pattern = ".tsv"))
                                    ), # update later
                          column(4,
                          verbatimTextOutput("deseqPath"),
                          selectInput("showdeseq",
                                    "Select DEG Table from DESeq results",
                                    choices = list.files("/transcriptr/workdir/results/deseq2_results/", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/deseq2_results", pattern = ".tsv"))
                                    ), # update later
                          column(4,
                          verbatimTextOutput("limmaPath"),
                          selectInput("showlimma",
                                    "Select DEG Table from limma results",
                                    choices = list.files("/transcriptr/workdir/results/limma_results/", pattern = ".tsv")) # NOTE::: required for Docker!!!!
                                    #choices = list.files("/mnt/WD1/test/biomartr/results/limma_results", pattern = ".tsv"))                          
                                    ), # update later
                          actionButton("intersectdeg", "Show intersect result") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptrintersecthelp") 
                        ),
                        box(
                          title = "Intersect Table",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 8,
                          status = "navy", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "ist",
                            status = "danger",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          DT::dataTableOutput("insettable") # update later
                        ),
                        box(
                          title = "Intesect Figure",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 4,
                          status = "navy", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "isf",
                            status = "danger",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          plotOutput("insetplot", height="400px", width = "100%") # update later
                        ),
                        box(
                      title = "Download Plot",
                      closable = FALSE,
                      collapsed = FALSE,
                      width = 4,
                      status = "navy",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                              radioButtons(inputId = "filetype_invenn",
                              label = "Choose file type to download the plot:",
                              inline = TRUE,
                              choices = list("PDF", "PNG", "SVG", "TIFF")),
                            conditionalPanel(
                              condition = "input.filetype_invenn == 'PNG'",

                            numericInput("invenn_pngheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 10,
                                step = 1),

                            numericInput("invenn_pngwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 12,
                                step = 1),

                            selectInput("invenn_pngunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "inch"
                                ),

                            numericInput("invenn_pngresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 150,
                                max = 400,
                                step = 75)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_invenn == 'TIFF'",

                            numericInput("invenn_tiffheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 10,
                                step = 1),

                            numericInput("invenn_tiffwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 12,
                                step = 1),

                            selectInput("invenn_tiffunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "inch"
                                ),

                            numericInput("invenn_tiffresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 300)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_invenn == 'PDF'",

                            numericInput("invenn_pdfheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 20,
                                value = 10,
                                step = 1),

                            numericInput("invenn_pdfwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 20,
                                value = 12,
                                step = 1)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_invenn == 'SVG'",

                            numericInput("invenn_svgheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 20,
                                value = 10,
                                step = 1),

                            numericInput("invenn_svgwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 20,
                                value = 12,
                                step = 1)
                            ),
                            downloadButton(outputId = "invennDownload", lable = "Download Plot")
                    )
                        ),
                        tabPanel("Session & logs",
                        box(
                          title = "GTF Information",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          verbatimTextOutput("geninfo")
                        ),
                        box(
                          title = "Reference Information",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 6,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          verbatimTextOutput("othinfo")
                        ),
                        box(
                          title = "Shiny logs",
                          closable = FALSE,
                          collapsed = TRUE,
                          width = 12,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          #bsButton("consoleOut", "Refresh Console"),
                          downloadButton("downloadlog", "Download logFile"),
                          verbatimTextOutput("consoleText")
                          ),
                        box(
                          title = "Session Information",
                          closable = FALSE,
                          collapsed = TRUE,
                          width = 12,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          verbatimTextOutput("seninfo")
                        )
                        ),
                        tabPanel("Download results",
                          div(style = "font-size:14px; color:black",
                          tags$ol(HTML("1. User can save individual analysis result in the UI.
                          Result folder will contains BAM, BedGraphs, FastQC, MultiQC, Subread and EdgeR results")),
                          tags$ol(HTML("2. User can also download the result as a zip file for future analysis.")),
                          tags$ol(HTML("3. The zip folder will contains FastQC, MultiQC, Subread, log files and DE analysis results")),
                          h4("Download transcriptR results as zip file", tipify(actionButton("download_btn", icon("download")),
                          "",
                          placement = "bottom"))
                          
                          ))
                  )
                )
              )

bodyQC <- tabItem(
  tabName = "qc",
  fluidRow(
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c()),
    box(
      title = "Upload FASTQ files",
      closable = FALSE,
      width = 6,
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
      
   tags$table(
     tags$tr(width = "100%",
       tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select sample type")),
       tags$td(width = "60%", 
       awesomeRadio("qc_endtype",
                 label = NULL,
                 choices = c(
                   "Single-end" = "qc_send",
                   "Paired-end" = "qc_pend"
                 ), 
                 selected = "qc_pend",
                 status = "success"),
       )),      
      # tags$tr(
      #   width = "100%",
      #   tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose R1 suffix")),
      #   tags$td(width = "60%", textInput("R1suffix1", label = NULL, value = "_R1_001.fastq.gz"))
      # ),
      # tags$tr(
      #   width = "100%",
      #   tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Choose R2 suffix")),
      #   tags$td(width = "60%", textInput("R2suffix1", label = NULL, value = "_R2_001.fastq.gz"))
      # ),
      tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select sample suffix/es")),
                              tags$td(width = "60%", 
                              conditionalPanel(
                            condition = "input.qc_endtype == 'qc_send' ",
                            textInput("R1suffix_qc_se", label = "Choose R1 suffix", value = "_R1_001.fastq.gz"),
                          ),
                          conditionalPanel(
                            condition = "input.qc_endtype == 'qc_pend' ",
                            textInput("R1suffix_qc_pe", label = "Choose R1 suffix", value = "_R1_001.fastq.gz"),
                            textInput("R2suffix_qc_pe", label = "Choose R2 suffix", value = "_R2_001.fastq.gz")
                          )
                          )),
      fileInput("file_transcriptome1",
        label = "Upload FASTQ files (datafiles)",
        accept = c(".fastq.gz", ".fq.gz", ".gz"),
        multiple = TRUE
      ),
      actionButton("showqc1", "Run Quality Analysis",
        width = "100%",
        icon("paper-plane"),
        style = "color: #fff; background-color: #000004; border-color: #2e6da4"
      ),
      verbatimTextOutput("fastq_files1")
    ),
    box(
      title = "Parameter Configurations",
      closable = FALSE,
      width = 6,
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
      verbatimTextOutput("jsonview1"),
      downloadButton("downloadjson", "Download JSON")
    ),
    box(
      title = "MultiQC result",
      closable = FALSE,
      collapsed = FALSE,
      width = 12,
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      label = boxLabel(
        text = "multiqc",
        status = "danger",
        style = "default"
      ),
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
      actionButton("showqc2", "Show MultiQC"),                    
      htmlOutput("qcplot1", height = "800px", width = "100%")
    ),
    box(
                          title = "Shiny logs",
                          closable = FALSE,
                          collapsed = TRUE,
                          width = 12,
                          status = "navy",
                          solidHeader = TRUE,
                          collapsible = TRUE,
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          #bsButton("consoleOut", "Refresh Console"),
                          verbatimTextOutput("consoleTextqc")
                          ),
  )
)


bodyIGV <- tabItem(tabName = "igviewer",
                  h2("Genome Viewer") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "transcriptrigvhelp"),

                    fluidRow(
                      useShinyjs(),
                      extendShinyjs(text = jscode, functions = c()),
                      box(
                          title = "Parameters",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 3,
                          status = "info", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "params",
                            status = "primary",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                          #img(src = "igv.png", width = "100%")
                          #actionButton("searchButton", "Search"),
                          selectInput("bedfile", label = "Select BED file", 
                  choices = c("",
                    list.files("./transcriptr/workdir/results/bedgraphs/", pattern = "*.bedgraph")),
                    selected = NULL,
                    multiple = FALSE),
                  # selectInput("igvGenome", "Select genome",
                  # choices = c("hg38", "hg19", "mm10", "tair10", "rhos"),
                  # selected = NULL,
                  # multiple = FALSE),
                  # actionButton("addBedTrackButton", "Add as Bed"), br(),
                  actionButton("addBedGraphTrackButton", "Add as BedGraph"), br(),
                  actionButton("removeUserTracksButton", "Remove User Tracks"), br(),
                  actionButton("getChromLocButton", "Get Region"),
                  actionButton("clearChromLocButton", "Clear Region"),
                  div(style="background-color: white; width: 200px; height:30px; padding-left: 5px;
                              margin-top: 10px; border: 1px solid blue;",
                      htmlOutput("chromLocDisplay")),
                  actionButton("submitigv", "Run IGV")
                    ),
                      box(
                          title = "Genome View",
                          closable = FALSE,
                          collapsed = FALSE,
                          width = 9,
                          status = "navy", 
                          solidHeader = TRUE, 
                          collapsible = TRUE,
                          label = boxLabel(
                            text = "igv",
                            status = "danger",
                            style = "default"
                          ),
                          dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                      # tags$table(
                      #   tags$tr(width = "100%",
                      #       tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select Human Genome")),
                      #       tags$td(width = "50%", selectInput("hugenome",
                      #               label = NULL,
                      #               choices = c(
                      #                 "HG38" = "hg38",
                      #                 "HG19" = "hg19"
                      #               ),
                      #               selected = TRUE,
                      #               multiple = FALSE))),
                      #   tags$tr(width = "100%",
                      #       tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select BED file")),
                      #       tags$td(width = "50%", selectInput("bedfile",
                      #               label = NULL,
                      #               choices = list.files("./transcriptr/workdir/results/bedgraphs/", pattern = "*.bedgraph")))),
                      #   tags$tr(width = "100%",
                      #       tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select Track Color")),
                      #       tags$td(width = "50%", colourInput("trackcolo",
                      #               label = NULL,
                      #               "blue")))
                      # ),                  
                      # actionButton("submitigview", "Show Genome"),
                      # plotOutput("igview", height = "140px", width = "100%")
                      igvShinyOutput('igvShiny_0')#,
      # igvShinyOutput('igvShiny_1'),
      #width=10
                    )
                  )
                )

#====================================================#
## MDA module ####
#====================================================#
bodyMDA <- tabItem(tabName = "pca",
                   h2("Multi-Dimensional Analysis") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "multidimshelp"),
                   
                   fluidRow(                    
                    box(
                       title = "Data Upload", 
                       width = 4, 
                       status="navy",
                       solidHeader = TRUE,
                        collapsible = TRUE,
                        collapsed = FALSE,
                        closable = FALSE,
                        dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                        br(), 
                        tags$tr(
                          width = "100%",
                          tags$td(width = "40%", tags$div(style = "font-size:14pX;", h3("Select analysis"))),
                          br(),
                          tags$td(
                            width = "60%",
                            awesomeRadio("mdanalysis",
                              label = NULL,
                              choices = c(
                                "MDS" = "mds",
                                "PCA" = "pca"
                              ),
                              selected = "mds",
                              status = "success"
                            ),
                            conditionalPanel(
                              condition = "input.mdanalysis == 'mds' ",
                              tags$table(width="100%",
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Number of most variables")),
                                tags$td(width = "60%", numericInput("numdms", label = NULL, value = 1000, max = 10000))                                
                              ),                              
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select size")),
                                tags$td(width = "60%", numericInput("cexmds", label = NULL, value = 1, max = 10))                                
                              ),                              
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select X-dimension")),
                                tags$td(width = "60%", numericInput("dim1mds", label = NULL, value = 1, max = 10))                                
                              ), 
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select Y-dimension")),
                                tags$td(width = "60%", numericInput("dim2mds", label = NULL, value = 2, max = 10))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "X-axis label")),
                                tags$td(width = "60%", textInput("xlabmds", label = NULL, value = ""))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Y-axis label")),
                                tags$td(width = "60%", textInput("ylabmds", label = NULL, value = ""))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Upload data file")),
                                tags$td(width = "60%", fileInput('file_mds',
                                        label = NULL,
                                        accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              'text/tab-separated-values',
                                              '.csv',
                                              '.tsv',
                                              '.txt'
                                            )))
                              )),
                              actionButton("runmda", "Run MDS Analysis",
                                        width = "75%",
                                        icon("paper-plane"), 
                                        style = "color: #fff; background-color: #078082; border-color: #2e6da4; position:relative; left:calc(15%);"),
                              HTML("<hr> <a href='mda_example_data.txt' target='_blank'> <i class='fa fa-download'> </i> example data</a> | "),
                              HTML("<a href='mda_example_group.txt' target='_blank'> <i class='fa fa-download'> </i> group data</a>"),
                              awesomeCheckbox("mdscolour", 
                                  label = "Plot Colour", value = FALSE,
                                  status = "success"),
                                conditionalPanel(
                                condition = "input.mdscolour == true",
                                  fileInput('file_groupmds',
                                          label = "Upload the group file",
                                          accept = c(
                                              ".csv",
                                              ".txt"
                                            )),
                                  awesomeRadio("ncolourmds", 
                                              label= "Select group number for color", 
                                              choices = c("2", "3", "4", "5", "6", "7", "8", "9", "10"), selected = "2", inline = TRUE),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 2",
                                    colourInput("mdscol21", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol22", label = "Select color 2", "firebricks")
                                    ), 
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 3",
                                    colourInput("mdscol31", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol32", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol33", label = "Select color 3", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 4",
                                    colourInput("mdscol41", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol42", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol43", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol44", label = "Select color 4", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 5",
                                    colourInput("mdscol51", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol52", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol53", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol54", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol55", label = "Select color 5", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 6",
                                    colourInput("mdscol61", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol62", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol63", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol64", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol65", label = "Select color 5", "firebricks"),
                                    colourInput("mdscol66", label = "Select color 6", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 7",
                                    colourInput("mdscol71", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol72", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol73", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol74", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol75", label = "Select color 5", "firebricks"),
                                    colourInput("mdscol76", label = "Select color 6", "firebricks"),
                                    colourInput("mdscol77", label = "Select color 7", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 8",
                                    colourInput("mdscol81", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol82", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol83", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol84", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol85", label = "Select color 5", "firebricks"),
                                    colourInput("mdscol86", label = "Select color 6", "firebricks"),
                                    colourInput("mdscol87", label = "Select color 7", "firebricks"),
                                    colourInput("mdscol88", label = "Select color 8", "firebricks")
                                    ),
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 9",
                                    colourInput("mdscol91", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol92", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol93", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol94", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol95", label = "Select color 5", "firebricks"),
                                    colourInput("mdscol96", label = "Select color 6", "firebricks"),
                                    colourInput("mdscol97", label = "Select color 7", "firebricks"),
                                    colourInput("mdscol98", label = "Select color 8", "firebricks"),
                                    colourInput("mdscol99", label = "Select color 9", "firebricks")
                                    ), 
                                  conditionalPanel(
                                    condition = "input.ncolourmds == 10",
                                    colourInput("mdscol101", label = "Select color 1", "firebricks"),
                                    colourInput("mdscol102", label = "Select color 2", "firebricks"),
                                    colourInput("mdscol103", label = "Select color 3", "firebricks"),
                                    colourInput("mdscol104", label = "Select color 4", "firebricks"),
                                    colourInput("mdscol105", label = "Select color 5", "firebricks"),
                                    colourInput("mdscol106", label = "Select color 6", "firebricks"),
                                    colourInput("mdscol107", label = "Select color 7", "firebricks"),
                                    colourInput("mdscol108", label = "Select color 8", "firebricks"),
                                    colourInput("mdscol109", label = "Select color 9", "firebricks"),
                                    colourInput("mdscol1010", label = "Select color 10", "firebricks")
                                    ), 
                              )),
                            conditionalPanel(
                              condition = "input.mdanalysis == 'pca' ",
                              tags$table(width="100%",
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Number of most variables")),
                                tags$td(width = "60%", numericInput("numpca", label = NULL, value = 1000, max = 10000))                                
                              ),                              
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select size")),
                                tags$td(width = "60%", numericInput("cexpca", label = NULL, value = 1, max = 10))                                
                              ),                              
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Main title")),
                                tags$td(width = "60%", textInput("maintitlepca", label = NULL, value = ""))                                
                              ), 
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Legend title")),
                                tags$td(width = "60%", textInput("legendtitlepca", label = NULL, value = ""))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "X-axis label")),
                                tags$td(width = "60%", textInput("xlabpca", label = NULL, value = ""))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Y-axis label")),
                                tags$td(width = "60%", textInput("ylabpca", label = NULL, value = ""))                                
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Requires Mean Point")),
                                tags$td(width = "60%", selectInput("meanpca", label = NULL, choices = c("TRUE", "FALSE")))
                              ), 
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", "Select Colour Palette")),
                                tags$td(width = "60%", selectInput("colorpca", label = NULL, 
                                          choices = c(
                                            "NPG" = "npg",
                                            "AAAS" = "aaas",
                                            "Lancet" = "lancet",
                                            "JCO" = "jco",
                                            "UCSC GB" = "ucscgb",
                                            "U Chiago" = "uchicago",
                                            "Simpsons" = "simpsons",
                                            "Rick And Morty" = "rickandmorty"
                                            )))
                              ),
                              #),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", h5("Upload the data file"))),
                                tags$td(width = "60%", fileInput('file_pca',
                                        label = NULL,
                                        accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              'text/tab-separated-values',
                                              '.csv',
                                              '.tsv',
                                              '.txt'
                                            )))
                              ),
                              tags$tr(
                                width = "100%",
                                tags$td(width = "40%", tags$div(style = "font-size:14pX;", h5("Upload group file"))),
                                tags$td(width = "60%", fileInput('group_pca',
                                        label = NULL,
                                        accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              'text/tab-separated-values',
                                              '.csv',
                                              '.tsv',
                                              '.txt'
                                            )))
                              )
                            ),
                              actionButton("runpca", "Run PCA Analysis",
                              width = "75%",
                              icon("paper-plane"),
                              style = "color: #fff; background-color: #078082; border-color: #2e6da4; position:relative; left:calc(15%);"),
                              HTML("<hr> <a href='mda_example_data.txt' target='_blank'> <i class='fa fa-download'> </i> example data</a> | "),
                              HTML("<a href='mda_example_group.txt' target='_blank'> <i class='fa fa-download'> </i> group data</a>")                             
                            )
                          )
                        )),
                              
                          box(
                                title = "Analysis Result", 
                                width = 8, 
                                status="navy",
                                solidHeader = TRUE,
                                  collapsible = TRUE,
                                  collapsed = FALSE,
                                  closable = FALSE,
                                dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),      
                                  # DT::dataTableOutput("mdatest1"),
                                  # DT::dataTableOutput("mdatest2"),                              
                                  plotOutput("mdsPlot", width = "100%", height = "100%"),
                                  plotlyOutput("pcaplot", width = "100%", height = "100%")                                 
                          ),
                                  conditionalPanel(
                                    condition = "input.mdanalysis == 'mds' ",                             
                                  box(
                                    title = "Zoom and download", 
                                    status="navy",
                                    width = 8,
                                    solidHeader = TRUE,
                                    collapsible = TRUE,
                                    collapsed = FALSE,
                                    closable = FALSE,
                                    dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ), 
                                    sliderInput(
                                    "mds_size",
                                    label = "Zoom in/out MDS diagram",
                                    value = 500,
                                    min = 200,
                                    max = 1200,
                                    ticks = TRUE,
                                    step = 10),
                                    radioButtons(
                                      inputId = "filetype_mds",
                                      label = "Choose file type to download the plot:",
                                      inline = TRUE,
                                      choices = list("PDF", "PNG", "SVG", "TIFF")
                                    ),
                                    downloadButton(outputId = "mdsDownload", lable = "Download Plot")))                              

                         ))
                            

#====================================================#
## HeatMap module ####
#====================================================#
bodyHeat <- tabItem(tabName = "heatmap",
                  h2("HeatMap with Annotation") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "heatmaphelp") ,
                  
                  fluidRow(
                    useShinyjs(),
                    extendShinyjs(text = jscode, functions = c()),

                    box(
                    title = "Input File",
                    closable = FALSE,
                    width = 6,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    fileInput("file1heat", 
                      label = "Upload your data file", 
                      accept = c(".csv",
                                  ".txt",
                                  ".xls",
                                  ".xlsx"),
                      multiple = FALSE),
                    HTML("<hr> <a href='heatmap_matrixTestData.csv' target='_blank'> <i class='fa fa-download'> </i> example data</a>"),
                    DT::dataTableOutput("heatmapdatafile")
                  ),
                  box(
                    title = "Input File",
                    closable = FALSE,
                    width = 6,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    fileInput("file2heat", 
                      label = "Upload your annotation file (datafile)", 
                      accept = c(".csv",
                                  ".txt",
                                  ".xls",
                                  ".xlsx"),
                      multiple = FALSE),
                    HTML("<hr> <a href='heatmap_annotation.csv' target='_blank'> <i class='fa fa-download'> </i> example data</a>"),                     
                    DT::dataTableOutput("heatmapannotationfile")
                  ),
                  box(
                    title = "Figure with Parameters setup",
                    closable = FALSE,
                    collapsed = TRUE,
                    width = 12,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    label = boxLabel(
                      text = "heatmap", 
                      status = "danger",
                      style = "default"
                    ),
                     dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    sidebar = boxSidebar(
                      startOpen = TRUE,
                      width = 34,
                      id = "heatmapParameters",
                      background = "#1D5E6B",
                      icon = shiny::icon("user-gear"),
                      actionButton(
                        "heatmapmainparams",
                        label = "1. HeatMap Main Parameters",
                        width = "400px",
                        icon("paper-plane"), 
                        style="color: #000000; background-color: #fcffa4; border-color: #2e6da4"
                      ),
                      br(),
                      br(),
                      actionButton("clusterheatmap",
                                  "2. HeatMap Cluster Settings",
                                  width = "400px",
                                  icon("paper-plane"), 
                        style="color: #000000; background-color: #fac228; border-color: #2e6da4"),
                      br(),
                      br(),
                      actionButton("legeparams",
                                  "3. Select Legend Parameters",
                                  width = "400px",
                                  icon("paper-plane"), 
                        style="color: #fff; background-color: #f57d15; border-color: #2e6da4"),
                      br(),
                      br(),
                      actionButton("rowannoparams",
                                  "4. Select Row Annotation Parameters",
                                  width = "400px",
                                  icon("paper-plane"), 
                        style="color: #fff; background-color: #d44842; border-color: #2e6da4"),
                      br(),
                      br(),                    
                      actionButton(
                        "coloheatmap",
                        "5. Color settings",
                        width = "400px",
                        icon("paper-plane"), 
                        style="color: #fff; background-color: #9f2a63; border-color: #2e6da4"
                      ),
                      br(),
                      br(),                    
                      actionButton(
                        "fontheatmap",
                        "6. Change Font Family",
                        width = "400px",
                        icon("paper-plane"), 
                        style="color: #fff; background-color: #65156e; border-color: #2e6da4"
                      ),
                      br(),
                      br(),
                      # Font Size settings
                      actionButton(
                        "fontsizeheatmap",
                        "7. Change Font Size",
                        width = "400px",
                        icon("paper-plane"), 
                        style="color: #fff; background-color: #280b53; border-color: #2e6da4"
                      ),
                      br(),
                      br(),
                      # Font Alpha settings
                      actionButton(
                        "fontalphaheatmap",
                        "8. Change Font Color Transparency",
                        width = "400px",
                        icon("paper-plane"), 
                        style="color: #fff; background-color: #000004; border-color: #2e6da4"
                      ),
                      br(),
                      br(),                      
                      awesomeCheckbox(
                        "annotationheatmap",
                        label = "HeatMap Annotation Settings",
                        value = FALSE,
                        status = "success",
                        width = "400px"
                      ),
                      conditionalPanel(
                        condition = "input.annotationheatmap == true",
                        h3("Group-wise (Character/Text-based) Annotations"),
                        h4("First Annotation Category"),
                        selectInput("groupannoval1",
                                    "Select Column for Group-wise Annotation:", ""),
                        actionButton("heatmapannoparams1",
                                  "Select Annotation Parameters"),
                        # character-based annotation-2
                        awesomeCheckbox("groupannosecond",
                                      label = "Choose Another Group-wise Annotation",
                                      value = FALSE, 
                                      status = "success",
                                      width = "400px"),
                        conditionalPanel(
                          condition = "input.groupannosecond == true",
                          selectInput("groupannoval2",
                                    "Select Another Column for Group-wise Annotation:", ""),
                          actionButton("heatmapannoparams2",
                                      "Select Annotation Parameters")
                        ),
                        # number-based annotation-1
                        h3("Numeric/Number-based Annotations"),
                        selectInput("annotationval1", "Select Column for Annotation (Numeric Value):", ""),
                        actionButton("heatmapannoparams3",
                                      "Select Annotation Parameters"),
                        #second annotation
                        awesomeCheckbox("annoagain1",
                                        label = "Want to Add Another Annotation",
                                        value = FALSE,
                                        status = "success",
                                        width = "300px"),
                        conditionalPanel(
                          condition = "input.annoagain1 == true",
                        selectInput("annotationval2", "Select Another Column for Annotation:", ""),
                        actionButton("heatmapannoparams4",
                                      "Select Annotation Parameters")),
                        # third annotation
                        awesomeCheckbox("annoagain2",
                                        label = "Want to Add Another Annotation",
                                        value = FALSE,
                                        status = "success",
                                        width = "300px"),
                        conditionalPanel(
                          condition = "input.annoagain2 == true",
                        selectInput("annotationval3", "Select Another Column for Annotation:", ""),
                        actionButton("heatmapannoparams5",
                                      "Select Annotation Parameters")
                        ),
                        # fourth annotation
                        awesomeCheckbox("annoagain3",
                                        label = "Want to Add Another Annotation",
                                        value = FALSE,
                                        status = "success",
                                        width = "300px"),
                        conditionalPanel(
                          condition = "input.annoagain3 == true",
                            selectInput("annotationval4", "Select Another Column for Annotation:", ""),
                        actionButton("heatmapannoparams6",
                                      "Select Annotation Parameters")
                      )),
                      br(),
                      br(),
                      actionButton("submitheatmap", "Submit")
                  ),
                  plotOutput("heatmapplot", height = "1400px", width = "1000px")
                  ),
                  box(
                      title = "Download Plot",
                      closable = FALSE,
                      collapsed = TRUE,
                      width = 12,
                      status = "navy",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                              radioButtons(inputId = "filetype_heat",
                              label = "Choose file type to download the plot:",
                              inline = TRUE,
                              choices = list("PDF", "PNG", "SVG", "TIFF")),
                            conditionalPanel(
                              condition = "input.filetype_heat == 'PNG'",                            
                            numericInput("heat_pngheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 50,
                                step = 1),
                            numericInput("heat_pngwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 35,
                                step = 1),
                            selectInput("heat_pngunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "cm"
                                ),
                            numericInput("heat_pngresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 150,
                                max = 400,
                                step = 75)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_heat == 'TIFF'",

                            numericInput("heat_tiffheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 10,
                                step = 1),

                            numericInput("heat_tiffwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 12,
                                step = 1),

                            selectInput("heat_tiffunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "inch"
                                ),

                            numericInput("heat_tiffresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 300)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_heat == 'PDF'",

                            numericInput("heat_pdfheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 30,
                                value = 20,
                                step = 1),

                            numericInput("heat_pdfwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 30,
                                value = 14,
                                step = 1)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_heat == 'SVG'",

                            numericInput("heat_svgheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 30,
                                value = 20,
                                step = 1),

                            numericInput("heat_svgwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 30,
                                value = 14,
                                step = 1)
                            ),
                            downloadButton(outputId = "heatDownload", lable = "Download Plot")
                    )
                  )
                  )



#====================================================#
## Venn module ####
#====================================================#
bodyVenn <- tabItem(tabName = "venn",
                    h2("Venn diagrams") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "VennDiagramhelp"),
                    
                    fluidRow(
                      # box(
                      #   title="Usage Instructions", 
                      #   width = 12, 
                      #   status = "primary",
                      #   solidHeader = TRUE,
                      #   collapsible = TRUE,
                      #   collapsed = TRUE,
                      #   closable = FALSE,                        
                      #   div(style="font-size:14px",
                      #       tags$ol("1. Please upload TEXT (.txt), Comma-separated file (.csv) file as input."),  
                      #       tags$ol(HTML("2. Each column represents a set, and each row represents an element (names/gene/SNPs). User can upload a file containing maximum of 6 sets")),
                      #       tags$ol(tags$a(href= "https://asntech.shinyapps.io/intervene/","3. Examples are taken from Intervene shiny app.", target="_blank")),
                      #       tags$ol(HTML("4. In addtion, we provided a small <i>question mark</i> or <i>tip</i> for parameter. Please check those (and the manual, if required) if you are running <b><i>transcriptR</i></b> for the first time.")),
                      #          tags$ol(HTML("5. If you are facing any issue with the pipeline, please contact us via <b><i>support</i></b> or <b><i>google-groups</i></b> or <b><i>email</i></b>."))
                      #   ),
                      #   h4("R packages used"),
                      #   div(style = "font-size:14px",
                      #       tags$ol(tags$a(href= "https://github.com/js229/Vennerable", "1. Vennerable", taget= "_blank")),
                      #       tags$ol(tags$a(href= "https://asntech.shinyapps.io/intervene/", "2. intervene", taget= "_blank"))
                      #       )
                      # ),
                      box( title = "Data Upload & Parameters Setup & settings", width = 4, 
                      status="navy",
                      solidHeader = TRUE,
                              collapsible = TRUE,
                              collapsed = FALSE,
                              closable = FALSE,
                      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ), 
                           tabBox(
                             id = "venntab", height = "100%", width = "100%",
                             
                             tabPanel("Upload",
                                      fileInput(
                                        'file_venn',
                                        label = h4("Upload file"),
                                        accept = c(
                                          'text/csv',
                                          'text/comma-separated-values',
                                          'text/tab-separated-values',
                                          '.csv',
                                          '.tsv'
                                        )
                                      ),
                                      checkboxInput('header_venn', label = 'Header', TRUE),
                                      radioButtons(
                                        'sep_venn',
                                        label = h4('Separator'),                                        
                                        choices = c(
                                          Comma = ',',
                                          Tab = '\t',
                                          Semicolon = ';'
                                        ),
                                        selected = ','
                                      ),
                                      br(),
                                      HTML("<hr> <a href='Whyte_et_al_2013_SEs_genes.csv'> <i class='fa fa-download'> </i> example data</a>")
                             ),
                             tabPanel("Settings",
                                      htmlOutput("venn_sets"),
                                      selectInput(
                                        "venn_type",
                                        label = h4("Venn type"),
                                        choices = list("Classical" = "Classical",
                                                       "Chow-Ruskey" = "ChowRuskey",
                                                       "Edwards" = "AWFE",
                                                       "Squares" = "squares",
                                                       "Battle" = "battle"
                                        ),
                                        selected = "Classical"),
                                      
                                      checkboxInput('doWeights', label = "Weighted", value = TRUE),
                                      checkboxInput('doEuler', label = "Eular", value = FALSE),
                                      
                                      sliderInput(
                                        "venn_lwd",
                                        label = h4("Border line width"),
                                        value = 2.0,
                                        min = 0.0,
                                        max = 10.0,
                                        ticks = TRUE,
                                        step = 0.5),
                                      
                                      selectInput('venn_lty', label=h4("Border line type"),
                                                  choices = list(
                                                    "Solid" = 1,
                                                    "Dashed" = 2,
                                                    "Dotted" = 3,
                                                    "Dot dash" = 4,
                                                    "Long dash" = 5,
                                                    "Two dash" = 6,
                                                    "Blank" = 0
                                                  ),
                                                  selected = "1"
                                      ),
                                      
                                      sliderInput(
                                        "venn_size",
                                        label = h4("Zoom in/out Venn diagram"),
                                        value = 500,
                                        min = 200,
                                        max = 1200,
                                        ticks = FALSE,
                                        step = 10)
                             ),
                             tabPanel("Font & Color",
                                      selectInput('venn_color_type', label= h4("Select color theme"),
                                                  choices = list(
                                                    Set1 = "Set1",                                          
                                                    Custom = "custom"
                                                  ),
                                                  selected = "Set1"
                                      ),
                                      conditionalPanel(condition = "input.venn_color_type=='custom'",
                                                       fluidRow(
                                                         column(2,colourpicker::colourInput("set1_color", label = "Set1", value = "#E41A1C")),
                                                         column(2, colourpicker::colourInput("set2_color", label = "Set2", value = "#377EB8")),
                                                         column(2,colourpicker::colourInput("set3_color", label = "Set3", value = "#4DAF4A")),
                                                         column(2, colourpicker::colourInput("set4_color", label = "Set4", value = "#984EA3")),
                                                         column(2,colourpicker::colourInput("set5_color", label = "Set5", value = "#FF7F00")),
                                                         column(2, colourpicker::colourInput("set6_color", label = "Set6", value = "#FFFF33")))
                                      ),
                                      numericInput(
                                        "venn_labelsize",
                                        label = h4("Label font size"),
                                        value = 15,
                                        min = 1,
                                        max = 50),
                                      
                                      numericInput(
                                        "venn_cex",
                                        label = h4("Number font size"),
                                        value = 1.5,
                                        min = 0.5,
                                        max = 20)
                             )
                           )
                      ),
                      box(
                        title = "Venn Diagram",
                        status="navy", 
                        width = 8,
                        solidHeader = TRUE,
                              collapsible = TRUE,
                              collapsed = FALSE,
                              closable = FALSE,
                        dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ), 
                        tabBox(
                          # The id lets us use input$tabset1 on the server to find the current tab
                          id = "tabset1", height = "100%", width = "100%",
                          tabPanel("Venn diagram",
                                   
                                   #htmlOutput('plot_text_p'),
                                   plotOutput("vennPlot", width = "100%", height = "100%"),
                                   box(title = "Download Plot",
                                     width = NULL,  status = "navy",
                                     dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                                     radioButtons(
                                       inputId = "filetype_venn",
                                       label = "Choose file type to download:",
                                       inline = TRUE,
                                       choices = list("PDF", "PNG","SVG","TIFF")),
                                     
                                     downloadButton(outputId = "VennDown", label = "Download Plot")
                                   )
                          )
                        )
                      )
                    )
)

#====================================================#
## UpSet module ####
#====================================================#
bodyUpSet <- tabItem(tabName = "upset", value="upset_plot",
                     h2("UpSet plots") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "UpSetplothelp"),
                     fluidRow(
                      #  box(
                      #    title="Usage Instructions", 
                      #    width =12, 
                      #    status = "primary",
                      #    solidHeader = TRUE,
                      #    collapsible = TRUE,
                      #    collapsed = TRUE,
                      #    closable = FALSE,                         
                      #    div(style="font-size:14px",
                      #        tags$ol("1. Please upload TEXT (.txt), Comma-separated file (.csv) file as input."),  
                      #        tags$ol(HTML("2. Each column represents a set, and each row represents an element (names/gene/SNPs).")),
                      #        tags$ol(tags$a(href= "https://asntech.shinyapps.io/intervene/","3. Examples are taken from Intervene shiny app.", target="_blank")),
                      #        tags$ol(HTML("4. In addtion, we provided a small <i>question mark</i> or <i>tip</i> for each parameter. Please check those (and the manual, if required) if you are running <b><i>transcriptR</i></b> for the first time.")),
                      #        tags$ol(HTML("5. If you are facing any issue with the pipeline, please contact us via <b><i>support</i></b> or <b><i>google-groups</i></b> or <b><i>email</i></b>."))
                      #    ),
                      #    h4("R packages used"),
                      #    div(style = "font-size:14px",
                      #        tags$ol(tags$a(href= "https://cran.r-project.org/web/packages/UpSetR/index.html", "1. UpSetR", taget= "_blank")),
                      #        tags$ol(tags$a(href= "https://asntech.shinyapps.io/intervene/", "2. intervene", taget= "_blank"))
                      #        )
                      #  ),
                       box(title = "Data Upload & Parameters Setup & settings",
                         status="navy",
                         width = 4,
                         solidHeader = TRUE,
                              collapsible = TRUE,
                              collapsed = FALSE,
                              closable = FALSE, 
                         dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),     
                         tabBox(
                           # The id lets us use input$tabset1 on the server to find the current tab
                           id = "upset_plot", height = "100%", width = "100%",
                           
                           tabPanel("Upload",
                                    radioButtons('upset_input_type',
                                                 label = h4('Input type '),
                                                 choices = c(
                                                   "List of Genes/SNPs" = 'list',
                                                   "Binary data (0 & 1)" = 'binary'
                                                 ),
                                                 selected = 'list'
                                    ),
                                    fileInput(
                                      'file1',
                                      label = h4("Upload file", tipify(actionButton("button11", icon("question")), 
                                      title = "Please check the file format as in the example file to prepare the right format to upload.",
                                      placement = "right")),
                                      accept = c(
                                        'text/csv',
                                        'text/comma-separated-values',
                                        'text/tab-separated-values',
                                        '.csv',
                                        '.tsv'
                                      )
                                    ),
                                    checkboxInput('header', label = 'Header', TRUE),
                                    radioButtons('sep',
                                                 label =h4('Separator'),
                                                 choices = c(
                                                   Comma = ',',
                                                   Semicolon = ';',
                                                   Tab = '\t'
                                                 ),
                                                 selected = ','
                                    ),                                    
                                    textAreaInput('upset_comb', label = "OR enter set combinations/expression", rows = 4, placeholder = "Enter combinations of sets to plot"),
                                    p("For example: A=3, B=3, C=2, A&B=1, A&C=2, B&C=1 ,A&B&C=1"),                                    
                                    
                                    HTML("<hr> <a href='Whyte_et_al_2013_SEs_genes.csv'> <i class='fa fa-download'> </i> List example data</a> | "),
                                    HTML("<a href='mutations_glioblastoma_TCGA.csv'> <i class='fa fa-download'> </i> Binary example data</a>")
                           ),
                           
                           tabPanel("Settings",
                                    #add content
                                    htmlOutput("sets"),
                                    numericInput(
                                      "nintersections",
                                      label = h4("Number of intersections to show"),
                                      value = 30,
                                      min = 1,
                                      max = 60),
                                    
                                    selectInput(
                                      "order",
                                      label = h4("Order intersections by"),
                                      choices = list("Degree" = "degree",
                                                     "Frequency" = "freq"),
                                      selected = "freq"),
                                    selectInput(
                                      "decreasing",
                                      label = h4("Increasing/Decreasing"),
                                      choices = list("Increasing" = "inc",
                                                     "Decreasing" = "dec"),
                                      selected = "dec"),
                                    
                                    selectInput(
                                      "scale.intersections",
                                      label = h4("Scale intersections"),
                                      choices = list("Original" = "identity",
                                                     "Log10" = "log10",
                                                     "Log2" = "log2"),
                                      selected = "identity"),
                                    
                                    selectInput(
                                      "scale.sets",
                                      label = h4("Scale sets"),
                                      choices = list("Original" = "identity",
                                                     "Log10" = "log10",
                                                     "Log2" = "log2"),
                                      selected = "identity"),
                                    
                                    
                                    sliderInput(
                                      "upset_width",
                                      label = h4("Plot width"),
                                      value = 650,
                                      min = 400,
                                      max = 1200,
                                      ticks = FALSE,
                                      step = 10),
                                    
                                    sliderInput(
                                      "upset_height",
                                      label = h4("Plot height"),
                                      value = 400,
                                      min = 300,
                                      max = 1000,
                                      ticks = FALSE,
                                      step = 10),
                                    
                                    sliderInput(
                                      "mbratio",
                                      label = h4("Bar matrix ratio"),
                                      value = 0.30,
                                      min = 0.20,
                                      max = 0.80,
                                      ticks = FALSE,
                                      step = 0.01),
                                    
                                    checkboxInput('show_numbers', label = h4("Show numbers on bars"), value = TRUE),
                                    
                                    sliderInput(
                                      "angle",
                                      label = h4("Angle of number on the bar"),
                                      min = -90,
                                      max = 90,
                                      value = 0,
                                      step = 1,
                                      ticks = F),
                                    
                                    checkboxInput('empty', label = h4("Show empty intersections"), value = FALSE),
                                    checkboxInput('keep.order', label = h4("Keep set order"), value = FALSE),
                                    
                                    numericInput(
                                      "pointsize",
                                      label = h4("Connecting point size"),
                                      value = 2.5,
                                      min = 1,
                                      max = 10),
                                    
                                    numericInput(
                                      "linesize",
                                      label = h4("Connecting line size"),
                                      value = 0.8,
                                      min = 0.5,
                                      max = 10)
                           ),
                           tabPanel("Font & Colors",
                                    
                                    colourpicker::colourInput("mbcolor", h4("Select main bar colour"), "#ea5d4e"),
                                    colourpicker::colourInput("sbcolor", h4("Select set bar colour"), "#317eab"),
                                    
                                    numericInput(
                                      "intersection_title_scale",
                                      label = h4("Font size of intersection size label"),
                                      value = 1.8,
                                      min = 0.5,
                                      max = 100),
                                    numericInput(
                                      "set_title_scale",
                                      label = h4("Set size label font"),
                                      value = 1.8,
                                      min = 0.5,
                                      max = 100),
                                    numericInput(
                                      "intersection_ticks_scale",
                                      label = h4("Intersection size ticks font"),
                                      value = 1.2,
                                      min = 0.5,
                                      max = 100),
                                    numericInput(
                                      "set_ticks_scale",
                                      label = h4("Set Size ticks font"),
                                      value = 1.5,
                                      min = 0.5,
                                      max = 100),
                                    numericInput(
                                      "intersection_size_numbers_scale",
                                      label = h4("Intersection Size Numbers font size"),
                                      value = 1.2,
                                      min = 0.5,
                                      max = 100),
                                    numericInput(
                                      "names_scale",
                                      label = h4("Set Names font size"),
                                      value = 1.5,
                                      min = 0.5,
                                      max = 100
                                    )
                           )
                         )
                       ),
                       box(
                        title = "UpSet Plot",
                         status="navy", 
                         width = 8, 
                         height = "100%",
                         solidHeader = TRUE,
                              collapsible = TRUE,
                              collapsed = FALSE,
                              closable = FALSE,
                         dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),      
                         tabBox(
                           id = "upsetplottab", height = "100%", width = "100%",
                           tabPanel("UpSet Plot",
                                    htmlOutput('plot_text'),
                                    plotOutput('plot1', width = "100%", height = "100%"),
                                    box(title = "Download Plot",
                                      width = NULL,  status = "navy",
                                      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                                      radioButtons(
                                        inputId = "filetype",
                                        label = "Choose file type to download:",
                                        inline = TRUE,
                                        choices = list("PDF", "PNG","SVG","TIFF")),
                                      
                                      downloadButton(outputId = "UpSetDown", label = "Download Plot")
                                    )
                           )
                         )
                       )
                     )
)

#========================================#
## Volcano plot
#========================================#

bodyVol <- tabItem(tabName = "volcano",
                h2("Volcano Plot") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "volcanohelp"),

                fluidRow(
                  useShinyjs(),
                  extendShinyjs(text = jscode, functions = c()),
                  box(
                    title = "Input File",
                    closable = FALSE,
                    width = 6,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    fileInput("file1_vol", 
                      label = "Upload your data file" , 
                      accept = c(".csv",
                                  ".txt",
                                  ".xls",
                                  ".xlsx"),
                      multiple = FALSE),
                    tags$table(width="100%",
                    
                    tags$tr(
                      tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select gene name column: "),
                      tags$td(width = "50%",          
                        selectInput("geneName_volc", label = NULL, "")
                      ))),
                    tags$tr(
                      tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select logFC column: "),
                      tags$td(width = "50%",          
                        selectInput("logfc_volc", label = NULL, "")
                      ))),
                    tags$tr(
                      tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select adjusted/ p-value column: "),
                      tags$td(width = "50%",          
                        selectInput("pvalue_volc", label = NULL, "")
                      )))
                      ),
                    HTML("example data: <a href='volcanoplot_testdata.xlsx'> test data</a>")                    
                  ),
                  box(
                    title = "Input File Viewer",
                    closable = FALSE,
                    width = 6,
                    status = "info",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    label = boxLabel(
                      text = "Volv", 
                      status = "primary",
                      style = "default"
                    ),
                    dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),                    
                    DT::dataTableOutput("volcanofile1")
                  ),
                  box(
                    title = "Static Volcano Plot with Parameters setup",
                    closable = FALSE,
                    collapsed = TRUE,
                    width = 6,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    label = boxLabel(
                      text = "Vols", 
                      status = "danger",
                      style = "default"
                    ),
                     dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    sidebar = boxSidebar(
                      startOpen = TRUE,
                      width = 34,
                      id = "volcanoParameters",
                      background = "#1D5E6B",
                      icon = shiny::icon("user-gear"),
                        numericInput("niVol",
                        label = "Select number of genes to show",
                        value = 20,
                        min = 10,
                        max = 1000,
                        step = 10),
                        numericInput("wordsize2", 
                                    label = "Change the genes Font size",
                                    value = 2,
                                    max = 5,
                                    step = 0.5),
                        numericInput("pointsize2", 
                                    label = "Change the points size",
                                    value = 2,
                                    #max = 5,
                                    step = 0.5),
                        numericInput("fc",
                                    label = "Select cut-off for log2FC",
                                    value = 1,
                                    min = 0.1,
                                    #max = 1,
                                    step = 0.1),
                        numericRangeInput("xlim_volc",
                                  h4("Choose X-axis range"),
                                  value = c(-10,10)),
                        numericRangeInput("ylim_volc",
                                  h4("Choose Y-axis range"),
                                  value = c(0,10)),
                        awesomeCheckbox("volcolchange",
                                        label = "Change Colors",
                                        value = FALSE,
                                        status = "success",
                                        width = "400px"
                                        ),
                        conditionalPanel(
                          condition = "input.volcolchange == true",
                            colourInput("highcol", "Select color for downregulated", "dodgerblue3"),
                            colourInput("lowcol", "Select color for upregulated", "firebrick3"),
                            colourInput("nonsigcol", "Select color for non-significant", "gray")
                            ),
                        awesomeCheckbox("volfontparams",
                                        label = "Change Font Size, Family and Color",
                                        value = FALSE,
                                        status = "success",
                                        width = "400px"
                                        ),
                        conditionalPanel(
                          condition = "input.volfontparams == true",
                        numericInput("wordsize", 
                                    label = "Change the Font size",
                                    value = 14),
                        colourInput("volfontcol", "Select for the Text", "gray"),
                        selectInput("fontfamily",
                                    label = "Change the Font Family",
                                    choices = c(
                                      "Arial" = "sans",
                                      "Times New Roman" = "Times",
                                      "Courier" = "mono"
                                    ),
                                    multiple = FALSE,
                                    selected = "Times")),
                        awesomeCheckbox("vollinesparams",
                                        label = "Adjust lines",
                                        value = FALSE,
                                        status = "success",
                                        width = "400px"
                                        ),
                        conditionalPanel(
                          condition = "input.vollinesparams == true",
                        colourInput("vlinecol", "Select color for log2FC cut-off line", "gray"),
                        colourInput("hlinecol", "Select color for p-value cut-off line", "gray"),
                        selectInput("linetype",
                                    label = "Change the line type",
                                    choices = c(
                                      "Blank" = "blank", 
                                      "Solid" = "solid", 
                                      "Dashed" = "dashed", 
                                      "Dotted" = "dotted", 
                                      "Dotdash" = "dotdash", 
                                      "Longdash" = "longdash", 
                                      "Twodash" = "twodash"),
                                    multiple = FALSE,
                                    selected = "dotted")),
                        awesomeCheckbox("volaxislabparams",
                                        label = "Lables and Title",
                                        value = FALSE,
                                        status = "success",
                                        width = "400px"
                                        ),
                        conditionalPanel(
                          condition = "input.volaxislabparams == true",
                        textInput("xaxislab_volc", "X-axis label", ""),
                        textInput("yaxislab_volc", "Y-axis label", ""),
                        textInput("title_volc", "Title", "") 
                    ), actionButton("submitvol", "Submit")),
                    plotOutput("volcanoplot", height = "800px", width = "100%")
                  ),
                    
                  box(
                    title = "Interactive Volcano Plot with Parameters setup",
                    closable = FALSE,
                    collapsed = TRUE,
                    width = 6,
                    status = "navy",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    label = boxLabel(
                      text = "Voli", 
                      status = "danger",
                      style = "default"
                    ),
                     dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"),
                          # boxDropdownItem("Parameters Details" %>%
                          # shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp")),
                          # dropdownDivider(),
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                    sidebar = boxSidebar(
                      startOpen = TRUE,
                      width = 34,
                      id = "volParams",
                      background = "#1D5E6B",
                      icon = shiny::icon("user-gear"),
                      numericInput("volp", 
                                  h4("Choose BH-corrected p-value cut-off"),
                                  value = 0.01),
                      numericInput("volfc",
                                  h4("Choose logFC cut-off"),
                                  value = 0.3, 
                                  step = 0.1),                     
                      actionButton("runVol", "Run Analysis"),
                         HTML("<hr> <a href='volcanoplot_example_data.txt' target='_blank'> <i class='fa fa-download'> </i> example data</a>"),
                          tipify(actionButton("button7d", icon("question")),
                                        title = "To download the example file, right click and save the file.",
                                        placement = "right"
                                        )
                    ),
                    plotlyOutput("volplot", height = "800px", width = "100%"),
                    downloadButton(outputId = "downloadvolcano", lable = "Download Plot")
                  ),
                  box(
                      title = "Download Static Plot",
                      closable = FALSE,
                      collapsed = TRUE,
                      width = 6,
                      status = "navy",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
                              radioButtons(inputId = "filetype_volc",
                              label = "Choose file type to download the plot:",
                              inline = TRUE,
                              choices = list("PDF", "PNG", "SVG", "TIFF")),
                            conditionalPanel(
                              condition = "input.filetype_volc == 'PNG'",

                            numericInput("volc_pngheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 10,
                                step = 1),

                            numericInput("volc_pngwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 12,
                                step = 1),

                            selectInput("volc_pngunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "inch"
                                ),

                            numericInput("volc_pngresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 150,
                                max = 400,
                                step = 75)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_volc == 'TIFF'",

                            numericInput("volc_tiffheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                value = 10,
                                step = 1),

                            numericInput("volc_tiffwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                value = 12,
                                step = 1),

                            selectInput("volc_tiffunits",
                                label = "Select units for figure resolution",
                                choices = c(
                                  "mm" = "mm",
                                  "cm" = "cm",
                                  "inch" = "in",
                                  "pixel" = "px"
                                ),
                                multiple = FALSE,
                                selected = "inch"
                                ),

                            numericInput("volc_tiffresol", 
                                label = "Select the resolution of the plot",
                                min = 75,
                                value = 300)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_volc == 'PDF'",

                            numericInput("volc_pdfheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 20,
                                value = 10,
                                step = 1),

                            numericInput("volc_pdfwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 20,
                                value = 12,
                                step = 1)
                            ),
                            conditionalPanel(
                              condition = "input.filetype_volc == 'SVG'",

                            numericInput("volc_svgheight", 
                                label = "Select the height of the plot",
                                min = 4,
                                max = 20,
                                value = 10,
                                step = 1),

                            numericInput("volc_svgwidth", 
                                label = "Select the width of the plot",
                                min = 4,
                                max = 20,
                                value = 12,
                                step = 1)
                            ),
                            downloadButton(outputId = "volcDownload", lable = "Download Plot")
                    )
                )
                )

#========================================#
## Gene Ontology (GO) analysis
#========================================#
bodyGo <- tabItem(
  tabName = "go",
  h2("Gene Ontology (GO) analysis") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "geneontologyhelp"),
  fluidRow(
    
    box(
      title = "Data Upload & Parameters Setup",
      width = 4,
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      closable = FALSE,
      label = boxLabel(
                      text = "params", 
                      status = "danger",
                      style = "default"
                    ),
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#shiny-tab-manual", icon = icon("file-pdf"))
                          ),
      # tabBox(
      #   id = "go_plot", height = "100%", width = "100%",
      tags$table(width = "100%",
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select organism")),
          tags$td(width = "50%", 
          selectInput(
            inputId = "goOrg",
            label = NULL,
            choices = list(
              "Human"="org.Hs.eg.db",
              "Mouse"="org.Mm.eg.db",
              "Rat"="org.Rn.eg.db",
              "Fly"="org.Dm.eg.db",
              "Arabidopsis"="org.At.tair.db",
              "Zebrafish"="org.Dr.eg.db",
              "Yeast"="org.Sc.sgd.db",
              "Worm"="org.Ce.eg.db",
              "Bovine"="org.Bt.eg.db",
              "Chicken"="org.Gg.eg.db",
              "Pig"="org.Ss.eg.db",
              "Rhesus"="org.Mmu.eg.db",
              "Canine"="org.Cf.eg.db",
              "E coli strain K12"="org.EcK12.eg.db",
              "Xenopus"="org.Xl.eg.db",
              "Chimp"="org.Pt.eg.db",
              "Anopheles"="org.Ag.eg.db",
              "E coli strain Sakai"="org.EcSakai.eg.db",
              "Myxococcus xanthus DK 1622"="org.Mxanthus.db"
            ),
          selected = "org.Hs.eg.db",
          multiple = FALSE
        ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select analysis type")),
          tags$td(width = "50%", 
          selectInput(
            inputId = "goType",
            label = NULL,
            choices = list(
              "Over representation" = 1,
              "GSEA" = 2
            ),
          selected = 1,
          multiple = FALSE
        ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select ontology class")),
          tags$td(width = "50%",
          selectInput("ontology",
            label = NULL,
            choices = list(
              "Cellular component" = "CC",
              "Molecular function" = "MF",
              "Biological process" = "BP"
            ),
            selected = "BP",
            multiple = FALSE
        ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select number of GO classes")),
          tags$td(width = "50%",
          numericInput("gonum",
            label = NULL,
            value = 20
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select p-value for correction")),
          tags$td(width = "50%",
          numericInput("numpval",
            label = NULL,
            value = 0.05
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select p-value adjustment method")),
          tags$td(width = "50%",          
          selectInput("pvaladjgo", 
          label = NULL,
            choices = list(
              "Benjamini-Hochberg" = "BH",
              "Benjamini-Yeketuli" = "BY",
              "Bonferroni" = "bonferroni",
              "Holm" = "holm",
              "Hommel" = "hommel",
              "Hochberg" = "hochberg",
              "FDR" = "fdr",
              "None" = "none"
            ),
            selected = "BH",
            multiple = FALSE
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select q-value")),
          tags$td(width = "50%",
          numericInput("numqval",
            label = NULL,
            value = 0.05
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", h6("Upload data file"))),
          tags$td(width = "50%",
          fileInput("file6",
            label = NULL,
            accept = c(
              ".csv",
              ".txt",
              ".xls",
              ".xlsx"
            )
          )))
      ),
      # actionButton("goDemo", "Use Demo Data"),
      actionButton("goRun", "Run Analysis",
      width = "75%",
        icon("paper-plane"),
        style = "color: #fff; background-color: #374806; border-color: #2e6da4; position:relative; left:calc(15%);"),
      HTML("<hr> <a href='pathway_example_data.txt' target='_blank'> <i class='fa fa-download'> </i> example data</a>"),
      tipify(actionButton("button7c", icon("question")),
        title = "To download the example file, right click and save the file.",
        placement = "right"
      )
    ),
    box(
      title = "Analysis Result",
      status = "navy",
      width = 8,
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      closable = FALSE,
      label = boxLabel(
                      text = "GO", 
                      status = "danger",
                      style = "default"
                    ),
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),
      tabBox(
        id = "go_tab", height = "100%", width = "100%",
        tabPanel(
          "Gene Ontology Plot",
          plotlyOutput("goplot", width = "1200px", height = "600px"),
          # conditionalPanel(
          # condition = "input.goRun == true",
          # box(
          downloadButton(outputId = "downloadGO", lable = "Download Gene Ontology Plot")
          # )
          # )
        ),
        tabPanel(
          "Gene Ontology data",
          div(style = "overflow-x: scroll", DT::dataTableOutput("tableGO", width = "80%")),
          #  br(),
          #  downloadButton('gotable', 'Download Gene Ontology data')
        )
      )
    )
  )
)

#========================================#
## Pathway enrichment analysis
#========================================#
bodyPathway <- tabItem(
  tabName = "pathway",
  h2("Pathway Enrichment Analysis") %>%
                          shinyhelper::helper(type = 'markdown', colour = "#07750c", icon = "info-circle", content = "PathwayEnrichmentAnalysishelp"),
  fluidRow(
    box(
      tags$style(
        HTML(".shiny-notification {
                            position:fixed;
                            top: calc(50%);
                            left: calc(50%);
                            background-color: black;
                            color: white;
                            font-size: 1.5em;
                            font-family: papyrus
                          }
                          ")
      ),
      tags$style(
        "#shiny-notification-progress_notif {
                            position:fixed;
                            top: calc(30%);
                            left: calc(50%);
                            background-color: black;
                            color: white;
                            font-size: 1.5em;
                            font-family: papyrus
                          }
                          "
      ),
      title = "Data Upload & Parameters Setup",
      width = 6,
      height = "500px",
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      closable = FALSE,
      label = boxLabel(
            text = "params", # Change to number when finished all package uploads
            status = "danger",
            style = "default"
          ),
          dropdownMenu = boxDropdown(
            icon = shiny::icon("file-pdf"),
            boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
          ),
      tags$table(width = "100%",
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select Database")),
          tags$td(width = "50%", selectInput("pathDB",
          label = NULL,
          #h4("Choose pathway database"),
          choices = list(
            "KEGGdb" = 2, #2
            "ReactomeDB" = 1, #1            
            "Wikipathways" = 3
          ),
          selected = 2,
          multiple = FALSE
        ))
        ),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select Organism")),
          tags$td(width = "50%", 
          conditionalPanel(
            condition = "input.pathDB == '1' ",
            selectInput("pathOrgReac",
              label = NULL,              
              choices = list(
                "Human" = "human",
                "Rat" = "rat",
                "Mouse" = "mouse",
                "C. elegans" = "celegans",
                "Yeast" = "yeast",
                "Zebrafish" = "zebrafish",
                "Fly" = "fly"
              ),
              selected = "human",
              multiple = FALSE
            )
          ),
          conditionalPanel(
            condition = "input.pathDB == '3' ",
            selectInput("pathOrgWiki",
              label = NULL,              
              choices = list(
                "Homo sapiens",
                "Mus musculus",
                "Rattus norvegicus",          
                "Arabidopsis thaliana",
                "Bos taurus",    
                "Caenorhabditis elegans",
                "Canis familiaris",  
                "Danio rerio",       
                "Drosophila melanogaster",
                #"Equus caballus",
                "Gallus gallus",          
                "Pan troglodytes",
                #"Populus trichocarpa",
                "Saccharomyces cerevisiae",
                #"Solanum lycopersicum",
                "Sus scrofa"
              ),
              selected = "Homo sapiens",
              multiple = FALSE
            )
          ),
          conditionalPanel(
            condition = "input.pathDB == '2' ",
            selectInput("pathOrgKegg",
              label = NULL,              
              choices = list(
                "Homo sapiens (human)" = "hsa",
                "Mus musculus (house mouse)" = "mmu",
                "Rattus norvegicus (rat)" = "rno",

                "Acanthisitta chloris (rifleman)" = "achl",
                "Acanthopagrus latus (yellowfin seabream)" = "alat",
                "Accipiter gentilis (Northern goshawk)" = "agen",
                "Acinonyx jubatus (cheetah)" = "aju",
                "Acipenser ruthenus (sterlet)" = "arut",
                "Ailuropoda melanoleuca (giant panda)" = "aml",
                "Alligator mississippiensis (American alligator)" = "amj",
                "Alligator sinensis (Chinese alligator)" = "asn",
                "Amphiprion ocellaris (clown anemonefish)" = "aoce",
                "Anas platyrhynchos (mallard)" = "apla",
                "Anguilla anguilla (European eel)" = "aang",
                "Anolis carolinensis (green anole)" = "acs",
                "Anoplopoma fimbria (sablefish)" = "afb",
                "Anser cygnoides domesticus (swan goose)" = "acyg",
                "Antechinus flavipes (yellow-footed antechinus)" = "afz",
                "Antrostomus carolinensis (chuck-will's-widow)" = "acar",
                "Apaloderma vittatum (bar-tailed trogon)" = "avit",
                "Aptenodytes forsteri (emperor penguin)" = "afor",
                "Apteryx mantelli mantelli (North Island brown kiwi)" = "aam",
                "Apteryx rowi (Okarito brown kiwi)" = "arow",
                "Aquila chrysaetos chrysaetos (golden eagle)" = "achc",
                "Artibeus jamaicensis (Jamaican fruit-eating bat)" = "ajm",
                "Arvicanthis niloticus (African grass rat)" = "anu",
                "Arvicola amphibius (Eurasian water vole)" = "aamp",
                "Astyanax mexicanus (Mexican tetra)" = "amex",
                "Athene cunicularia (burrowing owl)" = "acun",
                "Austrofundulus limnaeus (annual killifish)" = "alim",
                "Aythya fuligula (tufted duck)" = "aful",
                "Balaenoptera acutorostrata scammoni (minke whale)" = "bacu",
                "Balearica regulorum gibbericeps (East African grey crowned-crane)" = "breg",
                "Betta splendens (Siamese fighting fish)" = "bspl",
                "Bison bison bison (American bison)" = "bbis",
                "Boleophthalmus pectinirostris (great blue-spotted mudskipper)" = "bpec",
                "Bos indicus (zebu cattle)" = "biu",
                "Bos mutus (wild yak)" = "bom",
                "Bos taurus (cow)" = "bta",
                "Bubalus bubalis (water buffalo)" = "bbub",
                "Buceros rhinoceros silvestris (Rhinoceros hornbill)" = "brhi",
                "Budorcas taxicolor (takin)" = "btax",
                "Bufo bufo (common toad)" = "bbuf",
                "Bufo gargarizans (Asiatic toad)" = "bgar",
                "Callithrix jacchus (white-tufted-ear marmoset)" = "cjc",
                "Callorhinchus milii (elephant shark)" = "cmk",
                "Camelus bactrianus (Bactrian camel)" = "cbai",
                "Camelus dromedarius (Arabian camel)" = "cdk",
                "Camelus ferus (Wild Bactrian camel)" = "cfr",
                "Canis lupus dingo (dingo)" = "clud",
                "Canis lupus familiaris (dog)" = "cfa",
                "Capra hircus (goat)" = "chx",
                "Carassius auratus (goldfish)" = "caua",
                "Carassius gibelio (silver crucian carp)" = "cgib",
                "Cariama cristata (Red-legged seriema)" = "ccri",
                "Carlito syrichta (Philippine tarsier)" = "csyr",
                "Castor canadensis (American beaver)" = "ccan",
                "Cavia porcellus (domestic guinea pig)" = "cpoc",
                "Cebus imitator (Panamanian white-faced capuchin)" = "cimi",
                "Cercocebus atys (sooty mangabey)" = "caty",
                "Cervus canadensis (wapiti)" = "ccad",
                "Chaetura pelagica (chimney swift)" = "cpea",
                "Charadrius vociferus (killdeer)" = "cvf",
                "Cheilinus undulatus (humphead wrasse)" = "cud",
                "Chelonia mydas (green sea turtle)" = "cmy",
                "Chelonoidis abingdonii (Abingdon island giant tortoise)" = "cabi",
                "Chiloscyllium plagiosum (whitespotted bambooshark)" = "cpla",
                "Chlamydotis macqueenii (Macqueen's bustard)" = "cmac",
                "Chlorocebus sabaeus (green monkey)" = "csab",
                "Chrysemys picta (western painted turtle)" = "cpic",
                "Clupea harengus (Atlantic herring)" = "char",
                "Colius striatus (speckled mousebird)" = "csti",
                "Colobus angolensis palliatus (Angola colobus)" = "cang",
                "Colossoma macropomum (tambaqui)" = "cmao",
                "Columba livia (rock pigeon)" = "clv",
                "Coregonus clupeaformis (lake whitefish)" = "cclu",
                "Corvus brachyrhynchos (American crow)" = "cbrc",
                "Corvus cornix (hooded crow)" = "ccw",
                "Cottoperca gobio" = "cgob",
                "Coturnix japonica (Japanese quail)" = "cjo",
                "Cricetulus griseus (Chinese hamster)" = "cge",
                "Crocodylus porosus (Australian saltwater crocodile)" = "cpoo",
                "Crotalus tigris (Tiger rattlesnake)" = "ctig",
                "Ctenopharyngodon idella (grass carp)" = "cide",
                "Cuculus canorus (common cuckoo)" = "cuca",
                "Cyanistes caeruleus (blue tit)" = "ccae",
                "Cyclopterus lumpus (lumpfish)" = "clum",
                "Cynoglossus semilaevis (tongue sole)" = "csem",
                "Cyprinodon tularosa (White Sands pupfish)" = "ctul",
                "Cyprinodon variegatus (sheepshead minnow)" = "cvg",
                "Cyprinus carpio (common carp)" = "ccar",
                "Danio rerio (zebrafish)" = "dre",
                "Dasypus novemcinctus (nine-banded armadillo)" = "dnm",
                "Delphinapterus leucas (beluga whale)" = "dle",
                "Dermochelys coriacea (leatherback sea turtle)" = "dcc",
                "Desmodus rotundus (common vampire bat)" = "dro",
                "Dipodomys ordii (Ord's kangaroo rat)" = "dord",
                "Dipodomys spectabilis (banner-tailed kangaroo rat)" = "dsp",
                "Dromaius novaehollandiae (emu)" = "dne",
                "Dryobates pubescens (Downy woodpecker)" = "dpub",
                "Echinops telfairi (small Madagascar hedgehog)" = "etf",
                "Egretta garzetta (little egret)" = "egz",
                "Electrophorus electricus (electric eel)" = "eee",
                "Empidonax traillii (willow flycatcher)" = "etl",
                "Enhydra lutris kenyoni (northern sea otter)" = "elk",
                "Epinephelus fuscoguttatus (brown-marbled grouper)" = "efo",
                "Epinephelus lanceolatus (giant grouper)" = "ely",
                "Eptesicus fuscus (big brown bat)" = "efus",
                "Equus asinus (ass)" = "eai",
                "Equus caballus (horse)" = "ecb",
                "Equus przewalskii (Przewalski's horse)" = "epz",
                "Esox lucius (northern pike)" = "els",
                "Etheostoma cragini (Arkansas darter)" = "ecra",
                "Etheostoma spectabile (orangethroat darter)" = "esp",
                "Eublepharis macularius (Leopard gecko)" = "emc",
                "Eumetopias jubatus (Steller sea lion)" = "eju",
                "Eurypyga helias (sunbittern)" = "ehs",
                "Falco cherrug (Saker falcon)" = "fch",
                "Falco peregrinus (peregrine falcon)" = "fpg",
                "Felis catus (domestic cat)" = "fca",
                "Ficedula albicollis (collared flycatcher)" = "fab",
                "Fulmarus glacialis (Northern fulmar)" = "fga",
                "Galeopterus variegatus (Sunda flying lemur)" = "gvr",
                "Gallus gallus (chicken)" = "gga",
                "Gambusia affinis (western mosquitofish)" = "gaf",
                "Gasterosteus aculeatus (three-spined stickleback)" = "gat",
                "Gavia stellata (red-throated loon)" = "gste",
                "Gavialis gangeticus (Gharial)" = "ggn",
                "Gekko japonicus (Schlegel's Japanese gecko)" = "gja",
                "Geospiza fortis (medium ground-finch)" = "gfr",
                "Geotrypetes seraphini (Gaboon caecilian)" = "gsh",
                "Girardinichthys multiradiatus (darkedged splitfin)" = "gmu",
                "Gorilla gorilla gorilla (western lowland gorilla)" = "ggo",
                "Gracilinanus agilis (agile gracile opossum)" = "gas",
                "Gymnogyps californianus (California condor)" = "gcl",
                "Haliaeetus albicilla (white-tailed eagle)" = "hald",
                "Haliaeetus leucocephalus (bald eagle)" = "hle",
                "Hemicordylus capensis (graceful crag lizard)" = "hcg",
                "Heterocephalus glaber (naked mole-rat)" = "hgl",
                "Hippocampus comes (tiger tail seahorse)" = "hcq",
                "Hippoglossus hippoglossus (Atlantic halibut)" = "hhip",
                "Hippoglossus stenolepis (Pacific halibut)" = "hsp",
                "Hipposideros armiger (great roundleaf bat)" = "hai",
                "Hyaena hyaena (striped hyena)" = "hhv",
                "Hylobates moloch (silvery gibbon)" = "hmh",
                "Ictalurus punctatus (channel catfish)" = "ipu",
                "Kryptolebias marmoratus (mangrove rivulus)" = "kmr",
                "Labeo rohita (rohu)" = "lroh",
                "Lagopus muta (rock ptarmigan)" = "lmut",
                "Larimichthys crocea (large yellow croaker)" = "lco",
                "Lates calcarifer (barramundi perch)" = "lcf",
                "Latimeria chalumnae (coelacanth)" = "lcm",
                "Lemur catta (Ring-tailed lemur)" = "lcat",
                "Lepisosteus oculatus (spotted gar)" = "loc",
                "Leptonychotes weddellii (Weddell seal)" = "lww",
                "Leptosomus discolor (cuckoo roller)" = "ldi",
                "Leucoraja erinacea (little skate)" = "leri",
                "Lipotes vexillifer (Yangtze River dolphin)" = "lve",
                "Lonchura striata domestica (Bengalese finch)" = "lsr",
                "Loxodonta africana (African savanna elephant)" = "lav",
                "Lutra lutra (Eurasian river otter)" = "llv",
                "Lynx rufus (bobcat)" = "lruf",
                "Macaca fascicularis (crab-eating macaque)" = "mcf",
                "Macaca mulatta (rhesus monkey)" = "mcc",
                "Macaca nemestrina (pig-tailed macaque)" = "mni",
                "Macaca thibetana thibetana (Pere David's macaque)" = "mthb",
                "Malurus melanocephalus (red-backed fairy wren)" = "mmea",
                "Mandrillus leucophaeus (drill)" = "mleu",
                "Manis javanica (Malayan pangolin)" = "mjv",
                "Mastomys coucha (southern multimammate mouse)" = "mcoc",
                "Mauremys reevesii (Reeves's turtle)" = "mrv",
                "Maylandia zebra (zebra mbuna)" = "mze",
                "Megalobrama amblycephala (Wuchang bream)" = "mamb",
                "Meleagris gallopavo (turkey)" = "mgp",
                "Meriones unguiculatus (Mongolian gerbil)" = "mun",
                "Merops nubicus (carmine bee-eater)" = "mnb",
                "Mesitornis unicolor (brown roatelo)" = "mui",
                "Mesocricetus auratus (golden hamster)" = "maua",
                "Microcaecilia unicolor" = "muo",
                "Microcebus murinus (gray mouse lemur)" = "mmur",
                "Micropterus salmoides (largemouth bass)" = "msam",
                "Microtus fortis (reed vole)" = "mfot",
                "Microtus oregoni (creeping vole)" = "morg",
                "Miniopterus natalensis (Natal long-fingered bat)" = "mna",
                "Mirounga leonina (Southern elephant seal)" = "mlx",
                "Molossus molossus (Pallas's mastiff bat)" = "mmf",
                "Monodelphis domestica (gray short-tailed opossum)" = "mdo",
                "Monopterus albus (swamp eel)" = "malb",
                "Moschus berezovskii (Chinese forest musk deer)" = "mbez",
                "Mugil cephalus (flathead mullet)" = "mcep",
                "Mus caroli (Ryukyu mouse)" = "mcal",
                "Mus pahari (shrew mouse)" = "mpah",
                "Mustela putorius furo (domestic ferret)" = "mpuf",
                "Myotis brandtii (Brandt's bat)" = "myb",
                "Myotis davidii (David's myotis)" = "myd",
                "Myotis lucifugus (little brown bat)" = "mlf",
                "Myotis myotis (greater mouse-eared bat)" = "mmyo",
                "Myxocyprinus asiaticus (Chinese sucker)" = "masi",
                "Nannospalax galili (Upper Galilee mountains blind mole rat)" = "ngi",
                "Nanorana parkeri (Xizang Plateau frog)" = "npr",
                "Nematolebias whitei (Rio pearlfish)" = "nwh",
                "Neogale vison (American mink)" = "nvs",
                "Neomonachus schauinslandi (Hawaiian monk seal)" = "nsu",
                "Neophocaena asiaeorientalis asiaeorientalis (Yangtze finless porpoise)" = "nasi",
                "Nestor notabilis (Kea)" = "nnt",
                "Nipponia nippon (crested ibis)" = "nni",
                "Nomascus leucogenys (northern white-cheeked gibbon)" = "nle",
                "Notechis scutatus (mainland tiger snake)" = "nss",
                "Nothobranchius furzeri (turquoise killifish)" = "nfu",
                "Nothoprocta perdicaria (Chilean tinamou)" = "npd",
                "Notothenia coriiceps (black rockcod)" = "ncc",
                "Numida meleagris (helmeted guineafowl)" = "nmel",
                "Nyctereutes procyonoides (raccoon dog)" = "npo",
                "Ochotona princeps (American pika)" = "opi",
                "Odobenus rosmarus divergens (Pacific walrus)" = "oro",
                "Oenanthe melanoleuca (Eastern black-eared wheatear)" = "oma",
                "Oncorhynchus gorbuscha (pink salmon)" = "ogo",
                "Oncorhynchus kisutch (coho salmon)" = "oki",
                "Oncorhynchus mykiss (rainbow trout)" = "omy",
                "Oncorhynchus nerka (sockeye salmon)" = "one",
                "Oncorhynchus tshawytscha (Chinook salmon)" = "otw",
                "Onychostruthus taczanowskii (white-rumped snowfinch)" = "otc",
                "Opisthocomus hoazin (hoatzin)" = "oha",
                "Orcinus orca (killer whale)" = "oor",
                "Oreochromis aureus (blue tilapia)" = "oau",
                "Oreochromis niloticus (Nile tilapia)" = "onl",
                "Ornithorhynchus anatinus (platypus)" = "oaa",
                "Oryctolagus cuniculus (rabbit)" = "ocu",
                "Oryx dammah (scimitar-horned oryx)" = "oda",
                "Oryzias latipes (Japanese medaka)" = "ola",
                "Oryzias melastigma (Indian medaka)" = "oml",
                "Otolemur garnettii (small-eared galago)" = "oga",
                "Ovis aries (sheep)" = "oas",
                "Pan paniscus (bonobo)" = "pps",
                "Pan troglodytes (chimpanzee)" = "ptr",
                "Pangasianodon hypophthalmus (striped catfish)" = "phyp",
                "Panthera pardus (leopard)" = "ppad",
                "Panthera tigris altaica (Amur tiger)" = "ptg",
                "Panthera uncia (snow leopard)" = "puc",
                "Pantherophis guttatus (corn snake)" = "pgut",
                "Papio anubis (olive baboon)" = "panu",
                "Paralichthys olivaceus (Japanese flounder)" = "pov",
                "Paramormyrops kingsleyae" = "pki",
                "Parus major (Great Tit)" = "pmaj",
                "Passer montanus (Eurasian tree sparrow)" = "pmoa",
                "Pelecanus crispus (Dalmatian pelican)" = "pcri",
                "Pelodiscus sinensis (Chinese soft-shelled turtle)" = "pss",
                "Perca flavescens (yellow perch)" = "pflv",
                "Perognathus longimembris pacificus (Pacific pocket mouse)" = "plop",
                "Peromyscus leucopus (white-footed mouse)" = "pleu",
                "Phaethon lepturus (White-tailed tropicbird)" = "plet",
                "Phalacrocorax carbo (great cormorant)" = "pcao",
                "Phascolarctos cinereus (koala)" = "pcw",
                "Phasianus colchicus (Ring-necked pheasant)" = "pcoc",
                "Phocoena sinus (vaquita)" = "psiu",
                "Phodopus roborovskii (desert hamster)" = "prob",
                "Phyllostomus discolor (pale spear-nosed bat)" = "pdic",
                "Phyllostomus hastatus (greater spear-nosed bat)" = "phas",
                "Physeter catodon (sperm whale)" = "pcad",
                "Piliocolobus tephrosceles (Ugandan red Colobus)" = "pteh",
                "Pimephales promelas (fathead minnow)" = "pprm",
                "Pipistrellus kuhlii (Kuhl's pipistrelle)" = "pkl",
                "Plectropomus leopardus (leopard coralgrouper)" = "plep",
                "Podarcis muralis (common wall lizard)" = "pmua",
                "Poecilia formosa (Amazon molly)" = "pfor",
                "Poecilia latipinna (sailfin molly)" = "plai",
                "Poecilia mexicana (shortfin molly)" = "pmei",
                "Poecilia reticulata (guppy)" = "pret",
                "Pogona vitticeps (central bearded dragon)" = "pvt",
                "Polyodon spathula (Mississippi paddlefish)" = "pspa",
                "Polypterus senegalus (gray bichir)" = "psex",
                "Pongo abelii (Sumatran orangutan)" = "pon",
                "Prionailurus bengalensis (leopard cat)" = "pbg",
                "Propithecus coquereli (Coquerel's sifaka)" = "pcoq",
                "Protobothrops mucrosquamatus (Taiwan habu)" = "pmur",
                "Pseudonaja textilis (eastern brown snake)" = "ptex",
                "Pseudopodoces humilis (Tibetan ground-tit)" = "phi",
                "Pterocles gutturalis (yellow-throated sandgrouse)" = "pguu",
                "Pteropus alecto (black flying fox)" = "pale",
                "Pteropus giganteus (Indian flying fox)" = "pgig",
                "Pteropus vampyrus (large flying fox)" = "pvp",
                "Puma concolor (puma)" = "pcoo",
                "Puma yagouaroundi (jaguarundi)" = "pyu",
                "Pungitius pungitius (ninespine stickleback)" = "ppug",
                "Puntigrus tetrazona (Sumatra barb)" = "ptet",
                "Pygoscelis adeliae (Adelie penguin)" = "padl",
                "Pyrgilauda ruficollis (rufous-necked snowfinch)" = "pruf",
                "Python bivittatus (Burmese python)" = "pbi",
                "Rana temporaria (common frog)" = "rtem",
                "Rattus norvegicus (rat)" = "rno",
                "Rhincodon typus (whale shark)" = "rtp",
                "Rhinolophus ferrumequinum (greater horseshoe bat)" = "rfq",
                "Rhinopithecus bieti (black snub-nosed monkey)" = "rbb",
                "Rhinopithecus roxellana (golden snub-nosed monkey)" = "rro",
                "Rousettus aegyptiacus (Egyptian rousette)" = "ray",
                "Saimiri boliviensis boliviensis (Bolivian squirrel monkey)" = "sbq",
                "Salmo salar (Atlantic salmon)" = "sasa",
                "Salmo trutta (river trout)" = "stru",
                "Salvelinus namaycush (lake trout)" = "snh",
                "Salvelinus sp. IW2-2015 (Arctic char)" = "salp",
                "Sander lucioperca (pikeperch)" = "sluc",
                "Sarcophilus harrisii (Tasmanian devil)" = "shr",
                "Sceloporus undulatus (fence lizard)" = "sund",
                "Sciurus carolinensis (gray squirrel)" = "ncar",
                "Scleropages formosus (Asian bonytongue)" = "sfm",
                "Serinus canaria (common canary)" = "scan",
                "Seriola dumerili (greater amberjack)" = "sdu",
                "Seriola lalandi dorsalis (Yellowtail amberjack)" = "slal",
                "Silurus meridionalis (southern catfish)" = "smeo",
                "Siniperca chuatsi (mandarin fish)" = "schu",
                "Sinocyclocheilus anshuiensis" = "sanh",
                "Sinocyclocheilus grahami (golden-line barbel)" = "sgh",
                "Sinocyclocheilus rhinocerous" = "srx",
                "Solea senegalensis (Senegalese sole)" = "ssen",
                "Sorex araneus (European shrew)" = "sara",
                "Sphaerodactylus townsendi (Townsend's least gecko)" = "stow",
                "Strigops habroptila (Kakapo)" = "shab",
                "Struthio camelus australis (South African ostrich)" = "scam",
                "Sturnira hondurensis" = "shon",
                "Sturnus vulgaris (common starling)" = "svg",
                "Sus scrofa (pig)" = "ssc",
                "Syngnathus scovelli (Gulf pipefish)" = "sscv",
                "Tachysurus fulvidraco (yellow catfish)" = "tfd",
                "Taeniopygia guttata (zebra finch)" = "tgu",
                "Takifugu flavidus (sansaifugu)" = "tfs",
                "Takifugu rubripes (torafugu)" = "tru",
                "Talpa occidentalis (Iberian mole)" = "tod",
                "Tauraco erythrolophus (red-crested turaco)" = "teo",
                "Tetraodon nigroviridis (spotted green pufferfish)" = "tng",
                "Thamnophis sirtalis (common garter snake)" = "tsr",
                "Theropithecus gelada (gelada)" = "tge",
                "Tinamus guttatus (white-throated tinamou)" = "tgt",
                "Trachemys scripta elegans (red-eared slider)" = "tst",
                "Trachypithecus francoisi (Francois's langur)" = "tfn",
                "Trichechus manatus latirostris (Florida manatee)" = "tmu",
                "Triplophysa rosa" = "tros",
                "Tupaia chinensis (Chinese tree shrew)" = "tup",
                "Tympanuchus pallidicinctus (lesser prairie-chicken)" = "tpai",
                "Tyto alba (Barn owl)" = "tala",
                "Ursus americanus (American black bear)" = "uar",
                "Ursus arctos horribilis (North American brown bear)" = "uah",
                "Ursus maritimus (polar bear)" = "umr",
                "Varanus komodoensis (Komodo dragon)" = "vko",
                "Vicugna pacos (alpaca)" = "vpc",
                "Vulpes lagopus (Arctic fox)" = "vlg",
                "Vulpes vulpes (red fox)" = "vvp",
                "Xenopus laevis (African clawed frog)" = "xla",
                "Xenopus tropicalis (tropical clawed frog)" = "xtr",
                "Xiphias gladius (swordfish)" = "xgl",
                "Xiphophorus couchianus (Monterrey platyfish)" = "xco",
                "Xiphophorus hellerii (green swordtail)" = "xhe",
                "Xiphophorus maculatus (southern platyfish)" = "xma",
                "Zalophus californianus (California sea lion)" = "zca",
                "Zonotrichia albicollis (white-throated sparrow)" = "zab",
                "Zootoca vivipara (common lizard)" = "zvi"
              ),
              selected = "hsa",
              multiple = FALSE
            )
          ),
          )
        ),
         #),
        tags$tr(
          tags$td(width = "50%", 
          tags$div(style = "font-size:14px;", "Select analysis type")),
          tags$td(width = "50%",
        selectInput("pathType",
          label = NULL,
          choices = list(
            "Over representation" = 1,
            "GSEA" = 2
          ),
          selected = 1,
          multiple = FALSE
        ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select P-value cut-off for correction")),
          tags$td(width = "50%",
          numericInput("adjPvalue",
            label = NULL,
            value = 0.05
      ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select P-value correction method")),
          tags$td(width = "50%",      
          selectInput("pvaladj",
            label = NULL,
            choices = list(
              "Benjamini-Hochberg" = "BH",
              "Benjamini-Yeketuli" = "BY",
              "Bonferroni" = "bonferroni",
              "Holm" = "holm",
              "Hommel" = "hommel",
              "Hochberg" = "hochberg",
              "FDR" = "fdr",
              "None" = "none"
            ),
            selected = "BH"
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Upload data File"),
          tags$td(width = "50%",          
          fileInput("file7",
            label = NULL,
            accept = c(".csv",
                        ".txt",
                        ".xls",
                        ".xlsx")
          )))),
        ### select columns for logFC and gene symbol 
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select logFC column: "),
          tags$td(width = "50%",          
            selectInput("logfc_colm", label = NULL, "")
          ))),
        tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Select Gene Name column: "),
          tags$td(width = "50%",          
            selectInput("genesymbol_colm", label = NULL, "")
          ))), 
         verbatimTextOutput("genelist"),
         verbatimTextOutput("gns_col")
      ),
      conditionalPanel(
        condition = "input.pathType == 2",
        tags$table(width="100%",
          tags$tr(
            tags$td(width="50%", tags$div(style = "font-size:14pX;", "Minimum GeneSet size"),
          tags$td(width = "50%",          
          numericInput("mingssize", label = NULL, value = 10, min = 1)
          ))),
          tags$tr(
            tags$td(width="50%", tags$div(style = "font-size:14pX;", "Maximum GeneSet size"),
          tags$td(width = "50%",          
          numericInput("maxgssize", label = NULL, value = 500, min = 1)
          ))),
        )
      ),
      actionButton("Run", "Run Analysis",
      width = "75%",
        icon("paper-plane"),
        style = "color: #fff; background-color: #078082; border-color: #2e6da4; position:relative; left:calc(15%);"),
        br(),
        br(),
      HTML("<hr> <a href='pathway_example_data.txt' target='_blank'> <i class='fa fa-download'> </i> example data</a>"),
      tipify(actionButton("button7b", icon("question")),
        title = "To download the example file, right click and save the file.",
        placement = "right"
      )
    ),
    box(
      title = "Input File",
      closable = FALSE,
      width = 6,
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      label = boxLabel(
                      text = "Viewer", 
                      status = "primary",
                      style = "default"
                    ),
      dropdownMenu = boxDropdown(
                          icon = shiny::icon("file-pdf"), 
                          boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
                          ),      
      DT::dataTableOutput("pathwaydatafile", width = "100%", height = "500px")
    ),
    
    box(
      title = "Analysis Result - Figure",
      status = "navy",
      width = 6,
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = TRUE,
      closable = FALSE,
      label = boxLabel(
            text = "enrichF", # Change to number when finished all package uploads
            status = "danger",
            style = "default"
          ),
      dropdownMenu = boxDropdown(
        icon = shiny::icon("file-pdf"),
        boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
      ),
      sidebar = boxSidebar(
        startOpen = TRUE,
        width = 34,
        height = "300px",
        id = "pathfigparams",
        background = "#05282e",
        icon = shiny::icon("user-gear"),
        tags$head(tags$style(HTML("#pathfigparams+ div>.selectize-input {max-height: 300px}"))), 
        tags$table(width="90%",
          tags$tr(
            tags$td(width="50%", tags$div(style = "font-size:14pX;", "Type of figure"),
          tags$td(width = "50%",          
          selectInput("pathfigtype", label = NULL, choices = c(
            "Scatter" = "scatter",
            "Bar" = "bar"),
            selected = "scatter",
            multiple = FALSE)
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Color for high p-value"),
          tags$td(width = "50%",          
          colourpicker::colourInput("lowcolopath", label = NULL, "blue"),
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Color for mid p-value"),
          tags$td(width = "50%",          
          colourpicker::colourInput("midcolopath", label = NULL, "white"),
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Color for low p-value (towards 0)"),
          tags$td(width = "50%",          
          colourpicker::colourInput("highcolopath", label = NULL, "red"),
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Write x-axis title"),
          tags$td(width = "50%",          
          textInput("pathxaxistitle", label = NULL, value = "Gene ratio")
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Write y-axis title"),
          tags$td(width = "50%",          
          textInput("pathyaxistitle", label = NULL, value = "Description")
          ))),
          tags$tr(
          tags$td(width = "50%", tags$div(style = "font-size:14pX;", "Figure main title"),
          tags$td(width = "50%",          
          textInput("pathmaintitle", label = NULL, value = "TranscriptR pathway analysis")
          )))        
          
          )
          ),
        plotlyOutput("pathPlot", width = "100%", height = "960px")
      ),

    box(
      title = "Analysis Result - Table",
      status = "navy",
      width = 6,
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = TRUE,
      closable = FALSE,
      label = boxLabel(
            text = "enrichT", # Change to number when finished all package uploads
            status = "danger",
            style = "default"
          ),
          dropdownMenu = boxDropdown(
            icon = shiny::icon("file-pdf"),
            boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
          ),
          DTOutput("pathtable")
      ),
      conditionalPanel(
      condition = "input.pathDB == 2",
      box(
      title = "KEGG Pathway Map with Parameters setup",
      closable = FALSE,
      collapsed = TRUE,
      width = 12,
      height = "150px",
      status = "navy",
      solidHeader = TRUE,
      collapsible = TRUE,
      label = boxLabel(
        text = "pathmap", # Change to number when finished all package uploads
        status = "danger",
        style = "default"
      ),
      dropdownMenu = boxDropdown(
        icon = shiny::icon("file-pdf"),
        boxDropdownItem("Link to Manual", href = "#", icon = icon("file-pdf"))
      ),
      sidebar = boxSidebar(
        startOpen = TRUE,
        width = 34,
        #height = "500px",
        id = "pathmapParams",
        background = "#1D5E6B",
        icon = shiny::icon("user-gear"),
        #textInput("pathName", "KEGG Pathway ID"),
        #selectInput("pathName", "KEGG Pathway ID:", ""),
        colourpicker::colourInput("lowcolopathmap", "Select Color for down-regulated genes", "blue"),
        colourpicker::colourInput("midcolopathmap", "Select Color for non-significant genes", "white"),
        colourpicker::colourInput("highcolopathmap", "Select Color for up-regulated genes", "red"),
        actionButton("runPM", "Show Pathway Map"),
        br(),
        br(),
        h4("Save multiple pathway map together"),
        textInput("dirname", "Create directory name"),
        actionButton("copyfilepm", "Copy All Pathway Map files"),
        h4("Download all pathway maps as zip file", tipify(actionButton("download_pmdata", icon("download")),
          "Click the button only when the process finished.",
          placement = "bottom"
        )),
        actionButton("showpathmapgene",
                                  "Show gene details",
                                  width = "400px",
                                  icon("paper-plane"), 
                        style="color: #000000; background-color: #fac228; border-color: #2e6da4"),
      ),
      imageOutput("pathmap_result", height = "800px", width = "100%")
    ) 
  ))
)

#========================================#
## User Manual
#========================================#
bodyManual <- tabItem(tabName = "manual",

                      box(
                        title = "Online Manual",
                        status = "navy",
                        width = 12,
                        #height = "380px",
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        collapsed = FALSE,
                        closable = FALSE,
                        tags$a(href= "https://methylr.netlify.app/intro.html", "User Manual", target = "_blank")
                      ),
                      box(width = 12, 
                              status = "navy",
                              title = "PDF manual",
                              
                           solidHeader = TRUE,
                              collapsible = FALSE,
                              collapsed = FALSE,
                              closable = FALSE, 
                              div(style = "font-size:20px; color: red",
                                  tags$iframe(style="height:800px; width:100%", 
                                              src="https://drive.google.com/file/d/1za4rS344DWlEKZawidkNkLKgx-rAhxL7/preview")
                                  #tags$a(href= "https://methylr.netlify.app/intro.html", "User Manual", target = "_blank")
                              )
                          )
              )

rightsidebar = dashboardControlbar(
  skin = "dark",
  controlbarMenu(
    id = "menu",
    controlbarItem(
      "Tab 1",
      "Visitors",
      tags$body(HTML("<div style='display:inline-block;width:200px;'><script type='text/javascript' src='//rf.revolvermaps.com/0/0/6.js?i=5vfrhre4qej&amp;m=0&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33' async='async'></script></div>"))
    ),
    controlbarItem(
      "Tab 2",
      "Welcome to tab 2"
    )
  )
)




  ui = dashboardPage(
    #skin = "black", 
    header = dashboardHeader(
      title = "TranscriptR",    
      userOutput("user"),
      leftUi = tagList(
          dropdownBlock(
            id = "mydropdown",
            #icon = icon("vault"),
            #badgeStatus = "danger",
            title = "Disclaimer",
            div(style = "font-size:10px; font-weight: bold; color:red; font-family: monospace",
            HTML("WARNING - Please note we included the possibility of the human transcriptome data analysis from the sequence data. None of the owner of the code or data center is responsible for handling the human data. It is upto the user if they want to use the online version for human data analysis."), br(),
            tags$a(href = "https://hub.docker.com/jd21/transcriptr", "Otherwise the user can download the docker container and use it.", target= "_blank"), br(), 
            tags$a(href = "https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiK1t-O35uBAxUCQ_EDHaXdCWoQFnoECCsQAQ&url=https%3A%2F%2Fwww.phgfoundation.org%2Fmedia%2F123%2Fdownload%2Fgdpr-and-genomic-data-report.pdf%3Fv%3D1&usg=AOvVaw0UnPevgJgQUtlOm3SUuOx7&opi=89978449", "Please read more here", target= "_blank")
            )
          )
        )
      #)
      ),
      
    sidebar,
    rightsidebar,
    body = dashboardBody(
      tags$head(tags$style(HTML('
        .skin-black .main-header .logo {
          background-color: #012D5A;
        }
        .skin-black .main-header .logo:hover {
          background-color: #012D5A;
        }
      '))),
      tags$head(tags$style(HTML('
  .navbar-custom-menu>.navbar-nav>li>.dropdown-menu {
  width:400px;
  }
  '))),
      tabItems(
        bodyMain,
        bodyQC,     
        bodyTranscrip,
        bodyIGV,
        bodyMDA,
        bodyHeat,
        bodyVol,
        bodyGo,
        bodyPathway,
        bodyVenn,
        bodyUpSet,
        bodyManual
        
      )
    ),
    title = "TranscriptR",
    footer = dashboardFooter(
      left = "Copyright © 2023, Massimiliano Volpe and Jyotirmoy Das,",
      right = "Linköping Unversity, Sweden"
    )
  )
  

server = function(input, output, session) {
    options(shiny.maxRequestSize=3000*1024^2)
    
    ## helpfiles

    shinyhelper::observe_helpers(help_dir = "www/", withMathJax = TRUE)

    ## Next-Previous button,
    observeEvent(input$next1, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "sinfo")
  })

  observeEvent(input$next2, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "paramss")
  })

  observeEvent(input$back1, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "fqupload")
  })

  observeEvent(input$back2, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "sinfo")
  })

    shinyjs::runjs("window.scrollTo(0, 50)")

    runjs("
      $('.box').on('click', '.box-header h3', function() {
          $(this).closest('.box')
                 .find('[data-widget=collapse]')
                 .click();
      });")

    output$user <- renderUser({
      dashboardUser(
        name = "TranscriptR", 
        #image = "https://avatars.githubusercontent.com/u/15922543?v=4", 
        image = "https://images.unsplash.com/photo-1633167606207-d840b5070fc2?ixlib=rb-4.0.3&ixid=M3wxMjA3fDB8MHxzZWFyY2h8Mnx8ZG5hfGVufDB8fDB8fHww&auto=format&fit=crop&w=500&q=60",
        title = "TranscriptR",
        subtitle = "Shiny application",
        descriptionBlock(
              header = "Product License", 
              text = tags$a(href="https://www.gnu.org/licenses/gpl-3.0.en.html#license-text", target = "_blank", "Check Here!"), 
              #text = tags$a(href="https://hackmd.io/@JD2112/BJt_mRFk2", target = "_blank", "Check Here!"), 
              rightBorder = FALSE,
              marginBottom = TRUE
            ),        
        fluidRow(
          dashboardUserItem(
            width = 6,
            socialButton(
              href = "https://dropbox.com",
              icon = icon("dropbox")
            )
          ),
          dashboardUserItem(
            width = 6,
            socialButton(
              href = "https://github.com/JD2112/transcriptr",
              icon = icon("github")
            )
          )
        )
      )
    })

### human data disclaimer
observe({
if (req(input$transcriptrmenu == "qc" || input$transcriptrmenu == "transana")){
  message("Quality Analysis has been selected")
  showModal(modalDialog(
    title = HTML('<span style="color:red; font-size: 21px; font-weight:bold; font-family:sans-serif ">WARNING - Please note we included the possibility of the human transcriptome data analysis from the sequence data. None of the owner of the code or data center is responsible for handling the human data. It is upto the user if they want to use the online version for human data analysis.<span>
                                <button type = "button" class="close" data-dismiss="modal" ">
                                </button> '),
    tags$div(tags$a(href = "https://hub.docker.com/jd21/transcriptr", "Otherwise the user can download the docker container and use it.", target= "_blank")
    ),
    tags$div(tags$a(href = "https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiK1t-O35uBAxUCQ_EDHaXdCWoQFnoECCsQAQ&url=https%3A%2F%2Fwww.phgfoundation.org%2Fmedia%2F123%2Fdownload%2Fgdpr-and-genomic-data-report.pdf%3Fv%3D1&usg=AOvVaw0UnPevgJgQUtlOm3SUuOx7&opi=89978449", "Please read more on GDPR and other terms & conditons", target= "_blank")
    ),
                      footer = tagList(
                      actionButton("notagree", "Do not agree"),
                      actionButton("agree", "Agree to the terms & conditions")
                      ),

                      easyClose = FALSE

  ))

#   observeEvent(input$agree, {
#   removeModal()
#   md_out <- rmarkdown::render("www/TermsofUse.md")
#   showModal(modalDialog(
#     title = renderUI(HTML(readLines(md_out))),

#     footer = 
#                       #actionButton("notagree2", "Do not agree"),
#                       actionButton("agree2", "Agree to the terms & conditions")
#                       ,
#     easyClose = FALSE
#   ))
  
# })

observeEvent(input$agree, {
   removeModal()
  #updateSelectInput(session, "species", selected = 'homo_sapiens')
 })

observeEvent(input$notagree, {
  if(input$notagree){
    withProgress(
      message = "Thank you for your interest. Unfortunaly unable to process your request.",
      detail = "Session will reload in 5 seconds",
      value = 0,
      {        
        for (i in 1:15) {
          incProgress(1 / 15)
          Sys.sleep(0.25)
        }      
      }
    )
    session$reload()
  } else {}
   removeModal()

 })
}
})
###========

### main frame

output$frame <- renderUI(
  tags$iframe(src = "https://hackmd.io/@JD2112/HycD9CfWn", height = 380, width = "100%")
)

 #### ========Quality Control Analysis module start==================#####
  output$fastq_files1 <- renderText({
    # Test if file is selected
    if (!is.null(input$file_transcriptome1$datapath)) {
      filenames <- gsub(".fastq$", "", basename(input$file_transcriptome1$name))
      namefiles <- gsub("\\s", "\\\n]", filenames)

      file.copy(input$file_transcriptome1$datapath, paste0(getwd(),"/", input$file_transcriptome1$name), overwrite = TRUE)

      return(namefiles)
    } else {
      return(NULL)
    }
  })

  json_qc_pe <- reactive({
    toJSON(
      list(
        data = paste0(getwd(),"/"),
        workdir = "/transcriptr/workdir/",
        R1_suffix = input$R1suffix_qc_pe,
        R2_suffix = input$R2suffix_qc_pe
      ),
      pretty = TRUE
    )
  })

  json_qc_se <- reactive({
    toJSON(
      list(
        data = paste0(getwd(),"/"),
        workdir = "/transcriptr/workdir/",
        R1_suffix = input$R1suffix_qc_se
        
      ),
      pretty = TRUE
    )
  })

  output$jsonview1 <- renderPrint({
    if(input$qc_endtype == 'qc_send'){
      req(json_qc_se())} else {
      req(json_qc_pe())
      }
    
  })

  output$downloadjson <- downloadHandler(
    filename = function() {
      paste0("fastqc", ".json", sep = "")
    },
    content = function(file) {
      if(input$qc_endtype == 'qc_send'){write(json_qc_se(), file)} else {
      write(json_qc_pe(), file)
    }
    }
  )


  observeEvent(input$showqc1, {
    showModal(
      modalDialog(
        HTML('<span style="color:LightSeaGreen; font-size: 26px; font-weight:bold; font-family:sans-serif ">Are you sure? Please check that you click to send the JSON file to the server<span>
                    <button type = "button" class="close" data-dismiss="modal" ">
                    </button> '),
        footer = tagList(
          actionButton(inputId = "runqc", label = "Yes"),
          modalButton("No")
        ),
        easyClose = FALSE
      )
    )
  })

  observeEvent(input$runqc, {
    removeModal()
    updateActionButton(inputId = "showqc1", label = "Running...", icon = animateIcon("spinner", "spin", "fast"))

    if(input$qc_endtype == 'qc_send'){
      jsonqc1 <- gsub("\\[|\\]", "", json_qc_se())
      write(jsonqc1, "/srv/shiny-server/fastqc_single.json")
    } else {
      jsonqc1 <- gsub("\\[|\\]", "", json_qc_sp())
      write(jsonqc1, "/srv/shiny-server/fastqc.json")      
    }

    withProgress(
      message = "Running analysis... results will be displayed shortly.",
      detail = "This message disappears in 5 seconds",
      value = 0,
      {        
        for (i in 1:15) {
          incProgress(1 / 15)
          Sys.sleep(0.25)
        }      
      }
    )

    if(input$qc_endtype == 'qc_send'){
    system("snakemake -p --cores 12 -s /transcriptr/workflow/fastqc_single.smk --configfile /srv/shiny-server/fastqc_single.json")} else {
    system("snakemake -p --cores 12 -s /transcriptr/workflow/fastqc.smk --configfile /srv/shiny-server/fastqc.json")
    }
  })

observe({
file.copy("/transcriptr/workdir/results_qc/multiqc_report.html", "/srv/shiny-server/www/multiqc_report.html", overwrite = TRUE) # For Docker

#file.copy("/mnt/WD1/test/biomartr/results/multiqc_report.html", "www/multiqc_report.html", overwrite = TRUE)
})


  output$qcplot1 <- renderUI({
    if (input$showqc2) {
      tags$iframe(src = "multiqc_report.html", height = 800, width = "100%")
    } else {}
  })

# observe({
#   unlink("/transcriptr/workdir/results", recursive = TRUE)
# })

output$consoleTextqc <- renderPrint({
  message("Getting log, please wait")
    library(readr)
    library(stringr)
    library(dplyr)
    filepath <- list.files("/var/log/shiny-server/", full.names =T, pattern = ".log")
    #filepath <- list.files("/home/jyotirmoy/Downloads/logFile", full.names =T, pattern = ".log")
    txtfile <- readLines(textConnection(filepath[[2]]))
    read_file(txtfile) -> mylog
    tt1 <- unlist(strsplit(mylog, '\n') )
    ix1 <- tt1[grep("steps | failed", tt1)]
    #write(message, ix1, append = TRUE)
    return(ix1)

  })

output$downloadlog <- downloadHandler(
  filename = function(){
    paste("TranscriptR","_",input$downloadlog, sep="")
  },
  content = function(file){
    filepath <- list.files("/var/log/shiny-server/", full.names =T, pattern = ".log")
    #filepath <- list.files("/home/jyotirmoy/Downloads/logFile", full.names =T, pattern = ".log")
    txtfile <- readLines(textConnection(filepath[[2]]))
    read_file(txtfile) -> mylog
    
    writeLines(mylog)
  }
)
  #### ========Quality Control Analysis module End==================#####


#####============== Transcriptome module Start==============#######

transcriptR <- function(){

dir.create(paste0(tempdir(),"/transcriptrInput"))
dir.create(paste0(tempdir(),"/transcriptr_results"))
# Check sample navy file
file2trans = reactive({
    filetrans2 <- input$file_transmeta
    ext <- tools::file_ext(filetrans2$datapath)
    req(filetrans2)
    
    filetrans2 <- read.delim(filetrans2$datapath,
                             sep = "\t", 
                             header = T)

    file.copy(input$file_transmeta$datapath, paste0(getwd(),"/", input$file_transmeta$name), overwrite = TRUE)    

    return(filetrans2)
  }
  )


  output$transsamplefile <- DT::renderDataTable(
    file2trans(), rownames = FALSE, editable = TRUE,
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    )
    
  )


# Group Comparison file
file2comp = reactive({
    filecomp2 <- input$file_transcomp
    ext <- tools::file_ext(filecomp2$datapath)
    req(filecomp2)
    
    filecomp2 <- read.delim(filecomp2$datapath,
                             sep = "\t", 
                             header = T)

    file.copy(input$file_transcomp$datapath, paste0(getwd(),"/", input$file_transcomp$name), overwrite = TRUE) 

    return(filecomp2)
  }
  )

  output$transcompfile <- DT::renderDataTable(
    file2comp(), rownames = FALSE, editable = TRUE,
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    )
    
  )

# Input fastq files
  output$fastq_files <- renderText({
  # Test if file is selected
  if (!is.null(input$file_transcriptome$datapath)) {
      filenames = gsub(".fastq$", "", basename(input$file_transcriptome$name))
      namefiles = gsub("\\s","\\\n]", filenames)

      file.copy(input$file_transcriptome$datapath, paste0(getwd(),"/", input$file_transcriptome$name), overwrite = TRUE)

      return(namefiles)
  } else {
      return(NULL)
  }
})


# parameters settings JSON file

dbver = reactive({
    if(input$dbversion == 'latest'){paste(input$dbversion)} else {
      paste(as.character(input$customdbversion))
    }
})

# R1suffix = reactive({
#   if(input$endtype == 'send'){paste(input$R1suffix_se)} else {
#     paste(input$R1suffix_pe)
#   }
# })

# R2suffix = reactive({
#   if(input$endtype == 'pend'){paste(input$R2suffix_pe)} else {
#     NULL
#   }
# })

json_data_pe  <- reactive({
    toJSON(
      list(
        data = paste0(getwd(),"/"),
        workdir = "/transcriptr/workdir/",
        rearrangeFC = "/transcriptr/workflow/scripts/rearrangeFeatureCounts.R",
        edgeR = "/transcriptr/workflow/scripts/edger.R",
        DESeq2 = "/transcriptr/workflow/scripts/deseq2.R",
        limma = "/transcriptr/workflow/scripts/limma.R",
        refR = "/transcriptr/workflow/scripts/ref.R",
        sampleInfo = paste0(getwd(),"/",input$file_transmeta[1]),
        compsTab = paste0(getwd(),"/",input$file_transcomp[1]),
        R1_suffix = input$R1suffix_pe,
        R2_suffix = input$R2suffix_pe,      
        species = input$species,
        release = dbver(),
        overhang = as.character(input$overh),
        attribute = input$attribute,
        stranded = input$seqtype,
        cpm = as.character(input$cpm),
        nsamples = as.character(input$nsamples),
        logFC = as.character(input$logFC),
        FDR = as.character(input$FDR)#,
        #pval = as.character(0)
      ),
      pretty = TRUE
    )})

 json_data_se  <- reactive({
    toJSON(
      list(
        data = paste0(getwd(),"/"),
        workdir = "/transcriptr/workdir/",
        rearrangeFC = "/transcriptr/workflow/scripts/rearrangeFeatureCounts.R",
        edgeR = "/transcriptr/workflow/scripts/edger.R",
        DESeq2 = "/transcriptr/workflow/scripts/deseq2.R",
        limma = "/transcriptr/workflow/scripts/limma.R",
        refR = "/transcriptr/workflow/scripts/ref.R",
        sampleInfo = paste0(getwd(),"/",input$file_transmeta[1]),
        compsTab = paste0(getwd(),"/",input$file_transcomp[1]),
        R1_suffix = input$R1suffix_se,
        #R2_suffix = input$R2suffix_pe,      
        species = input$species,
        release = dbver(),
        overhang = as.character(input$overh),
        attribute = input$attribute,
        stranded = input$seqtype,
        cpm = as.character(input$cpm),
        nsamples = as.character(input$nsamples),
        logFC = as.character(input$logFC),
        FDR = as.character(input$FDR)#,
        #pval = as.character(0)
      ),
      pretty = TRUE
    )})

  output$jsonview <- renderPrint({
    if(input$endtype == 'send'){req(json_data_se())} else {
      req(json_data_pe())
    }
    
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0("testconfig", ".json", sep = "")
    },
    content = function(file) {
      if(input$endtype == 'send'){write(json_data_se(), file)} else {
        write(json_data_pe(), file)
      }
    }
  )
  
observeEvent(input$testtranscript, {
  showModal(
              modalDialog(
                  HTML('<span style="color:LightSeaGreen; font-size: 26px; font-weight:bold; font-family:sans-serif ">Are you sure? Please check that you click to send the JSON file to the server<span>
                    <button type = "button" class="close" data-dismiss="modal" ">
                    </button> '),
                  footer = tagList(
                      actionButton(inputId = "runtranscriptr", label = "Yes"),
                      modalButton("No")
                  ),
                  easyClose = FALSE
              )
          ) 
})

observeEvent(input$runtranscriptr,{
  removeModal()

  unlink("/transcriptr/workdir/results", recursive = TRUE) # remove qc result directory

  updateActionButton(inputId = "testtranscript", label = "Running...")
  # debounce to display the message
  if(input$endtype == 'send'){
      jsondata1 <- gsub("\\[|\\]", "", json_data_se())   
      write(jsondata1, "/srv/shiny-server/config_single.json")
    } else {
      jsondata1 <- gsub("\\[|\\]", "", json_data_pe())   
      write(jsondata1, "/srv/shiny-server/config.json")
    }
  

  withProgress(message = "Pipeline will run in the background... \nit will take few minutes to hours depending on number of samples.\n Once done results are available on other tabs", 
  detail = "This message disappears in 5 seconds", value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
  })
  
  if(input$endtype == 'send'){
      system("snakemake -p --cores 12 -s /transcriptr/workflow/transcriptr_single.smk --configfile /srv/shiny-server/config_single.json")
    } else {
      system("snakemake -p --cores 12 -s /transcriptr/workflow/transcriptr.smk --configfile /srv/shiny-server/config.json")
    }
  
  })

# }

# shinyCatch(transcriptR(), position = "bottom-right", blocking_level = "none", shiny = TRUE, prefix = "transcriptR", trace_back = spsOption("traceback"))

observe({
file.copy("/transcriptr/workdir/results/multiqc_report.html", "/srv/shiny-server/www/multiqc_report.html", overwrite = TRUE) # For Docker

#file.copy("/mnt/WD1/test/biomartr/results/multiqc_report.html", "www/multiqc_report.html", overwrite = TRUE)
})


output$qcplot <- renderUI(
  tags$iframe(src = "multiqc_report.html", height = 800, width = "100%")
)

# Read and display all PDF plots
observeEvent(input$showmdsplot, {
  library(qpdf)

  if (input$degplotsel == "edger"){
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/edgeR_results/MDS_Pre-Norm.pdf", "/transcriptr/workdir/results/edgeR_results/MDS_Post-Norm.pdf"), ## docker PATH
  #qpdf::pdf_combine(input = c("/mnt/WD1/test/biomartr/results/edgeR_results/MDS_Pre-Norm.pdf", "/mnt/WD1/test/biomartr/results/edgeR_results/MDS_Post-Norm.pdf"),
                  output = "www/MDS_edgeR.pdf")
  } else if(input$degplotsel == "deseq") {
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/deseq2_results/MDS_plot.pdf", "/transcriptr/workdir/results/deseq2_results/Dispersion_plot.pdf"),
                  output = "www/MDS_DESeq2.pdf")
  } else {
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/limma_results/MDS_Pre-Norm.pdf", "/transcriptr/workdir/results/limma_results/MDS_Post-Norm.pdf"),
                  output = "www/MDS_limma.pdf")
  }


})
output$mdsplottrans <- renderUI({
  if (input$degplotsel == "edger"){
  tags$iframe(style="height:700px; width:100%", src="MDS_edgeR.pdf")
  } else if (input$degplotsel == "deseq") {
    tags$iframe(style="height:700px; width:100%", src="MDS_DESeq2.pdf")
  } else {
    tags$iframe(style="height:700px; width:100%", src="MDS_limma.pdf")
  }
})


observeEvent(input$showboxplot, {
  library(qpdf)
  if (input$degplotsel == "edger"){
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/edgeR_results/Boxplot_Pre-Norm.pdf", "/transcriptr/workdir/results/edgeR_results/Boxplot_Post-Norm.pdf"), ## Docker PATH
                  output = "www/boxplot_edgeR.pdf")

  #qpdf::pdf_combine(input = c("/mnt/WD1/test/biomartr/results/edgeR_results/Boxplot_Pre-Norm.pdf", "/mnt/WD1/test/biomartr/results/edgeR_results/Boxplot_Post-Norm.pdf"), 
  #                output = "www/boxplot_edgeR.pdf")

  } else if (input$degplotsel == "deseq") {  
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/deseq2_results/Boxplot_Pre-Norm.pdf", "/transcriptr/workdir/results/deseq2_results/Boxplot_Post-Norm.pdf"),
                  output = "www/boxplot_deseq2.pdf")       
  } else {
  qpdf::pdf_combine(input = c("/transcriptr/workdir/results/limma_results/Boxplot_Pre-Norm.pdf", "/transcriptr/workdir/results/limma_results/Boxplot_Post-Norm.pdf"),
                  output = "www/boxplot_limma.pdf")  
  }

})
output$boxplottrans <- renderUI({
  if (input$degplotsel == "edger"){
  tags$iframe(style="height:700px; width:100%", src="boxplot_edgeR.pdf")
  } else if (input$degplotsel == "deseq") {
    tags$iframe(style="height:700px; width:100%", src="boxplot_deseq2.pdf")
  } else {
    tags$iframe(style="height:700px; width:100%", src="boxplot_limma.pdf")
  }
})

# Read and store DEG file for output
library(readxl)

filedeg = reactive({
      id_degtab <- showNotification("loading, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_degtab), add = TRUE)

    if(input$degtablesel == 'edger'){
      filedeg1 <- input$showdegtab1
    req(filedeg1)
                        
    filedeg1 <- read.delim(file.path("/transcriptr/workdir/results/edgeR_results",input$showdegtab1),
                             sep = "\t", 
                             header = T)

    #return(filedeg2)
    shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }

  degfile <- data.frame(filedeg1)
  filedeg_rows <- nrow(degfile)
  
  degfile$Gene_details <- shinyInput(actionButton, filedeg_rows, 'button_', label = "Show details", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})' )

    tableDEG  = data.frame(
          ID = as.character(degfile$id), 
          logFC = as.character(degfile$logFC),
          logCPM = as.character(degfile$logCPM),
          P_Value = as.character(degfile$PValue),
          FDR = as.character(degfile$FDR),          
          Gene_Symbol = as.character(degfile$name),
          Ensembl_Link = sapply(degfile$id, function(x)
            #toString(tags$a(href=paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", x),target = "_blank", x))), 
            toString(tags$a(href=paste0("http://www.ensembl.org/Multi/Search/Results?q=", x),target = "_blank", x))), 
          GeneCard_Link = sapply(degfile$name, function(x)
            toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", x),target = "_blank", x))), 
          Gene_Details = degfile$Gene_details,
          Chr = as.character(degfile$chr),
          Start = degfile$start,
          End = degfile$end,
          Strand = degfile$strand,
          Types = degfile$biotype 
          )

    return(tableDEG)

    } else if (input$degtablesel == 'deseq') {
      filedeg2 <- input$showdegtab2
    req(filedeg2)
                   
    filedeg2 <- read.delim(file.path("/transcriptr/workdir/results/deseq2_results",input$showdegtab2),
                             sep = "\t", 
                             header = T)

    shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }

  degfile <- data.frame(filedeg2)
  filedeg_rows <- nrow(degfile)
  
  degfile$Gene_details <- shinyInput(actionButton, filedeg_rows, 'button_', label = "Show details", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})' )

    tableDEG  = data.frame(
          ENSEMBL_geneID = sapply(degfile$id, function(x)
            #toString(tags$a(href=paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", x),target = "_blank", x))), 
            toString(tags$a(href=paste0("http://www.ensembl.org/Multi/Search/Results?q=", x),target = "_blank", x))), 
          baseMean = as.character(degfile$baseMean),  
          logFC = as.character(degfile$log2FoldChange),
          lfcSE = as.character(degfile$lfcSE),
          stat = as.character(degfile$stat),
          P_Value = as.character(degfile$pvalue),
          FDR = as.character(degfile$padj),
          Gene_Symbol = as.character(degfile$name),
          GeneCard_Link = sapply(degfile$name, function(x)
            toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", x),target = "_blank", x))), 
          Gene_Details = degfile$Gene_details,          
          Chr = as.character(degfile$chr),
          Start = degfile$start,
          End = degfile$end,
          Strand = degfile$strand,
          Types = degfile$biotype 
          )

    return(tableDEG)
    } else {
      filedeg3 <- input$showdegtab3
    req(filedeg3)
                       
    filedeg3 <- read.delim(file.path("/transcriptr/workdir/results/limma_results",input$showdegtab3),
                             sep = "\t", 
                             header = T)

    shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }



  #degfile <- data.frame(filedeg2)
  degfile <- filedeg3
  filedeg_rows <- nrow(degfile)
  
  degfile$Gene_details <- shinyInput(actionButton, filedeg_rows, 'button_', label = "Show details", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})' )

    tableDEG  = data.frame(
          ENSEMBL_geneID = sapply(degfile$id, function(x)
            #toString(tags$a(href=paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", x),target = "_blank", x))), 
            toString(tags$a(href=paste0("http://www.ensembl.org/Multi/Search/Results?q=", x),target = "_blank", x))), 
          #baseMean = as.character(degfile$baseMean),  
          logFC = as.character(degfile$logFC),
          AveExpr = as.character(degfile$AveExpr),
          t = as.character(degfile$t),
          P_Value = as.character(degfile$P.Value),
          FDR = as.character(degfile$adj.P.Val),
          B = as.character(degfile$B),
          Gene_Symbol = as.character(degfile$name),
          GeneCard_Link = sapply(degfile$name, function(x)
            toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", x),target = "_blank", x))), 
          Gene_Details = degfile$Gene_details,
          Chr = as.character(degfile$chr),
          Start = degfile$start,
          End = degfile$end,
          Strand = degfile$strand,
          Types = degfile$biotype 
          )

    return(tableDEG)
    }
    
  }
  )


### HGNC connection

output$geneframe <- renderUI({
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }

  row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
  defile <- data.frame(filedeg())
  genesdeg <- defile$Gene_Symbol[[row]]
  
  test1 <<- paste0("https://www.genenames.org/tools/search/#!/?query=",genesdeg)
  if (is.null(row) || row == '') {
    tags$iframe(src="www.google.com", height = "400px", width = "100%")
    } else {
    gene_test <- tags$iframe(src=test1, height = "400px", width = "100%")
    print(gene_test)
    gene_test
    }
})




### CPM plot


output$cpmplot <- renderPlotly({

  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
# dat <- read.delim("/transcriptr/workdir/results/edgeR_results/raw_counts_cpm_after_norm.txt", sep = "\t", header = T) # require for Docker
dat <- read.delim("/transcriptr/workdir/results/edgeR_results/raw_counts_cpm_after_norm.txt", sep = "\t", header = T) 

dat1 <- data.frame(dat)
dat.col <- ncol(dat1)
dat2 <- data.frame(rownames(dat1), dat1)
colnames(dat2)[1] <- "ENSGID"
rownames(dat2) <- NULL

#degfile1 <- data.frame(filedeg())
degfile1 <- filedeg()
colnames(degfile1)[1] <- "ENSGID"

matchdeg <- merge(dat2, degfile1, by = "ENSGID")
matchdeg1 <- data.frame(matchdeg[,1:dat.col+1])
rownames(matchdeg1) <- matchdeg[,1]
tdat <- data.frame(t(matchdeg1))
tdat$Samples <- rownames(tdat)

### Assign group info

    filetrans2 <- input$file_transmeta
    ext <- tools::file_ext(filetrans2$datapath)
    req(filetrans2)
    
    filetrans2 <- read.delim(filetrans2$datapath,
                             sep = "\t", 
                             header = T)

trfile <- data.frame(filetrans2)

tdat$groups <- trfile$condition

#tdat$groups <- 

row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
genesdeg <- tdat[[row]]

library(plotly)

tdat %>%
  plot_ly(x = ~groups,y = ~genesdeg, 
          color = ~groups, 
          type = "box",
          pointpos = 0,
          size = 12,
          hoverinfo = "text",
          text = ~paste0("Group: ",groups,                         
                        "<br>Sample Name: ",tdat$Samples),
          boxpoints = "all",
          showlegend = FALSE) %>% config(displaylogo = FALSE) %>%
          layout(yaxis = list(title = 'CPM value'))
})


### DEG table output
 output$degtable <- DT::renderDataTable(server = FALSE, {

     DT::datatable(
          filedeg(),
          callback = JS(c(
          "$('#gotopage').on('click', function(){",
          "  var page = parseInt($('#page').val())-1;",
          "  table.page(page).draw('page');",
          "});"
        )),
          extensions = c("FixedColumns","Buttons"),
          escape = FALSE,
          rownames = FALSE,
          filter = "bottom",
          style = "auto",
          selection = "single",
          #editable = TRUE,
          options = list(
            paging = TRUE,
            searching = TRUE,
            autowidth = TRUE,
            ordering = TRUE,
            dom = 'Bflrtip',
            scrollX = TRUE,
            lengthMenu = list(c(10, 50,100, -1), c('10', '50', '100','All')),
            fixedColumns = list(leftColumns = 2),
            buttons = list(
              list(extend = "excel", text = "Download Table", title = "", filename = paste0("TranscriptR_","DEG_Results_", input$showdegtab, sep=""),
                   exportOptions = list(
                     modifier = list(page = "current")
                   )
              )
            )
          ), class = "display"
        ) %>% 
                     formatStyle('logFC',
                                 color = styleInterval(c(0), c('blue', 'red'))
                     ) %>%
                     formatStyle(
                       columns = 'Gene_Symbol', 
                       valueColumns = 'logFC',
                       color = styleInterval(c(0), c('blue', 'red')))
    })

# Intersection output
  output$edgerPath <- renderPrint({
    "/transcriptr/workdir/results/edgeR_results"
  })


  output$deseqPath <- renderPrint({
    "/transcriptr/workdir/results/deseq2_results"
  })

    output$limmaPath <- renderPrint({
    "/transcriptr/workdir/results/limma_results"
  })

## Intersect table
interTab <- eventReactive(input$intersectdeg, {
  library(data.table)

  edgeR <- fread(file.path("/transcriptr/workdir/results/edgeR_results", input$showedger))
  DESeq2 <- fread(file.path("/transcriptr/workdir/results/deseq2_results", input$showdeseq))
  limma <- fread(file.path("/transcriptr/workdir/results/limma_results", input$showlimma))

  l <- list(edgeR, DESeq2, limma)
  dtm = Reduce(function(...) merge(..., all = TRUE), l)


  i <- list(
    edgeR = edgeR$id, 
    DESeq2 = DESeq2$id,
    limma = limma$id
  )

  ids <- intersect(intersect(edgeR$id, DESeq2$id), limma$id)
  t <- dtm[id %in% ids ,]

  colnames(t) <- c('id', 'name', 'chr', 'start', 'end', 'strand', 'biotype',
	'edgeR_logFC', 'edgeR_logCPM', 'edgeR_PValue', 'edgeR_FDR', 'edgeR_baseMean',	
	'DESeq2_log2FoldChange', 'DESeq2_lfcSE', 'DESeq2_stat', 'DESeq2_pvalue', 'DESeq2_padj',	
	'limma_logFC', 'limma_AveExpr', 'limma_t', 'limma_P.Value', 'limma_adj.P.Val', 'limma_B')

  t <- t[order(edgeR_FDR, decreasing=FALSE),]

  return(t)
})

#tableIntersect = data.frame(interTab())

output$insettable <- DT::renderDataTable(
  interTab(),
  extensions = 'Buttons',
          escape = FALSE,
          rownames = FALSE,
          filter = "bottom",
          style = "auto",
          selection = "single",
          options = list(
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autowidth = TRUE,
            ordering = TRUE,
            dom = 'Bflrtip',
            scrollX = TRUE,
            lengthMenu = list(c(10, 50,100, -1), c('10', '50', '100', 'All')),
            buttons = list(
              list(extend = "excel", text = "Download Table", title = "", filename = paste0("TranscriptR_","Intersect_Results_", input$showedger, sep=""),
                   exportOptions = list(
                     modifier = list(page = "current")
                   )
              )              
            )
          ), class = "display"
)

### Venn from intersect result
intersectVenn <- eventReactive(input$intersectdeg, {
  library(ggvenn)

  edgeR <- fread(file.path("/transcriptr/workdir/results/edgeR_results", input$showedger))
  DESeq2 <- fread(file.path("/transcriptr/workdir/results/deseq2_results", input$showdeseq))
  limma <- fread(file.path("/transcriptr/workdir/results/limma_results", input$showlimma))

  l <- list(edgeR, DESeq2, limma)
  dtm = Reduce(function(...) merge(..., all = TRUE), l)


  i <- list(
    edgeR = edgeR$id, 
    DESeq2 = DESeq2$id,
    limma = limma$id
  )

  return(i)
})
output$insetplot <- renderPlot(
  ggvenn(
  intersectVenn(), 
  fill_color = c("#0073C2FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
)

### Download Venn
          output$invennDownload <- downloadHandler(
          filename = function() {
            paste("DEGVennPlot",tolower(input$filetype_invenn), sep = ".")},
          content = function(file) {
          if(input$filetype_invenn == "PNG") 
          {png(file, 
            width = input$invenn_pngwidth , 
            height = input$invenn_pngheight, 
            units = input$invenn_pngunits, 
            res = input$invenn_pngresol, 
            type = "cairo")} 
          else if(input$filetype_invenn == "PDF")
          {pdf(file,
            width = input$invenn_pdfwidth , 
            height = input$invenn_pdfheight            
            )}
          else if(input$filetype_invenn == "SVG")
          {svg(file, 
            width = input$invenn_svgwidth , 
            height = input$invenn_svgheight
            )}
          else (tiff(file, 
            width = input$invenn_tiffwidth , 
            height = input$invenn_tiffheight, 
            units = input$invenn_tiffunits, 
            res = input$invenn_tiffresol,
          type = "cairo"))
            

          print(ggvenn(
            intersectVenn(), 
            fill_color = c("#0073C2FF", "#868686FF", "#CD534CFF"),
            stroke_size = 0.5, set_name_size = 4
          ))

          dev.off()
          }
          )

# SessionInfo ouput
output$seninfo <- renderPrint(
  capture.output(sessionInfo())
)

output$geninfo <- renderPrint({
  filepath <- list.files("/transcriptr/workdir/index/gtf/", full.names =T, pattern = ".txt")

  txtfile <- readLines(filepath)
  return(txtfile)
}
)

output$othinfo <- renderPrint({
  filepath <- list.files("/transcriptr/workdir/index/ref/", full.names =T, pattern = ".txt")
  txtfile <- readLines(filepath)
  return(txtfile)
}
)

output$consoleText <- renderPrint({
  message("Getting log, please wait")
    library(readr)
    library(stringr)
    library(dplyr)
    filepath <- list.files("/var/log/shiny-server/", full.names =T, pattern = ".log")
    #filepath <- list.files("/home/jyotirmoy/Downloads/logFile", full.names =T, pattern = ".log")
    txtfile <- readLines(textConnection(filepath[[2]]))
    read_file(txtfile) -> mylog
    tt1 <- unlist(strsplit(mylog, '\n') )
    ix1 <- tt1[grep("steps | failed", tt1)]
    #write(message, ix1, append = TRUE)
    return(ix1)

  })

                  #====================================================#
                    ## Download all transcriptR result as zip ###
                  #====================================================#

                  observeEvent(input$download_btn, {
                    showModal(modalDialog(
                      title = HTML('<span style="color:SteelBlue; font-size: 26px; font-weight:bold; font-family:sans-serif ">Proceed to download results<span>
                                <button type = "button" class="close" data-dismiss="modal" ">
                                </button> '),
                      HTML('<span style="color:LightCoral; font-size: 20px; font-weight:bold; font-family:sans-serif ">NOTE: It is a good idea to close the browser once you finished one analysis and run the next one. The cache memory will be cleared in this way ...<span>
                                <button type = "button" class="close" data-dismiss="modal" ">
                                </button> '),
                      footer = actionButton("ConfirmResult", "Confirm"),

                      easyClose = FALSE
                  ))
                  })

                  observeEvent(input$ConfirmResult, {
                    showModal(
                        modalDialog(
                            HTML('<span style="color:LightSeaGreen; font-size: 26px; font-weight:bold; font-family:sans-serif ">Are you sure?<span>
                              <button type = "button" class="close" data-dismiss="modal" ">
                              </button> '),
                            footer = tagList(
                                downloadButton(outputId = "download_tresult", "Yes"),
                                modalButton("No")
                            ),
                            easyClose = FALSE
                        )
                    )
                  })

                  output$download_tresult <- downloadHandler(
                  filename = function(){
                    paste("TranscriptR_result_", Sys.Date(), ".zip", sep = "")
                  },
                  content = function(file){
                  showNotification("zipping results, please wait...",
                                    duration = 20,
                                    type = "message",
                                    closeButton = TRUE)
                  removeModal()    
                  system('bash -c "zip -r /transcriptr/workdir/results/transcriptr_results.zip /transcriptr/workdir/results/edgeR_results/ /transcriptr/workdir/results/deseq2_results/ /transcriptr/workdir/results/limma_results/ /transcriptr/workdir/results/fastqc/ /transcriptr/workdir/results/logs/ /transcriptr/workdir/results/subread/ /transcriptr/workdir/results/multiqc_report.html /srv/shiny-server/.snakemake/log/*.log"')              
                  file.copy("/transcriptr/workdir/results/transcriptr_results.zip", file)


                  },
                    contentType = "application/zip"    
                  )

}

shinyCatch(transcriptR(), position = "bottom-right", blocking_level = "none", shiny = TRUE, prefix = "transcriptR", trace_back = spsOption("traceback"))


#####============== Transcriptome module end==============#######


#####============== IGviewer module Start==============#######

    if(!dir.exists("tracks"))
  dir.create("tracks")
addResourcePath("tracks", "tracks")

   tbl.bed5 <- reactive({

      inputBed <- input$bedfile
           table.bed <- read.delim(file.path("./transcriptr/workdir/results/bedgraphs/", inputBed), sep = "\t")

     tableBED <- data.frame(
      chr= as.character(table.bed[[1]]),
      start= table.bed[[2]],
      end = table.bed[[3]],
      value = table.bed[[4]],
      stringsAsFactors = FALSE)

      return(tableBED)
      
    
   })

observeEvent(input$submitigv, {


  observeEvent(input$searchButton, {
    printf("--- search")
    searchString = isolate(input$roi)
    if(nchar(searchString) > 0)
      showGenomicRegion(session, id="igvShiny_0", searchString)
  })
  
  observeEvent(input$addBedTrackButton, {
    showGenomicRegion(session, id="igvShiny_0", "chr1:7,426,231-7,453,241")
    loadBedTrack(session, id="igvShiny_0", trackName="bed5", tbl=tbl.bed5());
  })
  
  observeEvent(input$addBed9TrackButton, {
    showGenomicRegion(session, id="igvShiny_0", "chr1:161,199,757-161,201,277")
    loadBedTrack(session, id="igvShiny_0", trackName="bed9", tbl=tbl.bed9)
  })
  
  observeEvent(input$addBedGraphTrackButton, {
    showGenomicRegion(session, id="igvShiny_0", "chr1:7,426,231-7,453,241")
    loadBedGraphTrack(session, id="igvShiny_0", trackName="wig/bedGraph/local", tbl=tbl.bed5(),
                      color="blue", autoscale=TRUE)
  })
  
  observeEvent(input$addAutoscaledGroupBedGraphTrackButton, {
    showGenomicRegion(session, id="igvShiny_0", "chr1:7,432,868-7,433,167")
    loadBedGraphTrack(session, id="igvShiny_0", trackName="wig1a", tbl=tbl.wig, color="blue",
                      autoscale=TRUE, autoscaleGroup=1)
    tbl.wig1b <- tbl.wig
    tbl.wig1b$value <-tbl.wig1b$value * 10
    loadBedGraphTrack(session, id="igvShiny_0", trackName="wig1b", tbl=tbl.wig1b, color="brown",
                      autoscale=TRUE, autoscaleGroup=1)
  })
  
  #
  
  observeEvent(input$addBedGraphTrackFromURLButton, {
    showGenomicRegion(session, id="igvShiny_0", "chr1:154,946,914-155,080,475")
    url <- "https://www.encodeproject.org/files/ENCFF356YES/@@download/ENCFF356YES.bigWig"
    loadBedGraphTrackFromURL(session, id="igvShiny_0", trackName="bedGraph/remote",
                             url=url, color="brown",
                             trackHeight=50, autoscale=TRUE)
  })
  
  
  observeEvent(input$removeUserTracksButton, {
    printf("---- removeUserTracks")
    removeUserAddedTracks(session, id="igvShiny_0")
  })
  
  observeEvent(input$igvReady, {
    printf("--- igvReady")
    containerID <- input$igvReady
    printf("igv ready, %s", containerID)
    loadBedTrack(session, id=containerID, trackName="bed5 loaded on ready", tbl=tbl.bed5(), color="red");
  })
  
  observeEvent(input$trackClick, {
    printf("--- trackclick event")
    x <- input$trackClick
    print(x)
  })
  
  observeEvent(input[["igv-trackClick"]], {
    printf("--- igv-trackClick event")
    x <- input[["igv-trackClick"]]
    print(x)
    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions]
    attribute.values <- as.character(x)[attribute.value.positions]
    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)
    dialogContent <- renderTable(tbl)
    html <- HTML(dialogContent())
    showModal(modalDialog(html, easyClose=TRUE))
  })
  
  observeEvent(input$getChromLocButton, {
    # printf("--- getChromLoc event")
    # sends message to igv.js in browser; currentGenomicRegion.<id> event sent back
    # see below for how that can be captured and displayed
    getGenomicRegion(session, id="igvShiny_0")
  })
  
  observeEvent(input$clearChromLocButton, {
    output$chromLocDisplay <- renderText({" "})
  })
  
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    #printf("new chromLocString: %s", newLoc)
    output$chromLocDisplay <- renderText({newLoc})
  })
  
  # genomes <- c("hg38", "hg19", "mm10", "tair10", "rhos")
  # loci <- c("chr5:88,466,402-89,135,305",  "chr1:7,426,231-7,453,241", "MEF2C", "Mef2c",
  #           "1:7,432,931-7,440,395", "NC_007494.2:370,757-378,078",
  #           "chr1:6,575,383-8,304,088")
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    fasta.file <- file.path("./transcriptr/workdir/index/ref", "ref.fa")
    fastaIndex.file <- file.path("./transcriptr/workdir/index/ref", "ref.fa.fai")
    annotation.file <- file.path("./transcriptr/workdir/index/gff", "ref.gff3")
     
     genomeOptions <- parseAndValidateGenomeSpec(genomeName="local file",
                                                 initialLocus="all",
                                                 stockGenome=FALSE,
                                                 dataMode="localFiles",
                                                 fasta=fasta.file,
                                                 fastaIndex=fastaIndex.file,
                                                 genomeAnnotation=annotation.file
                                                 )
    x <- igvShiny(genomeOptions,
                  displayMode="SQUISHED",
                  tracks=list()
    )
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
})

#####============== IGviewer module end==============#######


#====================================================#
  ## Multi-Dimensional analysis ###
#====================================================#

  mds_size <- reactive({
    return(input$mds_size)
  })

  colorMDS <- reactive({
    groupMDS <- input$file_groupmds
    ext <- tools::file_ext(groupMDS$datapath)
    req(groupMDS)

    MDSgroup <- read.delim(groupMDS$datapath, header = T)

    if (input$ncolourmds == 2) {
      colr = c(input$mdscol21, input$mdscol22)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 3){
      colr = c(input$mdscol31, input$mdscol32, input$mdscol33)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 4){
      colr = c(input$mdscol41, input$mdscol42, input$mdscol43, input$mdscol44)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 5){
      colr = c(input$mdscol51, input$mdscol52, input$mdscol53, input$mdscol54, input$mdscol55)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 6){
      colr = c(input$mdscol61, input$mdscol62, input$mdscol63, input$mdscol64, input$mdscol65, input$mdscol66)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 7){
      colr = c(input$mdscol71, input$mdscol72, input$mdscol73, input$mdscol74, input$mdscol75, 
               input$mdscol76, input$mdscol77)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 8){
      colr = c(input$mdscol81, input$mdscol82, input$mdscol83, input$mdscol84, input$mdscol85, 
               input$mdscol86, input$mdscol87, input$mdscol88)[factor(MDSgroup$group)]
    } else if (input$ncolourmds == 9){
      colr = c(input$mdscol91, input$mdscol92, input$mdscol93, input$mdscol94, input$mdscol95, 
               input$mdscol96, input$mdscol97, input$mdscol98, input$mdscol99)[factor(MDSgroup$group)]
    } else {
      colr = c(input$mdscol101, input$mdscol102, input$mdscol103, input$mdscol104, input$mdscol105, 
               input$mdscol106, input$mdscol107, input$mdscol108, input$mdscol109, input$mdscol1010)[factor(MDSgroup$group)]
    }
  })
  
  #Input file view
  mdsFile = reactive({
    fileMDS <- input$file_mds
    ext <- tools::file_ext(fileMDS$datapath)
    req(fileMDS)

    fileMDSpath <- fileMDS$datapath
    MDSdata <- switch(ext,
               txt = readr::read_table(fileMDSpath),
               tsv = readr::read_tsv(fileMDSpath),
               csv = readr::read_csv(fileMDSpath),
               xls = readxl::read_xls(fileMDSpath),
               xlsx = readxl::read_xlsx(fileMDSpath, col_names = TRUE))
    
    #MDSdata2 <- as.data.frame(MDSdata)
    #rownames(MDSdata2) <- MDSdata[,1]
    return(MDSdata[,-1])
  })  

  output$mdatest1 <- DT::renderDataTable(
    mdsFile(), 
    #rownames = TRUE, 
    editable = TRUE,
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    ) 
  )

  mdsFileGroup = reactive({
    groupMDS <- input$file_groupmds
    ext <- tools::file_ext(groupMDS$datapath)
    req(groupMDS)

    groupMDSpath <- groupMDS$datapath
    MDSgroup <- switch(ext,
               txt = readr::read_table(groupMDSpath),
               tsv = readr::read_tsv(groupMDSpath),
               csv = readr::read_csv(groupMDSpath),
               xls = readxl::read_xls(groupMDSpath),
               xlsx = readxl::read_xlsx(groupMDSpath))
    return(MDSgroup)
  })  

  output$mdatest2 <- DT::renderDataTable(
    mdsFileGroup(), rownames = FALSE, editable = TRUE,
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    ) 
  )

  mdsPlot1 <- eventReactive(input$runmda, {
    id_mds <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_mds), add = TRUE)

    fileMDS <- input$file_mds
    ext <- tools::file_ext(fileMDS$datapath)
    req(fileMDS)

    # fileMDSpath <- fileMDS$datapath
    # MDSdata <- switch(ext,
    #            txt = readr::read_table(fileMDSpath),
    #            tsv = readr::read_tsv(fileMDSpath),
    #            csv = readr::read_csv(fileMDSpath),
    #            xls = readxl::read_xls(fileMDSpath),
    #            xlsx = readxl::read_xlsx(fileMDSpath, col_names = TRUE))    

    # #rownames(MDSdata) <- MDSdata[,1]
    MDSdata <- read.delim(fileMDS$datapath, sep = "\t", row.names = 1)
    samplenames <- colnames(MDSdata) # sample names added 

    groupMDS <- input$file_groupmds
    ext <- tools::file_ext(groupMDS$datapath)
    req(groupMDS)

    groupMDSpath <- groupMDS$datapath
    MDSgroup <- switch(ext,
               txt = readr::read_table(groupMDSpath),
               tsv = readr::read_tsv(groupMDSpath),
               csv = readr::read_csv(groupMDSpath),
               xls = readxl::read_xls(groupMDSpath),
               xlsx = readxl::read_xlsx(groupMDSpath))
    # #groups <- colnames(MDSdata) # sample names added 
    
    MDSdata2 <- as.matrix(MDSdata)    
    
    figureMDS <- plotMDS(MDSdata2,
              top = input$numdms,
              #labels = samplenames,
              pch = 19,
              cex = as.numeric(input$cexmds), 
              dim.plot = c(input$dim1mds,input$dim2mds),               
              #ndim = max(input$dimmds),
              #ndim = max(dim.plot), 
              #gene.selection = "pairwise", 
              xlab = input$xlabmds, 
              ylab = input$ylabmds,
              #col = c(input$mdscol1, input$mdscol2)[factor(MDSgroup$group)],
              col = colorMDS()
              
      )
    
    figureMDS
  })

  output$mdsPlot <- renderPlot({    
    mdsPlot1()
  },
    width = mds_size,
    height = mds_size,
    outputArgs = list()
  )

  output$mdsDownload <- downloadHandler(
    filename = function(){
      paste("MDSplot", tolower(input$filetype_mds), sep =".")
    },
    content = function(file)
    {
      width  <- mds_size()
      height <- mds_size()
      #width  <- session$clientData$output_plot_width
      #height <- ((session$clientData$output_plot_height)*1)
      #pixelratio <- session$clientData$pixelratio
      pixelratio <- 2

      if(input$filetype_mds == "PNG")
        png(file, width=4, height=4, units = "in", res=300)
      else if(input$filetype_mds == "SVG")
        svg(file, width=8, height=8)
      else if(input$filetype_mds == "TIFF")
        tiff(file, width=4, height=4, units = "in", res=300)
      else
        pdf(file, width = 8, height = 8)

      mdsPlot1()    

      dev.off()
    }
  )

#####============== PCA module Start ==============#######

  pcaPlot1 <- eventReactive(input$runpca,
  {
    #library(explor)
    id_pca <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_pca), add = TRUE)

    fileMDS <- input$file_pca
    ext <- tools::file_ext(fileMDS$datapath)
    req(fileMDS)

    fileMDSpath <- fileMDS$datapath
    # MDSdata <- switch(ext,
    #            txt = readr::read_table(fileMDSpath),
    #            tsv = readr::read_tsv(fileMDSpath),
    #            csv = readr::read_csv(fileMDSpath),
    #            xls = readxl::read_xls(fileMDSpath),
    #            xlsx = readxl::read_xlsx(fileMDSpath))  

    MDSdata <- read.delim(fileMDS$datapath, sep = "\t", row.names = 1)

    MDSdata1 <- MDSdata[1:input$numpca,]
    data.pca <- PCA(t(MDSdata1), graph = FALSE)
    eig.val <- get_eigenvalue(data.pca)
    ind <- fviz_pca_ind(data.pca)
    var <- fviz_pca_var(data.pca)

    fileGroup <- input$group_pca
    ext <- tools::file_ext(fileGroup$datapath)
    req(fileGroup)

    groupData <- read.delim(fileGroup$datapath, sep = "\t",header = T) 
    group <- colnames(groupData)
    
    pc <- fviz_pca_ind(data.pca,
                 fill.ind = groupData$group,
                 pch = 21, 
                 pointsize = input$cexpca,
                 palette = input$colorpca, #"jco",
                 title = input$maintitlepca, #"Principal Component Analysis",
                 xlab = input$xlabpca, 
                 ylab = input$ylabpca,
                 legend.title = input$legendtitlepca, #"Groups",
                 mean.point = input$meanpca, #FALSE,
                 repel = TRUE) %>% config(displaylogo = FALSE)

  }
  )

  output$pcaplot <- renderPlotly({
    pcaPlot1() %>% config(displaylogo = FALSE)
  })
#####============== PCA module end==============#######


#####============== Venn module Start ==============#######
  #====================================================#
  ## Venn module ####
#====================================================#
  venn_type <- reactive({
    return(input$venn_type)
  })

  doWeights <- reactive({
    return(input$doWeights)
  })

  doEuler <- reactive({
    return(input$doEuler)
  })

  venn_size <- reactive({
    return(input$venn_size)
  })

  venn_lwd <- reactive({
    return(as.numeric(input$venn_lwd))
  })

  venn_labelsize <- reactive({
    return(as.numeric(input$venn_labelsize))
  })

  venn_color_type <- reactive({
    return(input$venn_color_type)
  })


  venn_cex <- reactive({
    return(as.numeric(input$venn_cex))
  })

  venn_lty <- reactive({
    return(as.numeric(input$venn_lty))
  })

  set1_color <- reactive({
    return(input$set1_color)
  })
  set2_color <- reactive({
    return(input$set2_color)
  })
  set3_color <- reactive({
    return(input$set3_color)
  })
  set4_color <- reactive({
    return(input$set4_color)
  })
  set5_color <- reactive({
    return(input$set5_color)
  })
  set6_color <- reactive({
    return(input$set6_color)
  })



  venn_data <- reactive({
    inFile <- input$file_venn
    string <- ""

    if(is.null(inFile) == F)
    {
      data <- read_delim(input$file_venn$datapath, input$sep_venn , escape_double = FALSE, trim_ws = TRUE, col_names = input$header_venn)
      data <- lapply(data, function(x) x[!is.na(x)])
      return(data)
    }else{
      if (string != "")
      {
        string <- gsub("\n", "", string)
        if(string != ""){
          string <- as.list(unlist(strsplit(string, ",")))
          names <- lapply(string, function(x){x <- unlist(strsplit(x, "=")); x <- x[1]})
          names <- unlist(lapply(names, function(x){x <- gsub(" ", "", x)}))
          vennData <- as.numeric(unlist(lapply(string, function(x){x <- unlist(strsplit(x,"=")); x <- x[2]})))
          names(vennData) <- names
          return(vennData)
        }else{
          return(NULL)
        }
      }
      else{
        data <- read_delim(file = file.path('https://raw.githubusercontent.com/asntech/intervene-shiny/master/data/Whyte_et_al_2013_SEs_genes.csv'), ",", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
        return(lapply(data, function(x) x[!is.na(x)]))
      }
    }
  })

  set_names <- reactive({    
    names <- names(venn_data())
    return(names)
  })

  output$venn_sets <- renderUI({
    venn_sets <- selectInput('venn_sets', label = "Select sets",
                             choices = as.character(set_names()),
                             multiple = T, selectize = T, selected = as.character(set_names()[1:5]))
    return(venn_sets)
  })

  venn_selected_names <- reactive({
    venn_selected_names <- as.character(c(input$venn_sets))
  })

  venn_data_filtered <- reactive({

    data <- venn_data()
    if(is.null(input$venn_sets)){

    }else{
      data <- data[c(venn_selected_names())]
      return(data)
    }
    return(data)
  })

  venn_combinations <- reactive({
    string <- ""
    data <- venn_data_filtered()
    if (string !=""){
      return(data)
    }else
    {
      return(Venn(data))
    }
  })

  get_venn_gp <- reactive({
    venn_gp <- VennThemes(compute.Venn(venn_combinations()))
    venn_gp$SetText <- lapply(venn_gp$SetText,function(x) {x$fontsize<-venn_labelsize(); return(x)})
    venn_gp$FaceText <- lapply(venn_gp$FaceText,function(x) {x$cex<-venn_cex(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lwd<-venn_lwd(); return(x)})
    venn_gp$Set <- lapply(venn_gp$Set,function(x) {x$lty<-venn_lty(); return(x)})

    if (venn_color_type () == 'custom'){
      venn_gp$Set$Set1$col <- set1_color()
      venn_gp$Set$Set2$col <- set2_color()
      venn_gp$Set$Set3$col <- set3_color()
      venn_gp$Set$Set4$col <- set4_color()
      venn_gp$Set$Set5$col <- set5_color()
      venn_gp$Set$Set6$col <- set6_color()

      venn_gp$SetText$Set1$col <- set1_color()
      venn_gp$SetText$Set2$col <- set2_color()
      venn_gp$SetText$Set3$col <- set3_color()
      venn_gp$SetText$Set4$col <- set4_color()
      venn_gp$SetText$Set5$col <- set5_color()
      venn_gp$SetText$Set6$col <- set6_color()

    }

    return(venn_gp)
  })

  data_size <- reactive({
    return(length(venn_data_filtered()))
  })

  get_venn_type <- reactive({
    if (venn_type() == 'Classical'){
      if(data_size() < 4)
        return("circles")
      else
        return("ellipses")
    }else if (venn_type() == 'ChowRuskey' && data_size() < 3){
      return("circles")
    }
    else{
      return(venn_type())
    }
  })

  output$vennPlot <- renderPlot({
    plot(compute.Venn(venn_combinations(), doWeights = doWeights(), doEuler = doEuler(), type = get_venn_type()),
         gp = get_venn_gp(),
         show = list(Universe = FALSE)
    )
  },
    width = venn_size,
    height = venn_size,
    outputArgs = list()
  )

  output$VennDown <- downloadHandler(
    filename = function(){
      paste("Venn_diagram", tolower(input$filetype_venn), sep =".")
    },
    content = function(file){
      width  <- venn_size()
      height <- venn_size()
      pixelratio <- 2

      if(input$filetype_venn == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio, units = "px", res=72*pixelratio)
      else if(input$filetype_venn == "SVG")
        svg(file, width=8, height=8)
      else if(input$filetype_venn == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = 8, height = 8)

      plot(venn_combinations(),
           doWeights = doWeights(),
           type = get_venn_type(),
           doEuler = doEuler(),
           show = list(Universe = FALSE)           
      )
      dev.off()
    }
  )

#====================================================#
  ## UpSet module ####
#====================================================#
  #Some of the code for upset module is taken from
  #https://github.com/hms-dbmi/UpSetR-shiny

  output$plot_text <- renderUI({
    if(is.null(My_data()) == T){
      h5("There is no data entered. Please upload your data to draw UpSet plot here!")
    }
    else{
      HTML(" ")
    }
  })

  My_dat <- reactive({
    inFile <- input$file1
    input_type = input$upset_input_type
    if (is.null(inFile) == T){

      My_dat<- fromExpression(c('H3K4me2&H3K4me3'=16321,'H3K4me2&H3K4me3&H3K27me3'=5756,'H3K27me3'=25174,'H3K4me3&H3K27me3'=15539,'H3K4me3'=32964,'H3K4me2&H3K27me3'=19039,'H3K4me2'=60299,'H3K27ac&H3K4me2&H3K4me3&H3K27me3'=7235,'H3K27ac&H3K4me2&H3K4me3'=17505,'H3K27ac&H3K4me2'=21347,'H3K27ac&H3K4me2&H3K27me3'=1698,'H3K27ac&H3K4me3'=8134,'H3K27ac&H3K4me3&H3K27me3'=295,'H3K27ac&H3K27me3'=7605,'H3K27ac'=42164))
      return(My_dat)
    }
    else if(is.null(inFile) == F && input_type == 'binary'){
      read.csv(inFile$datapath, header = input$header,
               sep = input$sep, quote = input$quote)
    }else if (is.null(inFile) == F && input_type == 'list'){
      My_dat <- read_delim(inFile$datapath, input$sep , escape_double = FALSE, trim_ws = TRUE, col_names = input$header)
      My_dat <- fromList(lapply(as.list(My_dat), function(x) x[!is.na(x)]))

      return(My_dat)
    }else{
      return(NULL)
    }
  })

  venneulerData <- reactive({
    string <- input$upset_comb
    string <- gsub("\n", "", string)
    if(string != ""){
      string <- as.list(unlist(strsplit(string, ",")))
      names <- lapply(string, function(x){x <- unlist(strsplit(x, "=")); x <- x[1]})
      names <- unlist(lapply(names, function(x){x <- gsub(" ", "", x)}))
      values <- as.numeric(unlist(lapply(string, function(x){x <- unlist(strsplit(x,"=")); x <- x[2]})))
      names(values) <- names
      venneuler <- fromExpression(values)
      return(venneuler)
    }
  })

  My_data <- reactive({
    string <- input$upset_comb
    if(string != ""){
      My_data <- venneulerData()
    }
    else {
      My_data <- My_dat()
    }
    return(My_data)
  })

  FindStartEnd <- function(data){
    startend <- c()
    for(i in 1:ncol(data)){
      column <- data[, i]
      column <- (levels(factor(column)))
      if((column[1] == "0") && (column[2] == "1" && (length(column) == 2))){
        startend[1] <- i
        break
      }
      else{
        next
      }
    }
    for(i in ncol(data):1){
      column <- data[ ,i]
      column <- (levels(factor(column)))
      if((column[1] == "0") && (column[2] == "1") && (length(column) == 2)){
        startend[2] <- i
        break
      }
      else{
        next
      }
    }
    return(startend)
  }

  startEnd <- reactive({
    startEnd <- FindStartEnd(My_data())
  })

  setSizes <- reactive({
    if(is.null(My_data()) != T){
      sizes <- colSums(My_data()[startEnd()[1]:startEnd()[2]])
      sizes <- sizes[order(sizes, decreasing = T)]

      names <- names(sizes); sizes <- as.numeric(sizes);
      maxchar <- max(nchar(names))
      total <- list()
      for(i in 1:length(names)){
        spaces <- as.integer((maxchar - nchar(names[i]))+1)
        spaces <- paste(rep(" ", each=spaces), collapse = "")
        total[[i]] <- paste(paste(names[i], ":", sep=""), spaces, sizes[i], "\n", sep="")
      }
      total <- unlist(total)
      total <- paste(total, collapse = " ")
      return(total)
    }
    else{
      return(NULL)
    }
  })

  output$setsizes <- renderText({
    if(is.null(setSizes()) != T){
      paste("---Set Sizes---\n", setSizes())
    }
    else{
      paste("---Set Sizes---\n", "\n No Data Entered")
    }
  })

  Specific_sets <- reactive({
    Specific_sets <- as.character(c(input$upset_sets))
  })

  output$sets <- renderUI({
    if(is.null(My_data()) == T){
      sets <-  selectInput('upset_sets', label="Select at least two sets ",
                           choices = NULL,
                           multiple=TRUE, selectize=TRUE, selected = Specific_sets())
    }
    else{
      data <- My_data()[startEnd()[1]:startEnd()[2]]
      topfive <- colSums(data)
      topfive <- as.character(head(names(topfive[order(topfive, decreasing = T)]), 5))
      sets <- selectInput('upset_sets', label="Select sets ",
                          choices = as.character(colnames(My_data()[ , startEnd()[1]:startEnd()[2]])),
                          multiple=TRUE, selectize=TRUE, selected = topfive)
    }
    return(sets)
  })


  mat_prop <- reactive({
    mat_prop <- input$mbratio
  })
  upset_width <- reactive({
    return(input$upset_width)
  })
  upset_height <- reactive({
    return(input$upset_height)
  })

  bar_prop <- reactive({
    bar_prop <- (1 - input$mbratio)
  })

  orderdat <- reactive({
    orderdat <- as.character(input$order)
    if(orderdat == "degree"){
      orderdat <- c("degree")
    }
    else if(orderdat == "freq"){
      orderdat <- "freq"
    }
    return(orderdat)
  })

  show_numbers <- reactive({
    show_numbers <- input$show_numbers
    if(show_numbers){
      show_numbers <- "yes"
      return(show_numbers)
    }
    else{
      show_numbers <- FALSE
      return(show_numbers)
    }

  })

  main_bar_color <- reactive({
    mbcolor <- input$mbcolor
    return(mbcolor)
  })
  sets_bar_color <- reactive({
    sbcolor <- input$sbcolor
    return(sbcolor)
  })


  decrease <- reactive({
    decrease <- as.character(input$decreasing)
    if(decrease == "inc"){
      decrease <- FALSE
    }
    else if(decrease == "dec"){
      decrease <- TRUE
    }
    return(decrease)
  })

  number_angle <- reactive({
    angle <- input$angle
    return(angle)
  })

  line_size <- reactive({
    line_size <- input$linesize
    return(line_size)
  })

  emptyIntersects <- reactive({
    if(isTRUE(input$empty)){choice <- "on"
      return(choice)
    }
    else{
      return(NULL)
    }
  })

  scale.intersections <- reactive({
    return(input$scale.intersections)
  })

  scale.sets <- reactive({
    return(input$scale.sets)
  })

  keep.order <- reactive({
    return(input$keep.order)
  })

  # A plot of fixed size
  output$plot1 <- renderPlot({

    if(length(My_data()) == 0){stop()}
    if(length(Specific_sets()) == 1){
      stop()
    }
    upset(data = My_data(),
          nintersects = input$nintersections,
          point.size = input$pointsize,
          line.size = line_size(),
          sets = Specific_sets(),
          order.by = orderdat(),
          main.bar.color= main_bar_color(),
          sets.bar.color= sets_bar_color(),
          decreasing = c(decrease()),
          show.numbers = show_numbers(),
          number.angles = number_angle(),
          scale.intersections = scale.intersections(),
          scale.sets = scale.sets(),
          keep.order = keep.order(),
          mb.ratio = c(as.double(bar_prop()), as.double(mat_prop())),
          empty.intersections = emptyIntersects(),
          text.scale = c(input$intersection_title_scale, input$intersection_ticks_scale,
                         input$set_title_scale, input$set_ticks_scale, input$names_scale,
                         input$intersection_size_numbers_scale))},
    width = upset_width,
    height = upset_height
  )

  output$UpSetDown <- downloadHandler(

    filename = function(){
      paste("UpSet_plot", tolower(input$filetype), sep =".")
    },
    content = function(file){
      width <- upset_width()
      height <- upset_height()
      pixelratio <- 2

      if(input$filetype == "PNG")
        png(file, width=width*pixelratio, height=height*pixelratio,
            res=72*pixelratio, units = "px")
      else if(input$filetype == "SVG")
        svg(file, width = width/100, height = height/100)
      else if(input$filetype == "TIFF")
        tiff(file, width=width*pixelratio, height=height*pixelratio, units = "px")
      else
        pdf(file, width = width/100, height = height/100, onefile=FALSE)

      upset(data = My_data(),
            nintersects = input$nintersections,
            point.size = input$pointsize,
            line.size = line_size(),

            sets = Specific_sets(),
            order.by = orderdat(),
            main.bar.color= main_bar_color(),
            sets.bar.color= sets_bar_color(),
            decreasing = c(decrease()),
            number.angles = number_angle(),
            show.numbers = show_numbers(),
            scale.intersections = scale.intersections(),
            scale.sets = scale.sets(),
            keep.order = keep.order(),
            mb.ratio = c(as.double(bar_prop()), as.double(mat_prop())),
            empty.intersections = emptyIntersects(),
            text.scale = c(input$intersection_title_scale, input$intersection_ticks_scale,
                           input$set_title_scale, input$set_ticks_scale, input$names_scale,
                           input$intersection_size_numbers_scale))

      dev.off()
    }
  )
#####============== UpSet module end==============#######

####============= Pathway analysis module start ====#######
pathwayAnalysis <- function(){

  filePath <- reactive({
    filepath1 <- input$file7
    ext <- tools::file_ext(filepath1$datapath)
    req(filepath1)    

    filepathway <- filepath1$datapath
    filepath1 <- switch(ext,
               txt = readr::read_table(filepathway),
               tsv = readr::read_tsv(filepathway),
               csv = readr::read_csv(filepathway),
               xls = readxl::read_xls(filepathway),
               xlsx = readxl::read_xlsx(filepathway))

    # filepath1 <- read.delim(filepath1$datapath,
    #   sep = "\t",
    #   header = TRUE
    # )

    updateSelectInput(session, "logfc_colm", choices = names(filepath1))
    updateSelectInput(session, "genesymbol_colm", choices = names(filepath1))

    return(filepath1)
  })

  output$pathwaydatafile <- DT::renderDataTable(
    filePath(),
    extensions = "Buttons",
          escape = FALSE,
          rownames = FALSE,
          filter = "bottom",
          style = "auto",
          selection = "single",
          options = list(
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autowidth = TRUE,
            ordering = TRUE,
            dom = "Bflrtip",
            scrollX = TRUE,
            lengthMenu = list(c(6, 20, 100, -1), c("6", "20", "100", "All"))
    )
  )

 


  pathwayDataTable <- eventReactive(input$Run,{

    id_pathplot1 <- showNotification("analysing, please wait...", duration = NULL, closeButton = FALSE)
     on.exit(removeNotification(id_pathplot1), add = TRUE)

    pathwaydata <- input$file7
    ext <- tools::file_ext(pathwaydata$datapath)
    
    req(pathwaydata)
    #validate(need(ext == "csv", "Please upload a csv file"))
    
    filepathway <- pathwaydata$datapath
    pathdata <- switch(ext,
               txt = readr::read_table(filepathway),
               tsv = readr::read_tsv(filepathway),
               csv = readr::read_csv(filepathway),
               xls = readxl::read_xls(filepathway),
               xlsx = readxl::read_xlsx(filepathway))

    #pathdata <- read.delim(pathwaydata$datapath, sep = "\t")
    pathdata1 <- data.frame(pathdata[[input$genesymbol_colm]], pathdata[[input$logfc_colm]])
    colnames(pathdata1) <- c("genesymbol", "logFC")

            # purpose: 1) to get the entrezID, use mygene library. 2) mygene library needs taxa ID, 3) taxa ID matched with KEGG database and saved in www folder

    OrgDbReactome <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgReac == 'human'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$pathOrgReac == 'mouse') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgReac == 'rat') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgReac == 'celegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgReac == 'yeast') {
        namesdb <- print("org.Sc.sgd.db")}
      else if (input$pathOrgReac == 'zebrafish') {
        namesdb <- print("org.Dr.eg.db")}
      else {
        namesdb <- print("org.Dm.eg.db")}
      
      return(namesdb)
    #})
    }

    OrgDbWiki <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgWiki == 'Homo sapiens'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$ppathOrgWiki == 'Mus musculus') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgWiki == 'Rattus norvegicus') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgWiki == 'Arabidopsis thaliana') {
        namesdb <- print("org.At.tair.db")}
      else if (input$pathOrgWiki == 'Bos taurus') {
        namesdb <- print("org.Bt.eg.db")}
      else if (input$pathOrgWiki == 'Caenorhabditis elegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgWiki == 'Canis familiaris') {
        namesdb <- print("org.Cf.eg.db")}
      else if (input$pathOrgWiki == 'Danio rerio') {
        namesdb <- print("org.Dr.eg.db")}
      else if (input$pathOrgWiki == 'Drosophila melanogaster') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Gallus gallus') {
        namesdb <- print("org.Gg.eg.db")}
      else if (input$pathOrgWiki == 'Pan troglodytes') {
        namesdb <- print("org.Pt.eg.db")}
      else if (input$pathOrgWiki == 'Populus trichocarpa') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Saccharomyces cerevisiae') {
        namesdb <- print("org.Sc.sgd.db")}
      else {
        namesdb <- print("org.Ss.eg.db")}
      return(namesdb)
    #})
    }

    gene <- function() {
        if(input$pathDB == '1') {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbReactome())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")          

          return(geneList)

        } else if (input$pathDB == '2') {
          taxaData <- reactive({
            dbTable <- read.csv("www/kegg_organisms_taxonomyID.csv", header = F, sep = ",")
            matchID <- dbTable[dbTable$V2 == input$pathOrgKegg, ]
            return(matchID)
          }) #taxa ID save in taxaData()$V3

          library(mygene)
          geneIDconv <- queryMany(pathdata1$genesymbol, scopes="symbol", fields="entrezgene", 
                                  species=taxaData()$V3)
          
          new2 = data.frame(
              SYMBOL = geneIDconv$query,
              ENTREZID = geneIDconv$entrezgene
            )
          geneIDconvTable <- new2
          
          colnames(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconvTable, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          return(geneList)
        } else {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbWiki())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")
          
          return(geneList)
        }
    }         


    shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }


    if(input$pathType == "1") {
      if(input$pathDB == "1") {
        message("starting Reactome analysis")
        genelists <- names(gene())

        pathdata.reactome.ora.run <- enrichPathway(gene = genelists,
                                                   pvalueCutoff=input$adjPvalue,
                                                   pAdjustMethod = input$pvaladj,
                                                   readable = T)
        message("Finished Reactome ORA run, creating Figure and Table")
        pathdata.reactome.ora1 <- as.data.frame(pathdata.reactome.ora.run)
        
        tmp1 <- pathdata.reactome.ora1 %>%
          separate(GeneRatio, c("GeneR1", "GeneR2"), "/")
        tmp1$GeneRatio <- as.numeric(tmp1$GeneR1)/as.numeric(tmp1$GeneR2)
        pathdata.reactome.ora1$GeneRatio <- tmp1$GeneRatio
        tmp2 <- pathdata.reactome.ora1 %>%
          separate(BgRatio, c("bg1", "bg2"), "/")
        tmp2$BgRatio <- as.numeric(tmp2$bg1)/as.numeric(tmp2$bg2)
        pathdata.reactome.ora1$BgRatio <- tmp2$BgRatio
        
        
        tableReactomeORA  = data.frame(
          ID = sapply(pathdata.reactome.ora1$ID, function(x)
            toString(tags$a(href=paste0("https://reactome.org/content/query?q=", x),target = "_blank", x))),
          Description = pathdata.reactome.ora1$Description,
          GeneRatio = pathdata.reactome.ora1$GeneRatio,
          BgRatio = pathdata.reactome.ora1$BgRatio,
          pvalue = pathdata.reactome.ora1$pvalue,
          p.adjust = pathdata.reactome.ora1$p.adjust,
          qval = pathdata.reactome.ora1$qvalue,          
          geneID = pathdata.reactome.ora1$geneID,
          Count = pathdata.reactome.ora1$Count) 
        
      } else if(input$pathDB == "2"){
        message("Starting KEGG ORA")

        genelists <- names(gene())
        pathdata.kegg.ora.run <- enrichKEGG(gene = genelists,
                                            organism = input$pathOrgKegg,
                                            pvalueCutoff=input$adjPvalue,
                                            pAdjustMethod = input$pvaladj
                                            #readable = T
        )

        message("Finished KEGG ORA run, creating Figure and Table")

        pathdata.kegg.ora1 <- as.data.frame(pathdata.kegg.ora.run)
        
        tmp1 <- pathdata.kegg.ora1 %>%
          separate(GeneRatio, c("GeneR1", "GeneR2"), "/")
        tmp1$GeneRatio <- as.numeric(tmp1$GeneR1)/as.numeric(tmp1$GeneR2)
        pathdata.kegg.ora1$GeneRatio <- tmp1$GeneRatio
        tmp2 <- pathdata.kegg.ora1 %>%
          separate(BgRatio, c("bg1", "bg2"), "/")
        tmp2$BgRatio <- as.numeric(tmp2$bg1)/as.numeric(tmp2$bg2)
        pathdata.kegg.ora1$BgRatio <- tmp2$BgRatio

        filepath_rows <- nrow(pathdata.kegg.ora1)

        pathdata.kegg.ora1$pathmap <- shinyInput(actionButton, filepath_rows, "button_", label = "Show map", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})')
        
        tablekeggORA  = data.frame(
          KEGGID = pathdata.kegg.ora1$ID,
          PathMap = pathdata.kegg.ora1$pathmap,
          ID = sapply(pathdata.kegg.ora1$ID, function(x)
            toString(tags$a(href=paste0("https://www.kegg.jp/pathway/", x),target = "_blank", x))),
          Description = pathdata.kegg.ora1$Description,
          GeneRatio = pathdata.kegg.ora1$GeneRatio,
          BgRatio = pathdata.kegg.ora1$BgRatio,
          pvalue = pathdata.kegg.ora1$pvalue,
          p.adjust = pathdata.kegg.ora1$p.adjust,
          qval = pathdata.kegg.ora1$qvalue,          
          geneID = pathdata.kegg.ora1$geneID,
          Count = pathdata.kegg.ora1$Count)

          
        
      } else {
        message("Starting WikiPathways ORA")

        genelists <- names(gene())
        pathdata.wiki.ora.run <- enrichWP(gene = genelists,
                                          organism = input$pathOrgWiki,
                                          pvalueCutoff = input$adjPvalue,
                                          pAdjustMethod = input$pvaladj
        )

        message("Finished WikiPathways ORA run, creating Figure and Table")

        pathdata.wiki.ora1 <- as.data.frame(pathdata.wiki.ora.run)
        
        tmp1 <- pathdata.wiki.ora1 %>%
          separate(GeneRatio, c("GeneR1", "GeneR2"), "/")
        tmp1$GeneRatio <- as.numeric(tmp1$GeneR1)/as.numeric(tmp1$GeneR2)
        pathdata.wiki.ora1$GeneRatio <- tmp1$GeneRatio
        tmp2 <- pathdata.wiki.ora1 %>%
          separate(BgRatio, c("bg1", "bg2"), "/")
        tmp2$BgRatio <- as.numeric(tmp2$bg1)/as.numeric(tmp2$bg2)
        pathdata.wiki.ora1$BgRatio <- tmp2$BgRatio
        
        tablewikiORA  = data.frame(
          #WikiID = pathdata.wiki.ora1$ID,
          ID = sapply(pathdata.wiki.ora1$ID, function(x)
            toString(tags$a(href=paste0("https://www.wikipathways.org//index.php?query=", x, "&title=Special%3ASearchPathways&doSearch=1&sa=Search"),target = "_blank", x))),
          Description = pathdata.wiki.ora1$Description,
          GeneRatio = pathdata.wiki.ora1$GeneRatio,
          BgRatio = pathdata.wiki.ora1$BgRatio,
          pvalue = pathdata.wiki.ora1$pvalue,
          p.adjust = pathdata.wiki.ora1$p.adjust,
          qval = pathdata.wiki.ora1$qvalue,
          geneID = pathdata.wiki.ora1$geneID,
          Count = pathdata.wiki.ora1$Count)        
      }
    } else {
      if(input$pathDB == "1") {
        geneList <- gene()
        pathdata.reactome.gsea.run <- gsePathway(geneList,
                                                 pvalueCutoff = input$adjPvalue,
                                                 organism = "human",
                                                 exponent = 1,
                                                 minGSSize = input$mingssize,
                                                 maxGSSize = input$maxgssize,
                                                 eps = 0,
                                                 pAdjustMethod = input$pvaladj,
                                                 verbose = FALSE,
                                                 by = "fgsea")
        # GSEA table
        pathdata.reactome.gsea.tmp <- as.data.frame(pathdata.reactome.gsea.run)
        
        tableReactomeGSEA  = data.frame(
          #ReactomeID = pathdata.reactome.gsea.tmp$ID,
          ReactomeID = sapply(pathdata.reactome.gsea.tmp$ID, function(x)
            toString(tags$a(href=paste0("https://reactome.org/content/query?q=", x),target = "_blank", x))),
          Description = pathdata.reactome.gsea.tmp$Description,
          setSize = pathdata.reactome.gsea.tmp$setSize,
          enrichmentScore = pathdata.reactome.gsea.tmp$enrichmentScore,
          pvalue = pathdata.reactome.gsea.tmp$pvalue,
          p.adjust = pathdata.reactome.gsea.tmp$p.adjust,
          #qval = pathdata.reactome.gsea.tmp$qvalues,
          rank = pathdata.reactome.gsea.tmp$rank,
          leading_edge = pathdata.reactome.gsea.tmp$leading_edge,
          core_enrichment = pathdata.reactome.gsea.tmp$core_enrichment
        )
              
      } else if(input$pathDB == "2"){
        geneList <- gene()
        pathdata.kegg.gsea.run <- gseKEGG(geneList,
                                          pvalueCutoff = input$adjPvalue,
                                          organism = input$pathOrgKegg,
                                          exponent = 1,
                                          minGSSize = input$mingssize,
                                          maxGSSize = input$maxgssize,
                                          eps = 0,
                                          pAdjustMethod = input$pvaladj,
                                          verbose = FALSE,
                                          by = "fgsea")
        # GSEA table
        pathdata.kegg.gsea.tmp <- as.data.frame(pathdata.kegg.gsea.run)

        filepath_rows <- nrow(pathdata.kegg.gsea.tmp)

        pathdata.kegg.gsea.tmp$pathmap <- shinyInput(actionButton, filepath_rows, "button_", label = "Show map", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})')
        
        tablekeggGSEA  = data.frame(
          KEGGID = pathdata.kegg.gsea.tmp$ID,
          PathMap = pathdata.kegg.gsea.tmp$pathmap,
          ID = sapply(pathdata.kegg.gsea.tmp$ID, function(x)
            toString(tags$a(href=paste0("https://www.kegg.jp/pathway/", x),target = "_blank", x))),
          Description = pathdata.kegg.gsea.tmp$Description,
          setSize = pathdata.kegg.gsea.tmp$setSize,
          enrichmentScore = pathdata.kegg.gsea.tmp$enrichmentScore,
          NES = pathdata.kegg.gsea.tmp$NES,
          pvalue = pathdata.kegg.gsea.tmp$pvalue,
          p.adjust = pathdata.kegg.gsea.tmp$p.adjust,
          #qval = pathdata.kegg.gsea.tmp$qvalues,
          rank = pathdata.kegg.gsea.tmp$rank,
          leading_edge = pathdata.kegg.gsea.tmp$leading_edge,
          core_enrichment = pathdata.kegg.gsea.tmp$core_enrichment
        )
                
      } else {
        geneList <- gene()
        pathdata.wiki.gsea.run <- gseWP(geneList,
                                        pvalueCutoff = input$adjPvalue,
                                        organism = "Homo sapiens",
                                        #exponent = 1,
                                        minGSSize = input$mingssize,
                                        maxGSSize = input$maxgssize,
                                        #eps = 0,
                                        pAdjustMethod = input$pvaladj,
                                        #verbose = FALSE,
                                        #by = "fgsea")
                                        )
        # GSEA figure
        pathdata.wiki.gsea.tmp <- as.data.frame(pathdata.wiki.gsea.run)
        
        tablewikiGSEA  = data.frame(
          #WikiID = pathdata.wiki.gsea.tmp$ID,
          WikiID = sapply(pathdata.wiki.gsea.tmp$ID, function(x)
            toString(tags$a(href=paste0("https://www.wikipathways.org//index.php?query=", x, "&title=Special%3ASearchPathways&doSearch=1&sa=Search"),target = "_blank", x))),
          Description = pathdata.wiki.gsea.tmp$Description,
          setSize = pathdata.wiki.gsea.tmp$setSize,
          enrichmentScore = pathdata.wiki.gsea.tmp$enrichmentScore,
          NES = pathdata.wiki.gsea.tmp$NES,
          pvalue = pathdata.wiki.gsea.tmp$pvalue,
          p.adjust = pathdata.wiki.gsea.tmp$p.adjust,
          #qval = pathdata.wiki.gsea.tmp$qvalues,
          rank = pathdata.wiki.gsea.tmp$rank,
          leading_edge = pathdata.wiki.gsea.tmp$leading_edge,
          core_enrichment = pathdata.wiki.gsea.tmp$core_enrichment
        )      
      }
    }
    #}
    
  })
  
  shared_df <- SharedData$new(pathwayDataTable)

  output$pathPlot <- renderPlotly({
    js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
    
    if(input$pathType == "1"){
    plot_ly() %>%
          add_trace(data = shared_df,
                    x = ~GeneRatio,
                    y = ~Description,
                    size = ~Count,
                    type = input$pathfigtype,
                    color = ~p.adjust,                    
                    colors = c(input$highcolopath, input$midcolopath, input$lowcolopath),
                    alpha = 1,
                    text = ~paste('</br> Pathway = ', Description,
                                  '</br> GeneRatio = ', GeneRatio,
                                  '</br> Count = ', Count,
                                  '</br> FDR = ', p.adjust),
                    name = "",
                    hoverinfo = 'text',
                    texttemplate = '%{y: .3s}', textposition = 'outside',
                    mode = "markers",
                    #customdata = urls_reactome_ora,
                    marker = list(sizeref = 0.3, color = ~factor(p.adjust))) %>%
          layout(title = input$pathmaintitle,
                 xaxis = list(title = input$pathxaxistitle), #showgrid = F
                 yaxis = list(title = input$pathyaxistitle),#showgrid = F
                 showlegend = TRUE,
                 legend=list(title = list(text='p.adjust'))
                 ) %>%
          config(displaylogo = FALSE) #%>% onRender(js)
    } else {
      plot_ly() %>%
          add_trace(data = shared_df,
                    x = ~enrichmentScore,
                    y = ~Description,
                    size = ~setSize,
                    type = input$pathfigtype,
                    color = ~p.adjust,                    
                    colors = c(input$highcolopath, input$midcolopath, input$lowcolopath),
                    alpha = 1,
                    text = ~paste('</br> Pathway = ', Description,
                                  '</br> enrichmentScore = ', enrichmentScore,
                                  '</br> Count = ', setSize,
                                  '</br> FDR = ', p.adjust),
                    name = "",
                    hoverinfo = 'text',
                    texttemplate = '%{y: .3s}', textposition = 'outside',
                    mode = "markers",
                    #customdata = urls_kegg_gsea,
                    marker = list(sizeref = 0.3, color = ~factor(p.adjust))) %>%
          layout(title = input$pathmaintitle,
                 xaxis = list(title = input$pathxaxistitle), #showgrid = F
                 yaxis = list(title = input$pathyaxistitle),#showgrid = F
                 showlegend = TRUE,
                 legend=list(title = list(text='p.adjust'))
                 ) %>%
          config(displaylogo = FALSE)
    }

})

output$pathtable <- renderDT({
  DT::datatable(
    shared_df,
    extensions = 'Buttons',
    escape = FALSE,
    filter = "bottom",
    style = "auto",
    selection = "single",
    rownames = FALSE,
    options = list(
      paging = TRUE,
      digits = 5,
      searching = TRUE,
      fixedColumns = TRUE,
      fixedColumns = list(leftColumns = 3),
      autowidth = TRUE,
      ordering = TRUE,
      dom = 'Bflrtip',
      scrollX = TRUE,
      lengthMenu = list(c(12, 20, 50, -1), c('12', '20', '50', 'All')),
      buttons = list(
        list(extend = "excel", text = "Download table", title= "", filename = "transcriptR_Pathway_Analysis_page",
             exportOptions = list(
               modifier = list(page = "current")
             )
        )
        # ,
        # list(extend = "excel", text = "Download full result", filename = "transcriptR_Pathway_Reactome_ORA_Alldata",
        #      exportOptions = list(
        #        modifier = list(page="all")
        #      ))
      )
    ), class = "display"
  )
}, server = FALSE)


  observe({
    shinyInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))
      }
      inputs
    }

    pathfile <- data.frame(pathwayDataTable())
    filepath_rows <- nrow(pathfile)

    pathfile$pathmap <- shinyInput(actionButton, filepath_rows, "button_", label = "Show map", onclick = 'Shiny.setInputValue(\"select_button\", this.id, {priority: \"event\"})')
  })

  observeEvent(input$runPM, {
    id_pathmap1 <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_pathmap1), add = TRUE)
    library(pathview)

    # Get the differential analysis table (gene name, logFC and p.value)
    pmdata <- input$file7
    ext <- tools::file_ext(pmdata$datapath)
    req(pmdata)

    filepathwaymap <- pmdata$datapath
    data.pm <- switch(ext,
               txt = readr::read_table(filepathwaymap),
               tsv = readr::read_tsv(filepathwaymap),
               csv = readr::read_csv(filepathwaymap),
               xls = readxl::read_xls(filepathwaymap),
               xlsx = readxl::read_xlsx(filepathwaymap))

    # data.pm <- read.delim(pmdata$datapath,
    #   sep = "\t",
    #   header = T
    # )

    data.pm1 <- data.frame(data.pm[[input$genesymbol_colm]], data.pm[[input$logfc_colm]])
    colnames(data.pm1) <- c("gene", "logFC")
    anaData <- data.pm1

    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    my.symbols <- anaData$gene
    # my.symbols <- anaData$external_gene_name

    gene2entrez <- bitr(
      my.symbols,
      OrgDb = hs,
      fromType = "SYMBOL",
      toType = "ENTREZID"
    )

    colnames(gene2entrez)[1] <- "gene"
    #colnames(anaData)[which(names(anaData) == "name")] <- "gene"
    mergeData <- merge(gene2entrez, anaData, by = "gene")
    test.fc <- mergeData$logFC

    names(test.fc) <- mergeData$ENTREZID

    pathwaydata <- input$file7
    ext <- tools::file_ext(pathwaydata$datapath)
    
    req(pathwaydata)
    #validate(need(ext == "csv", "Please upload a csv file"))
    
    filepathway <- pathwaydata$datapath
    pathdata <- switch(ext,
               txt = readr::read_table(filepathway),
               tsv = readr::read_tsv(filepathway),
               csv = readr::read_csv(filepathway),
               xls = readxl::read_xls(filepathway),
               xlsx = readxl::read_xlsx(filepathway))

    #pathdata <- read.delim(pathwaydata$datapath, sep = "\t")
    pathdata1 <- data.frame(pathdata[[input$genesymbol_colm]], pathdata[[input$logfc_colm]])
    colnames(pathdata1) <- c("genesymbol", "logFC")

            # purpose: 1) to get the entrezID, use mygene library. 2) mygene library needs taxa ID, 3) taxa ID matched with KEGG database and saved in www folder

    OrgDbReactome <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgReac == 'human'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$pathOrgReac == 'mouse') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgReac == 'rat') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgReac == 'celegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgReac == 'yeast') {
        namesdb <- print("org.Sc.sgd.db")}
      else if (input$pathOrgReac == 'zebrafish') {
        namesdb <- print("org.Dr.eg.db")}
      else {
        namesdb <- print("org.Dm.eg.db")}
      
      return(namesdb)
    #})
    }

    OrgDbWiki <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgWiki == 'Homo sapiens'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$ppathOrgWiki == 'Mus musculus') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgWiki == 'Rattus norvegicus') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgWiki == 'Arabidopsis thaliana') {
        namesdb <- print("org.At.tair.db")}
      else if (input$pathOrgWiki == 'Bos taurus') {
        namesdb <- print("org.Bt.eg.db")}
      else if (input$pathOrgWiki == 'Caenorhabditis elegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgWiki == 'Canis familiaris') {
        namesdb <- print("org.Cf.eg.db")}
      else if (input$pathOrgWiki == 'Danio rerio') {
        namesdb <- print("org.Dr.eg.db")}
      else if (input$pathOrgWiki == 'Drosophila melanogaster') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Gallus gallus') {
        namesdb <- print("org.Gg.eg.db")}
      else if (input$pathOrgWiki == 'Pan troglodytes') {
        namesdb <- print("org.Pt.eg.db")}
      else if (input$pathOrgWiki == 'Populus trichocarpa') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Saccharomyces cerevisiae') {
        namesdb <- print("org.Sc.sgd.db")}
      else {
        namesdb <- print("org.Ss.eg.db")}
      return(namesdb)
    #})
    }

    gene <- function() {
        if(input$pathDB == '1') {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbReactome())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")          

          return(geneList)

        } else if (input$pathDB == '2') {
          taxaData <- reactive({
            dbTable <- read.csv("www/kegg_organisms_taxonomyID.csv", header = F, sep = ",")
            matchID <- dbTable[dbTable$V2 == input$pathOrgKegg, ]
            return(matchID)
          }) #taxa ID save in taxaData()$V3

          library(mygene)
          geneIDconv <- queryMany(pathdata1$genesymbol, scopes="symbol", fields="entrezgene", 
                                  species=taxaData()$V3)
          
          new2 = data.frame(
              SYMBOL = geneIDconv$query,
              ENTREZID = geneIDconv$entrezgene
            )
          geneIDconvTable <- new2
          
          colnames(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconvTable, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          return(geneList)
        } else {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbWiki())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")
          
          return(geneList)
        }
    }  

    geneList <- gene()

    row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    mappath <- pathwayDataTable()$KEGGID[[row]]

    pathview(
      gene.data = geneList,
      species = input$pathOrgKegg,
      pathway.id = mappath,
      low = list(gene = input$lowcolopathmap),
      high = list(gene = input$highcolopathmap),
      mid = list(gene = input$midcolopathmap),
      kegg.native = T,
      discrete = list(gene = T),
      map.symbol = T,
      pdf.size = c(12, 10),
      same.layer = T,
      # limit=c(round(max(test.fc),2), round(min(test.fc), 2)),
      node.sum = "mean"
    )

  })


  pathmapSize <- reactive({
    return(input$pathmap_size)
  })


  map_pic <- reactiveVal(FALSE)
  observeEvent(input$runPM, map_pic(TRUE), ignoreInit = TRUE)
  # observeEvent(input$pathName, map_pic(TRUE), ignoreInit = TRUE)

  output$pathmap_result <- renderImage(
    {
        id_pathmap2 <- showNotification("still working, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_pathmap2), add = TRUE)
    #   # Get the differential analysis table (gene name, logFC and p.value)
      pmdata <- input$file7
      ext <- tools::file_ext(pmdata$datapath)
      req(pmdata)

      filepathwaymap <- pmdata$datapath
      data.pm <- switch(ext,
               txt = readr::read_table(filepathwaymap),
               tsv = readr::read_tsv(filepathwaymap),
               csv = readr::read_csv(filepathwaymap),
               xls = readxl::read_xls(filepathwaymap),
               xlsx = readxl::read_xlsx(filepathwaymap))

      # data.pm <- read.delim(pmdata$datapath,
      #   sep = "\t",
      #   header = T
      # )

    data.pm1 <- data.frame(data.pm[[input$genesymbol_colm]], data.pm[[input$logfc_colm]])
    colnames(data.pm1) <- c("gene", "logFC")
    anaData <- data.pm1

    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    my.symbols <- anaData$gene

      gene2entrez <- bitr(
        my.symbols,
        OrgDb = hs,
        fromType = "SYMBOL",
        toType = "ENTREZID"
      )

      colnames(gene2entrez)[1] <- "gene"
      #colnames(anaData)[which(names(anaData) == "name")] <- "gene"
      mergeData <- merge(gene2entrez, anaData, by = "gene")
      test.fc <- mergeData$logFC

      names(test.fc) <- mergeData$ENTREZID

      pathwaydata <- input$file7
    ext <- tools::file_ext(pathwaydata$datapath)
    
    req(pathwaydata)
    #validate(need(ext == "csv", "Please upload a csv file"))
    
    filepathway <- pathwaydata$datapath
    pathdata <- switch(ext,
               txt = readr::read_table(filepathway),
               tsv = readr::read_tsv(filepathway),
               csv = readr::read_csv(filepathway),
               xls = readxl::read_xls(filepathway),
               xlsx = readxl::read_xlsx(filepathway))

    #pathdata <- read.delim(pathwaydata$datapath, sep = "\t")
    pathdata1 <- data.frame(pathdata[[input$genesymbol_colm]], pathdata[[input$logfc_colm]])
    colnames(pathdata1) <- c("genesymbol", "logFC")

            # purpose: 1) to get the entrezID, use mygene library. 2) mygene library needs taxa ID, 3) taxa ID matched with KEGG database and saved in www folder

    OrgDbReactome <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgReac == 'human'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$pathOrgReac == 'mouse') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgReac == 'rat') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgReac == 'celegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgReac == 'yeast') {
        namesdb <- print("org.Sc.sgd.db")}
      else if (input$pathOrgReac == 'zebrafish') {
        namesdb <- print("org.Dr.eg.db")}
      else {
        namesdb <- print("org.Dm.eg.db")}
      
      return(namesdb)
    #})
    }

    OrgDbWiki <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgWiki == 'Homo sapiens'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$ppathOrgWiki == 'Mus musculus') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgWiki == 'Rattus norvegicus') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgWiki == 'Arabidopsis thaliana') {
        namesdb <- print("org.At.tair.db")}
      else if (input$pathOrgWiki == 'Bos taurus') {
        namesdb <- print("org.Bt.eg.db")}
      else if (input$pathOrgWiki == 'Caenorhabditis elegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgWiki == 'Canis familiaris') {
        namesdb <- print("org.Cf.eg.db")}
      else if (input$pathOrgWiki == 'Danio rerio') {
        namesdb <- print("org.Dr.eg.db")}
      else if (input$pathOrgWiki == 'Drosophila melanogaster') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Gallus gallus') {
        namesdb <- print("org.Gg.eg.db")}
      else if (input$pathOrgWiki == 'Pan troglodytes') {
        namesdb <- print("org.Pt.eg.db")}
      else if (input$pathOrgWiki == 'Populus trichocarpa') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Saccharomyces cerevisiae') {
        namesdb <- print("org.Sc.sgd.db")}
      else {
        namesdb <- print("org.Ss.eg.db")}
      return(namesdb)
    #})
    }

    gene <- function() {
        if(input$pathDB == '1') {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbReactome())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")          

          return(geneList)

        } else if (input$pathDB == '2') {
          taxaData <- reactive({
            dbTable <- read.csv("www/kegg_organisms_taxonomyID.csv", header = F, sep = ",")
            matchID <- dbTable[dbTable$V2 == input$pathOrgKegg, ]
            return(matchID)
          }) #taxa ID save in taxaData()$V3

          library(mygene)
          geneIDconv <- queryMany(pathdata1$genesymbol, scopes="symbol", fields="entrezgene", 
                                  species=taxaData()$V3)
          
          new2 = data.frame(
              SYMBOL = geneIDconv$query,
              ENTREZID = geneIDconv$entrezgene
            )
          geneIDconvTable <- new2
          
          colnames(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconvTable, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          return(geneList)
        } else {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbWiki())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")
          
          return(geneList)
        }
    }  

      geneList <- gene()

      row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
      mappath <- pathwayDataTable()$KEGGID[[row]]

      pathview(
        gene.data = geneList,
        species = input$pathOrgKegg,
        #pathway.id = input$pathName,
        pathway.id = mappath,
        low = list(gene = input$lowcolopathmap),
        high = list(gene = input$highcolopathmap),
        mid = list(gene = input$midcolopathmap),
        kegg.native = T,
        discrete = list(gene = T),
        map.symbol = T,
        pdf.size = c(12, 10),
        same.layer = T,
        # limit=c(round(max(test.fc),2), round(min(test.fc), 2)),
        node.sum = "mean"
      )      
      

      req(map_pic())
      list(
        #src = normalizePath(file.path(paste0(input$pathName, ".pathview.png"))),
        src = normalizePath(file.path(paste0(mappath, ".pathview.png"))),
        contentType = "image/png" #,
        # width = pathmapSize,
        # height = pathmapSize,
        # outputArgs = list()
      )
    },
    deleteFile = FALSE
  )

  ## Gene information
  #=================================

  observeEvent(input$showpathmapgene, {

  showModal(modalDialog(
  title = "Gene Information",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
  output$geneinfor <- DT::renderDataTable(server = FALSE,{
    id_pathmapgene <- showNotification("generating information table, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_pathmapgene), add = TRUE)
      # Get the differential analysis table (gene name, logFC and p.value)
      # get user input data file for logFC
      pmdata <- input$file7
      ext <- tools::file_ext(pmdata$datapath)
      req(pmdata)

      filepathwaymap <- pmdata$datapath
      data.pm <- switch(ext,
               txt = readr::read_table(filepathwaymap),
               tsv = readr::read_tsv(filepathwaymap),
               csv = readr::read_csv(filepathwaymap),
               xls = readxl::read_xls(filepathwaymap),
               xlsx = readxl::read_xlsx(filepathwaymap))

      # data.pm <- read.delim(pmdata$datapath,
      #   sep = "\t",
      #   header = T
      # )

    data.pm1 <- data.frame(data.pm[[input$genesymbol_colm]], data.pm[[input$logfc_colm]])
    colnames(data.pm1) <- c("gene", "logFC")
    anaData <- data.pm1

    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    my.symbols <- anaData$gene

      gene2entrez <- bitr(
        my.symbols,
        OrgDb = hs,
        fromType = "SYMBOL",
        toType = "ENTREZID"
      )

      colnames(gene2entrez)[1] <- "gene"
      #colnames(anaData)[which(names(anaData) == "name")] <- "gene"
      mergeData <- merge(gene2entrez, anaData, by = "gene")
      test.fc <- mergeData$logFC

      names(test.fc) <- mergeData$ENTREZID
      pathwaydata <- input$file7
    ext <- tools::file_ext(pathwaydata$datapath)
    
    req(pathwaydata)
    #validate(need(ext == "csv", "Please upload a csv file"))
    
    filepathway <- pathwaydata$datapath
    pathdata <- switch(ext,
               txt = readr::read_table(filepathway),
               tsv = readr::read_tsv(filepathway),
               csv = readr::read_csv(filepathway),
               xls = readxl::read_xls(filepathway),
               xlsx = readxl::read_xlsx(filepathway))

    #pathdata <- read.delim(pathwaydata$datapath, sep = "\t")
    pathdata1 <- data.frame(pathdata[[input$genesymbol_colm]], pathdata[[input$logfc_colm]])
    colnames(pathdata1) <- c("genesymbol", "logFC")

            # purpose: 1) to get the entrezID, use mygene library. 2) mygene library needs taxa ID, 3) taxa ID matched with KEGG database and saved in www folder

    OrgDbReactome <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgReac == 'human'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$pathOrgReac == 'mouse') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgReac == 'rat') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgReac == 'celegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgReac == 'yeast') {
        namesdb <- print("org.Sc.sgd.db")}
      else if (input$pathOrgReac == 'zebrafish') {
        namesdb <- print("org.Dr.eg.db")}
      else {
        namesdb <- print("org.Dm.eg.db")}
      
      return(namesdb)
    #})
    }

    OrgDbWiki <- function(){
      #observeEvent(input$get_sp, {
      if(input$pathOrgWiki == 'Homo sapiens'){
        namesdb <- print("org.Hs.eg.db")} 
      else if (input$ppathOrgWiki == 'Mus musculus') {
        namesdb <- print("org.Mm.eg.db")}
      else if (input$pathOrgWiki == 'Rattus norvegicus') {
        namesdb <- print("org.Rn.eg.db")}
      else if (input$pathOrgWiki == 'Arabidopsis thaliana') {
        namesdb <- print("org.At.tair.db")}
      else if (input$pathOrgWiki == 'Bos taurus') {
        namesdb <- print("org.Bt.eg.db")}
      else if (input$pathOrgWiki == 'Caenorhabditis elegans') {
        namesdb <- print("org.Ce.eg.db")}
      else if (input$pathOrgWiki == 'Canis familiaris') {
        namesdb <- print("org.Cf.eg.db")}
      else if (input$pathOrgWiki == 'Danio rerio') {
        namesdb <- print("org.Dr.eg.db")}
      else if (input$pathOrgWiki == 'Drosophila melanogaster') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Gallus gallus') {
        namesdb <- print("org.Gg.eg.db")}
      else if (input$pathOrgWiki == 'Pan troglodytes') {
        namesdb <- print("org.Pt.eg.db")}
      else if (input$pathOrgWiki == 'Populus trichocarpa') {
        namesdb <- print("org.Dm.eg.db")}
      else if (input$pathOrgWiki == 'Saccharomyces cerevisiae') {
        namesdb <- print("org.Sc.sgd.db")}
      else {
        namesdb <- print("org.Ss.eg.db")}
      return(namesdb)
    #})
    }

    gene <- function() {
        if(input$pathDB == '1') {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbReactome())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")          

          return(geneList)

        } else if (input$pathDB == '2') {
          taxaData <- reactive({
            dbTable <- read.csv("www/kegg_organisms_taxonomyID.csv", header = F, sep = ",")
            matchID <- dbTable[dbTable$V2 == input$pathOrgKegg, ]
            return(matchID)
          }) #taxa ID save in taxaData()$V3

          library(mygene)
          geneIDconv <- queryMany(pathdata1$genesymbol, scopes="symbol", fields="entrezgene", 
                                  species=taxaData()$V3)
          
          new2 = data.frame(
              SYMBOL = geneIDconv$query,
              ENTREZID = geneIDconv$entrezgene
            )
          geneIDconvTable <- new2
          
          colnames(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconvTable, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          return(geneList)
        } else {
          geneIDconv <- bitr(pathdata1$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDbWiki())
          names(pathdata1)[1] <- "SYMBOL"
          pathdata2 <- merge(pathdata1, geneIDconv, by = "SYMBOL")
          geneList <- pathdata2[,2]
          names(geneList) <- as.character(pathdata2[,3])
          geneList = sort(geneList, decreasing = TRUE)
          geneData <- names(geneList)[abs(geneList) > 0]

          message("finished Reactome GeneID Conversion")
          
          return(geneList)
        }
    }  

      geneList <- gene()

      row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
      mappath <- pathwayDataTable()$KEGGID[[row]]

      pv.data <- pathview(
        gene.data = geneList,
        species = input$pathOrgKegg,
        #pathway.id = input$pathName,
        pathway.id = mappath,
        low = list(gene = input$lowcolopathmap),
        high = list(gene = input$highcolopathmap),
        mid = list(gene = input$midcolopathmap),
        kegg.native = T,
        discrete = list(gene = T),
        map.symbol = T,
        pdf.size = c(12, 10),
        same.layer = T,
        # limit=c(round(max(test.fc),2), round(min(test.fc), 2)),
        node.sum = "mean"
      )   
      

    pv.data.res <- pv.data$plot.data.gene[c(1,2,9)]

    colnames(pv.data.res) <- c("KEGG_name", "Gene_symbol", "log2FC")

    pathgenelist3 <- data.frame(
      Gene_symbol = sapply(pv.data.res$Gene_symbol, function(x) {
        toString(tags$a(href = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", x), target = "_blank", x))
      }),
      log2FC = pv.data.res$log2FC
    )

    #  return(entrez2gene)
    DT::datatable(
      pathgenelist3,
      extensions = c("FixedColumns", "Buttons"),
      escape = FALSE,
      rownames = FALSE,
      filter = "bottom",
      style = "auto",
      selection = "single",
      options = list(
        paging = TRUE,
        searching = TRUE,
        autowidth = TRUE,
        ordering = TRUE,
        dom = "Bflrtip",
        scrollX = TRUE,
        lengthMenu = list(c(10, 50, 100, -1), c("10", "50", "100", "All")),
        buttons = list(
          list(
            extend = "excel", text = "Download Table", title = "", filename = "Pathway Gene List",
            exportOptions = list(
              modifier = list(page = "current")
            )
          )
        )
        ), class = "display"
      ) %>%
      formatStyle("log2FC",
        color = styleInterval(c(0), c("blue", "red"))
      ) %>%
      formatStyle(
        columns = "Gene_symbol",
        valueColumns = "log2FC",
        color = styleInterval(c(0), c("blue", "red"))
      )
  }),
  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitpmgeneinfo", "OK")
              )
            )
          )

  })

  observeEvent(input$submitpmgeneinfo, {
  removeModal()
})

## download pathmaps
  observeEvent(input$copyfilepm, {
    dir.create(input$dirname)

    filelist <- list.files(getwd(), "[.]pathview[.]png$", full.names = TRUE)
    file.copy(filelist, input$dirname, overwrite = TRUE)

    unlink(filelist, recursive = FALSE, force = FALSE)
  } 
  )

observeEvent(input$download_pmdata, {
        showModal(modalDialog(
          title = HTML('<span style="color:SteelBlue; font-size: 26px; font-weight:bold; font-family:sans-serif ">Proceed to download results<span>
                    <button type = "button" class="close" data-dismiss="modal" ">
                    </button> '),
          HTML('<span style="color:LightCoral; font-size: 20px; font-weight:bold; font-family:sans-serif ">NOTE: It is a good idea to close the browser once you finished one analysis and run the next one. The cache memory will be cleared in this way ...<span>
                    <button type = "button" class="close" data-dismiss="modal" ">
                    </button> '),
          footer = actionButton("ConfirmChAMP", "Confirm"),

          easyClose = FALSE
        ))
        })

        observeEvent(input$ConfirmChAMP, {
          showModal(
              modalDialog(
                  HTML('<span style="color:LightSeaGreen; font-size: 26px; font-weight:bold; font-family:sans-serif ">Are you sure?<span>
                    <button type = "button" class="close" data-dismiss="modal" ">
                    </button> '),
                  footer = tagList(
                      downloadButton(outputId = "download_data", "Yes"),
                      modalButton("No")
                  ),
                  easyClose = FALSE
              )
          )
        })

      
  output$download_data <- downloadHandler(  
        filename = function(){
          paste("TranscriptR_results", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){   
        showNotification("zipping results, please wait...", 
                          duration = 5,
                          type = "message", 
                          closeButton = TRUE)
        removeModal()
          
          results <- input$dirname

          filesToSave <- list.files(results, full.names= TRUE)
      
          system2("zip", args=(paste(file,filesToSave,sep=" ")))

        },
          contentType = "application/zip"    
      )


}

shinyCatch(pathwayAnalysis(), position = "bottom-right", blocking_level = "none", shiny = TRUE, prefix = "transcriptR", trace_back = spsOption("traceback"))

####===============Pathway analysis module end =====#######


#####============== Volcano plot module Start==============#######

          knitr_table <- function(x) {
            x %>% 
              knitr::kable(format = "html", digits = Inf, 
                          format.args = list(big.mark = ",")) %>%
              kableExtra::kable_styling(font_size = 15)
          }
          
          volpl = eventReactive(input$submitvol, {
            volplotpath <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
            on.exit(removeNotification(volplotpath), add = TRUE)

          dat.volc <- file1volc()

          data <- data.frame(dat.volc[[input$geneName_volc]], dat.volc[[input$logfc_volc]], dat.volc[[input$pvalue_volc]])

          colnames(data)[1] <- "ID"
          colnames(data)[2] <- "logFC"
          colnames(data)[3] <- "adjPValue"
          
          data <- data %>% 
            mutate(
              Expression = case_when(logFC >= (input$fc) & adjPValue <= 0.05 ~ "upregulated",
                                    logFC <= -(input$fc) & adjPValue <= 0.05 ~ "downregulated",
                                    TRUE ~ "nonsignificant")
            )        
          
          data <- data %>% 
            mutate(
              Significance = case_when(
                abs(logFC) >= (input$fc) & adjPValue <= 0.05 & adjPValue > 0.001 ~ "adjPValue 0.001", 
                abs(logFC) >= (input$fc) & adjPValue <= 0.001 & adjPValue > 0.0001 ~ "adjPValue 0.0001",
                abs(logFC) >= (input$fc) & adjPValue <= 0.0001 ~ "adjPValue 0.00001", 
                TRUE ~ "Unchanged")
            )

          top <- input$niVol
          top_genes <- bind_rows(
            data %>% 
              filter(Expression == 'upregulated') %>% 
              arrange(adjPValue, desc(abs(logFC))) %>% 
              head(top),
            data %>% 
              filter(Expression == 'downregulated') %>% 
              arrange(adjPValue, desc(abs(logFC))) %>% 
              head(top)
          )

          p2 <- ggplot(data, aes(logFC, -log(adjPValue,10))) +
            geom_point(aes(color = Expression), size = input$pointsize2) +
            #xlab(expression("log"[2]),input$xaxislab_volc) + 
            #ylab(expression("-log"[10]),input$yaxislab_volc) +
            xlab(bquote("log"[2] ~ .(input$xaxislab_volc))) + 
            ylab(bquote("-log"[10] ~ .(input$yaxislab_volc))) +
            ggtitle(input$title_volc) +
            scale_color_manual(values = c(input$highcol, input$nonsigcol, input$lowcol)) +
            guides(colour = guide_legend(override.aes = list(size=2)))
          
          p3 <- p2 +
            xlim(input$xlim_volc) + ylim(input$ylim_volc) +
            geom_label_repel(data = top_genes, label.size = NA,
                            mapping = aes(logFC, -log(adjPValue,10), label = ID),
                            size = input$wordsize2,
                            max.overlaps = 3000) +
            theme_minimal() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grids
            theme(axis.line = element_line(colour = "black")) + # add axes lines
            theme(legend.title = element_blank(), text = element_text(size = input$wordsize)) +
            geom_vline(xintercept=c(-input$fc, input$fc), col=input$vlinecol, 
                      linetype = input$linetype) + # add vline/logFC cut-off line
            geom_hline(yintercept=-log10(0.05), col=input$hlinecol, 
                      linetype = input$linetype) + # add adjPVal cut-off line
            theme(text = element_text(size = input$wordsize,
                                      family = input$fontfamily,
                                      color = input$volfontcol))
          
          p3
          }
          )
          
          output$volcanoplot <- renderPlot({
            volpl()
          })

file1volc = reactive({
    filevolc1 <- input$file1_vol
    ext <- tools::file_ext(filevolc1$datapath)
    req(filevolc1)
    
    filepath <- filevolc1$datapath
    filevolc1 <- switch(ext,
               txt = readr::read_table(filepath),
               tsv = readr::read_tsv(filepath),
               csv = readr::read_csv(filepath),
               xls = readxl::read_xls(filepath),
               xlsx = readxl::read_xlsx(filepath))
    
    updateSelectInput(session, "geneName_volc", choices = names(filevolc1))
    updateSelectInput(session, "logfc_volc", choices = names(filevolc1))
    updateSelectInput(session, "pvalue_volc", choices = names(filevolc1))

    return(filevolc1)
  
  }
  )

  output$volcanofile1 <- DT::renderDataTable( 
    file1volc(), extensions = "Buttons",
          escape = FALSE,
          rownames = FALSE,
          filter = "bottom",
          style = "auto",
          selection = "single",
          options = list(
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autowidth = TRUE,
            ordering = TRUE,
            dom = "Bflrtip",
            scrollX = TRUE,
            lengthMenu = list(c(6, 20, 100, -1), c("6", "20", "100", "All"))
    )
  )
          output$volcDownload <- downloadHandler(
          filename = function() {
            paste("VolcanoPlot", Sys.Date(), tolower(input$filetype_volc), sep = ".")},
          content = function(file) {
          if(input$filetype_volc == "PNG") 
          {png(file, 
            width = input$volc_pngwidth , 
            height = input$volc_pngheight, 
            units = input$volc_pngunits, 
            res = input$volc_pngresol, 
            type = "cairo")} 
          else if(input$filetype_volc == "PDF")
          {pdf(file,
            width = input$volc_pdfwidth , 
            height = input$volc_pdfheight            
            )}
          else if(input$filetype_volc == "SVG")
          {svg(file, 
            width = input$volc_svgwidth , 
            height = input$volc_svgheight
            )}
          else (tiff(file, 
            width = input$volc_tiffwidth , 
            height = input$volc_tiffheight, 
            units = input$volc_tiffunits, 
            res = input$volc_tiffresol,
          type = "cairo"))
            

          print(volpl())

          dev.off()
          }
          )

# interactive Volcano plot
volplot1 <- eventReactive(input$runVol, {
    library(dplyr)
    library(stringr) 
    id_intervol1 <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_intervol1), add = TRUE)

    voldata1 <- as.data.frame(file1volc())
    colnames(voldata1) <- c("ID", "logFC", "Pvalue")
    voldata1["group"] <- "NotSignificant"
    voldata1[which(voldata1['Pvalue'] < input$volp & abs(voldata1['logFC']) < input$volfc), "group"] <- "Significant"
    voldata1[which(voldata1['Pvalue'] > input$volp & abs(voldata1['logFC']) > input$volfc), "group"] <- "FoldChange"
    voldata1[which(voldata1['Pvalue'] < input$volp & abs(voldata1['logFC']) > input$volfc), "group"] <- "Significant&FoldChange"
    sfc <- grep("Significant&FoldChange", voldata1$group)
    top_peaks <- voldata1[with(voldata1, sfc),]

    a <- list()
    for(i in seq_len(nrow(top_peaks))){
      m <- top_peaks[i, ]
      a[[i]] <- list(
        x = m[["logFC"]],
        y = -log10(m[["Pvalue"]]),
        text = m[["ID"]],
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40
      )
    }

    urls_vol <- str_c("https://www.ncbi.nlm.nih.gov/gene/?term=", voldata1$ID)

    volpl <- plot_ly(data = voldata1,
                     x = voldata1$logFC,
                     y = -log10(voldata1$Pvalue),
                     text = voldata1$ID,
                     mode = "markers",
                     color = voldata1$group) %>%
      layout(title = "Volcano plot") %>%
      layout(annotation = a) %>% config(displaylogo = FALSE)

    js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"  

    volpl %>% onRender(js)
  })

output$volplot <- renderPlotly({    
    volplot1()
    
  })

  output$downloadvolcano <- downloadHandler(
    filename = function() { paste0("Volcano_Plot", ".html", sep='') },
    content = function(figure2) {
      htmlwidgets::saveWidget(as_widget(volplot1()), figure2, selfcontained = TRUE)
    })

#####============== Volcano plot module end==============#######


#####============== Gene Ontology module Start==============#######
##========================================####
  #6. Gene Ontology analysis
##========================================####

  #6.1 GO enrichment plot
  goPlot1 <- eventReactive(input$goRun, {
    library(dplyr)
    library(stringr)
    id_go1 <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_go1), add = TRUE)
    
    file6 <- input$file6
    ext <- tools::file_ext(file6$datapath)
    
    req(file6)
    
    filegopath <- file6$datapath
    godata <- switch(ext,
               txt = readr::read_table(filegopath),
               tsv = readr::read_tsv(filegopath),
               csv = readr::read_csv(filegopath),
               xls = readxl::read_xls(filegopath),
               xlsx = readxl::read_xlsx(filegopath))

    #godata <- read.delim(file6$datapath, sep = "\t", row.names = 1)
    
    # prepare the input gene list for pathway analysis
    godata1 <- cbind(data.frame(godata$gene), data.frame(godata$logFC))
    geneIDconv <- bitr(godata1$godata.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    names(godata1)[1] <- "SYMBOL"
    godata2 <- merge(godata1, geneIDconv, by = "SYMBOL")
    geneList <- godata2[,2]
    names(geneList) <- as.character(godata2[,3])
    geneList = sort(geneList, decreasing = TRUE)
    gene <- names(geneList)[abs(geneList) > 0]
    
    datatype <- reactive({
      switch(input$goType,
             "Over representation analysis" = 1,
             "Gene set enrichment analysis" = 2)
    })
    
    if(input$goType == "1") {
      godata.ora.run <- enrichGO(gene = gene,
                                 #universe= names(geneList),
                                 OrgDb = org.Hs.eg.db,
                                 ont= "ALL",
                                 pAdjustMethod = input$pvaladjgo,
                                 pvalueCutoff  = input$numpval,
                                 qvalueCutoff  = input$numqval,
                                 #eps = 0,
                                 readable = TRUE)
      
      godata.ora.tmp <- data.frame(godata.ora.run@result)
      
      tmp1 <- godata.ora.tmp %>%
        separate(GeneRatio, c("GeneR1", "GeneR2"), "/")
      tmp1$GeneRatio <- as.numeric(tmp1$GeneR1)/as.numeric(tmp1$GeneR2)
      godata.ora.tmp$GeneRatio <- tmp1$GeneRatio
      tmp2 <- godata.ora.tmp %>%
        separate(BgRatio, c("bg1", "bg2"), "/")
      tmp2$BgRatio <- as.numeric(tmp2$bg1)/as.numeric(tmp2$bg2)
      godata.ora1 <- godata.ora.tmp
      
      data <- data.frame(godata.ora1$ID, godata.ora1$ONTOLOGY, godata.ora1$p.adjust, godata.ora1$Description, stringsAsFactors = FALSE)
      data$godata.ora1.ONTOLOGY <- sort(data$godata.ora1.ONTOLOGY, decreasing = TRUE)
      
      if(input$ontology == "BP") {
        data.bp <- data[data$godata.ora1.ONTOLOGY == "BP",]
        data.bp.num <- data.bp[1:input$gonum,]
        urls_bp <- str_c("http://amigo.geneontology.org/amigo/term/",data.bp.num$godata.ora1.ID)
        
        figure_bp <- plot_ly(data.bp.num, 
                             x= ~godata.ora1.p.adjust, 
                             y = ~godata.ora1.Description, 
                             type = "bar", 
                             customdata = urls_bp,
                             marker = list(color = 'rgb(0,73,83')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Biological Process'))
        
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_bp %>% onRender(js)
        
      } else if (input$ontology == "CC"){
        data.cc <- data[data$godata.ora1.ONTOLOGY == "CC",]
        data.cc.num <- data.cc[1:input$gonum,]
        
        urls_cc <- str_c("http://amigo.geneontology.org/amigo/term/",data.cc.num$godata.ora1.ID)
        
        figure_cc <- plot_ly(data.cc.num, x= ~godata.ora1.p.adjust, y = ~godata.ora1.Description, 
                             type = "bar",customdata = urls_cc, marker = list(color = 'rgb(36,33,36')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Cellular Component'))
        
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_cc %>% onRender(js)
        
      } else {
        data.mf <- data[data$godata.ora1.ONTOLOGY == "MF",]
        data.mf.num <- data.mf[1:input$gonum,]
        
        urls_mf <- str_c("http://amigo.geneontology.org/amigo/term/",data.mf.num$godata.ora1.ID)
        
        figure_mf <- plot_ly(data.mf.num, x= ~godata.ora1.p.adjust, y = ~godata.ora1.Description, 
                             type = "bar",customdata = urls_mf, marker = list(color = 'rgb(54,69,79')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Molecular Function'))
        
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_mf %>% onRender(js)
        
      }
    } else {
      godata.gsea.run <- gseGO(geneList = geneList,
                               OrgDb = org.Hs.eg.db,
                               ont = "ALL",
                               minGSSize = 10,
                               maxGSSize = 500,
                               eps = 0,
                               pvalueCutoff = input$numpval,
                               verbose = FALSE)
      # GSEA figure
      godata.gsea.tmp <- as.data.frame(godata.gsea.run@result)
      # visualization of the pathway analysis result
      godata.gsea1 <- godata.gsea.tmp#[1:input$num,]
      
      data.gsea <- data.frame(godata.gsea1$ID, godata.gsea1$ONTOLOGY, godata.gsea1$p.adjust, godata.gsea1$Description)
      data.gsea$godata.gsea1.result.ONTOLOGY <- sort(data.gsea$godata.gsea1.result.ONTOLOGY, decreasing = TRUE)
      
      if(input$ontology == "BP") {
        data.gsea.bp <- data.gsea[data.gsea$godata.gsea1.ONTOLOGY == "BP",]
        data.gsea.bp.num <- data.gsea.bp[1:input$gonum,]
        
        urls_gsea_bp <- str_c("http://amigo.geneontology.org/amigo/term/",data.gsea.bp.num$godata.gsea1.ID)
        
        figure_bp_gsea <- plot_ly(data.gsea.bp.num, x= ~godata.gsea1.p.adjust, y = ~godata.gsea1.Description, 
                                  type = "bar",customdata = urls_gsea_bp, marker = list(color = 'rgb(0,73,83')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Biological Process'))
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_bp_gsea %>% onRender(js)
      } else if (input$ontology == "CC"){
        data.gsea.cc <- data.gsea[data.gsea$godata.gsea1.ONTOLOGY == "CC",]
        data.gsea.cc.num <- data.gsea.cc[1:input$gonum,]
        urls_gsea_cc <- str_c("http://amigo.geneontology.org/amigo/term/",data.gsea.cc.num$godata.gsea1.ID)
        figure_cc_gsea <- plot_ly(data.gsea.cc.num, x= ~godata.gsea1.p.adjust, y = ~godata.gsea1.Description, 
                                  type = "bar", customdata = urls_gsea_cc, marker = list(color = 'rgb(36,33,36')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Cellular Component'))
        
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_cc_gsea %>% onRender(js)
      } else {
        data.gsea.mf <- data.gsea[data.gsea$godata.gsea1.ONTOLOGY == "MF",]
        data.gsea.mf.num <- data.gsea.mf[1:input$gonum,]
        urls_gsea_mf <- str_c("http://amigo.geneontology.org/amigo/term/",data.gsea.mf.num$godata.gsea1.ID)
        
        figure_mf_gsea <- plot_ly(data.gsea.mf.num, x= ~godata.gsea1.p.adjust, y = ~godata.gsea1.Description, 
                                  type = "bar", customdata = urls_gsea_mf, marker = list(color = 'rgb(54,69,79')) %>%
          layout(xaxis=list(title = 'adjusted p-value',type='category', 
                            nticks = 10),
                 yaxis = list(title = 'Gene Ontology: Molecular Function'))
        js <- "
            function(el, x) {
                el.on('plotly_click', function(d) {
                    var point = d.points[0];
                    var url = point.data.customdata[point.pointIndex];
                window.open(url);
              });
            }"
        
        figure_mf_gsea %>% onRender(js)
        
      }
    }
  })



  output$goplot <- renderPlotly({
    # if(input$goRun) {goPlot1()}
    # else {goPlot2()}
    goPlot1()
  })

  output$downloadGO <- downloadHandler(
    filename = function() { paste0("Figure_GO_Enrichment_Analysis", ".html", sep='') },
    content = function(figure) {
      htmlwidgets::saveWidget(as_widget(goPlot1()), figure, selfcontained = TRUE)
    })

  output$tableGO <- DT::renderDataTable(server=FALSE,{
    id_gotab2 <- showNotification("generating table, please wait...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id_gotab2), add = TRUE)
    
    file6 <- input$file6
    ext <- tools::file_ext(file6$datapath)
    
    req(file6)
    
    filegopath <- file6$datapath
    godata <- switch(ext,
               txt = readr::read_table(filegopath),
               tsv = readr::read_tsv(filegopath),
               csv = readr::read_csv(filegopath),
               xls = readxl::read_xls(filegopath),
               xlsx = readxl::read_xlsx(filegopath))

    #godata <- read.delim(file6$datapath, sep = "\t", row.names = 1)
    
    # prepare the input gene list for pathway analysis
    godata1 <- cbind(data.frame(godata$gene), data.frame(godata$logFC))
    geneIDconv <- bitr(godata1$godata.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    names(godata1)[1] <- "SYMBOL"
    godata2 <- merge(godata1, geneIDconv, by = "SYMBOL")
    geneList <- godata2[,2]
    names(geneList) <- as.character(godata2[,3])
    geneList = sort(geneList, decreasing = TRUE)
    gene <- names(geneList)[abs(geneList) > 0]
    
    datatype <- reactive({
      switch(input$goType,
             "Over representation analysis" = 1,
             "Gene set enrichment analysis" = 2)
    })
    
    if(input$goType == "1") {
      
      godata.ora.run <- enrichGO(gene = gene,
                                 #universe= names(geneList),
                                 OrgDb = org.Hs.eg.db,
                                 ont= "ALL",
                                 pAdjustMethod = input$pvaladjgo,
                                 pvalueCutoff  = input$numpval,
                                 qvalueCutoff  = input$numqval,
                                 #eps = 0,
                                 readable = TRUE)
      
      godata.ora.tmp <- data.frame(godata.ora.run@result)
      #godata.reactome.ora1 <- as.data.frame(godata.reactome.ora.tmp)
      tmp1 <- godata.ora.tmp %>%
        separate(GeneRatio, c("GeneR1", "GeneR2"), "/")
      tmp1$GeneRatio <- as.numeric(tmp1$GeneR1)/as.numeric(tmp1$GeneR2)
      godata.ora.tmp$GeneRatio <- tmp1$GeneRatio
      tmp2 <- godata.ora.tmp %>%
        separate(BgRatio, c("bg1", "bg2"), "/")
      tmp2$BgRatio <- as.numeric(tmp2$bg1)/as.numeric(tmp2$bg2)
      godata.ora.tmp$BgRatio <- tmp2$BgRatio
      
      tableGOORA  = data.frame(
        #GOID = godata.ora.tmp$ID,
        GOID = sapply(godata.ora.tmp$ID, function(x)
          toString(tags$a(href=paste0("http://amigo.geneontology.org/amigo/term/", x),target = "_blank", x))),
        GOclassification = godata.ora.tmp$ONTOLOGY,
        GOdescription = godata.ora.tmp$Description,
        GeneRatio = godata.ora.tmp$GeneRatio,
        BgRatio = godata.ora.tmp$BgRatio,
        pvalue = godata.ora.tmp$pvalue,
        p.adjust = godata.ora.tmp$p.adjust,
        qval = godata.ora.tmp$qvalue,
        geneID = godata.ora.tmp$geneID,
        GeneCount = godata.ora.tmp$Count)
      
      DT::datatable(
        tableGOORA,
        escape = FALSE,
        extensions = 'Buttons',
        rownames = FALSE,
        filter = "bottom",
        style = "auto",
        selection = "single",
        options = list(
          paging = TRUE,
          searching = TRUE,
          fixedColumns = TRUE,
          autowidth = TRUE,
          ordering = TRUE,
          dom = 'Bflrtip',
          scrollX = TRUE,
          lengthMenu = list(c(20, 50,100, -1), c('20', '50', '100', 'All')),
          buttons = list(
            list(extend = "excel", text = "Download Table", title = "", filename = "transcriptR_GO_ORA_page",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            )
          )
        ), class = "display"
      )} else {
        godata.gsea.run <- gseGO(geneList = geneList,
                                 OrgDb = org.Hs.eg.db,
                                 ont = "ALL",
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 pvalueCutoff = input$numpval,
                                 keyType = "ENTREZID",
                                 exponent = 1,
                                 eps = 0,
                                 seed = FALSE,
                                 by = "fgsea",
                                 verbose = FALSE)
        # GSEA figure
        godata.gsea.tmp <- data.frame(godata.gsea.run@result)
        
        godata.gsea.tmp$ID = sapply(godata.gsea.tmp$ID, function(x)
          toString(tags$a(href=paste0("http://amigo.geneontology.org/amigo/term/", x),target = "_blank", x)))
        
        DT::datatable(
          #tableGSEA.go,
          godata.gsea.tmp,
          escape = FALSE,
          extensions = 'Buttons',
          rownames = FALSE,
          options = list(
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autowidth = TRUE,
            ordering = TRUE,
            dom = 'Bflrtip',
            scrollX = TRUE,
            lengthMenu = list(c(20, 50,100, -1), c('20', '50', '100', 'All')),
            buttons = list(
              list(extend = "excel", text = "Download Table", title = "", filename = "transcriptR_GO_GSEA_page", 
                   exportOptions = list(
                     modifier = list(page = "current")
                   )
              )              
            )
          ), class = "display"
        )}
    #}
    
  })
#####============== Gene Ontology module end==============#######

#####============== heatmap module start==============#######

file1heat = reactive({
    fileheat1 <- input$file1heat
    ext <- tools::file_ext(fileheat1$datapath)
    req(fileheat1)
        
    fileHeatpath <- fileheat1$datapath
    fileheat1 <- switch(ext,
               txt = readr::read_table(fileHeatpath),
               tsv = readr::read_tsv(fileHeatpath),
               csv = readr::read_csv(fileHeatpath),
               xls = readxl::read_xls(fileHeatpath),
               xlsx = readxl::read_xlsx(fileHeatpath))

    # fileheat1 <- read.delim(fileheat1$datapath,
    #                         sep = ",",
    #                         header = TRUE,
    #                         row.names = 1)
    
    fileheat11 <- as.matrix(fileheat1)
    rownames(fileheat11) = rownames(fileheat1)

    return(fileheat11)
  
  }
  )
    
  output$heatmapdatafile <- DT::renderDataTable( 
    file1heat(), rownames = FALSE,
    filter = "bottom",
        style = "auto",
        selection = "single",
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    )
  )
  
  file2heat = reactive({
    fileheat2 <- input$file2heat
    ext <- tools::file_ext(fileheat2$datapath)
    req(fileheat2)
    
    fileheatpath2 <- fileheat2$datapath
    fileheat2 <- switch(ext,
               txt = readr::read_table(fileheatpath2),
               tsv = readr::read_tsv(fileheatpath2),
               csv = readr::read_csv(fileheatpath2),
               xls = readxl::read_xls(fileheatpath2),
               xlsx = readxl::read_xlsx(fileheatpath2))

    # fileheat2 <- read.delim(fileheat2$datapath,
    #                          sep = ",", 
    #                          header = T)

    fileheat21 <- fileheat2
  
    updateSelectInput(session, "annotationval1", choices = names(fileheat21))
    updateSelectInput(session, "annotationval2", choices = names(fileheat21))
    updateSelectInput(session, "annotationval3", choices = names(fileheat21))
    updateSelectInput(session, "annotationval4", choices = names(fileheat21))
    updateSelectInput(session, "groupannoval1", choices = names(fileheat21))
    updateSelectInput(session, "groupannoval2", choices = names(fileheat21))
    return(fileheat2)
  }
  )
  output$heatmapannotationfile <- DT::renderDataTable(
    file2heat(), rownames = FALSE,
    filter = "bottom",
        style = "auto",
        selection = "single",
    options = list(
      autowidth = TRUE,
      lengthMenu = list(c(10, 20, 50, -1), c('10', '20', '50', 'All')),
      searching = TRUE
    )
    
  )

observeEvent(input$legeparams, {
  showModal(modalDialog(
    title = "Select Legend Parameters",
    tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
    selectInput(
                "legendside",
                label = "Select Legend Side",
                choices = c(
                          "Right" = "right",
                          "Left" = "left",
                          "Top" = "top",
                          "Bottom" = "bottom"
                        )
                      ),
                      textInput("legetitle",
                                label = "Enter Legend Title",
                                value = "bVals"),
                      selectInput("legedir",
                                  label = "Select Legend Direction",
                                  choices = c(
                                    "Vertical" = "vertical",
                                    "Horizontal" = "horizontal"
                                  )),
                      selectInput("legetitlpos",
                                  label = "Legend Title Position",
                                  choices = c(
                                    "Top Left" = "topleft",
                                    "Top Center" = "topcenter",
                                    "Left Center" = "leftcenter",
                                    "Left Top" = "lefttop",
                                    "Left Center Rotated" = "leftcenter-rot",
                                    "Left Top Rotated" = "lefttop-rot"
                                  )),
                      numericInput("legelabrot",
                                  label = "Legend Label Rotation",
                                  value = 0),
                      selectInput("legetype",
                                  label = "Legend Type",
                                  choices = c(
                                    "Grid" = "grid",
                                    "Points" = "points",
                                    "Lines" = "lines",
                                    "Boxplot" = "boxplot"
                                  )),
                      numericInput("legefontsize",
                                  label = "Legend Font Size",
                                  value = 12)
  ),
              easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitlegeparams", "OK")
              )
  ))
})

observeEvent(input$submitlegeparams, {
  removeModal()
})

#Group-wise Annotation-1 modal
observeEvent(input$heatmapannoparams1, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
  textInput("groupannoname1",
            label = "Change the Group-wise Annotation Name",
            value = "Gender"),
  colourInput("groupannocolo1",
            h5("Select Group-1 First Element color"),
            value = "#B0C4DE"), 
  colourInput("groupannocolo2",
            h5("Select Group-1 Second Element color"),
            value = "#E6E6FA"),
  numericInput("annoheight",
              label = "Choose Annotation Height",
              min = 3,
              value = 4.5),
  selectInput("annoheightunit",
              label = "Select Annotation Height Unit",
              choices = c(
                "mm" = "mm",
                "cm" = "cm",
                "inch" = "in"
              ),
              selected = "cm",
              multiple = FALSE)),
  
  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams1", "OK")
              )
  ))
})

observeEvent(input$submitannotationparams1, {
  removeModal()
})

#Group-wise Annotation-2 modal
observeEvent(input$heatmapannoparams2, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
  
  textInput("groupannoname2",
              label = "Change the Second Group-wise Annotation Name",
              value = "Sample Group"),
    colourInput("groupannocolo3",
                h4("Select Group-2 First Element color"),
                value = "#C0C0C0"),
    colourInput("groupannocolo4",
                h4("Select Group-2 Second Element color"),
                value = "#B8860B"),
    colourInput("groupannocolo5",
                h4("Select Group-2 Third Element color"),
                value = "#FFD700"),
  numericInput("annoheight",
              label = "Choose Annotation Height",
              min = 3,
              value = 4.5),
  selectInput("annoheightunit",
              label = "Select Annotation Height Unit",
              choices = c(
                "mm" = "mm",
                "cm" = "cm",
                "inch" = "in"
              ),
              selected = "cm",
              multiple = FALSE),
  
  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams2", "OK")
              )
  ))
})

observeEvent(input$submitannotationparams2, {
  removeModal()
})

#Numer-based Annotation-1 modal
observeEvent(input$heatmapannoparams3, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
    textInput("annoval1name",
            label = "Change the Annotation Name",
            value = "Age"),
    selectInput(
      "annopointtype",
      "Select Point Types",
      choices = c(
        "Square" = "15",
        "Triangle" = "17",
        "Circle" = "19",
        "Diamond" = "18"
      ),
      selected = "19",
      multiple = FALSE
    ), 
    numericInput(
      "annopointsize",
      "Choose Annotation Point Size",
      value = 2
    ),
    colourInput(
      "annopointcolo1",
      "Choose Point Color",
      value = "#3680ac"
    ),

  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams3", "OK")
              )
    )  )
})

observeEvent(input$submitannotationparams3, {
  removeModal()
})


#Numer-based Annotation-2 modal
observeEvent(input$heatmapannoparams4, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
    textInput("annoval2name",
          label = "Change the Annotation Name",
          value = "BMI"),
    selectInput(
      "annopointtype2",
      "Select Point Types",
      choices = c(
        "Square" = "15",
        "Triangle" = "17",
        "Circle" = "19",
        "Diamond" = "18"
      ),
      selected = "19",
      multiple = FALSE
    ), 
    numericInput(
      "annopointsize2",
      "Choose Annotation Point Size",
      value = 2
    ),
    colourInput(
      "annopointcolo2",
      "Choose Point Color",
      value = "#3680ac"
    ),

  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams4", "OK")
              )
    )  )
})

observeEvent(input$submitannotationparams4, {
  removeModal()
})

#Numer-based Annotation-3 modal
observeEvent(input$heatmapannoparams5, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
    textInput("annoval3name",
          label = "Change the Annotation Name",
          value = "CRP"),
    selectInput(
      "annopointtype3",
      "Select Point Types",
      choices = c(
        "Square" = "15",
        "Triangle" = "17",
        "Circle" = "19",
        "Diamond" = "18"
      ),
      selected = "19",
      multiple = FALSE
    ), 
    numericInput(
      "annopointsize3",
      "Choose Annotation Point Size",
      value = 2
    ),
    colourInput(
      "annopointcolo3",
      "Choose Point Color",
      value = "#3680ac"
    ),

  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams5", "OK")
              )
    )  )
})

observeEvent(input$submitannotationparams5, {
  removeModal()
})


#Numer-based Annotation-4 modal
observeEvent(input$heatmapannoparams6, {
  showModal(modalDialog(
  title = "Select Annotation Parameters",
    textInput("annoval4name",
          label = "Change the Annotation Name",
          value = "LDH"),
    selectInput(
      "annopointtype4",
      "Select Point Types",
      choices = c(
        "Square" = "15",
        "Triangle" = "17",
        "Circle" = "19",
        "Diamond" = "18"
      ),
      selected = "19",
      multiple = FALSE
    ), 
    numericInput(
      "annopointsize4",
      "Choose Annotation Point Size",
      value = 2
    ),
    colourInput(
      "annopointcolo4",
      "Choose Point Color",
      value = "#3680ac"
    ),

  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitannotationparams6", "OK")
              )
    )  )
})

observeEvent(input$submitannotationparams6, {
  removeModal()
})

#HeatMap Clustering modal
observeEvent(input$clusterheatmap, {
  showModal(modalDialog(
  title = "Select Clustering Parameters",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
  h4("Column-wise clustering"),
  checkboxInput("colmdend",
                label = "Show/Hide Column Dendrogram",
                value = TRUE,
                width = "400px"
                ),
  selectInput("colmdend1",
                label = "Place Column Dendrogram",
                choices = c(
                  "Top" = "top",
                  "Bottom" = "bottom"
                ),
                selected = "top",
                multiple = FALSE),
  h4("Row-wise clustering"),
  checkboxInput("rowdend",
                label = "Show/Hide Row Dendrogram",
                value = TRUE,
                width = "400px"
                ),
    selectInput("rowdend1",
                label = "Place Row Dendrogram",
                choices = c(
                  "Right" = "right",
                  "Left" = "left"
                ),
                selected = "left",
                multiple = FALSE)),
    # numericInput("rowclustnum",
    #             "Select number of cluster",
    #             value = 1),

  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitclusterheatmap", "OK")
              )
    )  )
})

observeEvent(input$submitclusterheatmap, {
  removeModal()
})

# Main parameters Modal
observeEvent(input$heatmapmainparams, {
  showModal(modalDialog(
  tabsetPanel(
    id = "selpar1",
    tabPanel(
      title = "Select Parameters",
        tags$head(
            tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                      .inline .form-group{display: table-row;}")
            ),
          
        tags$div(class = "inline",
        numericInput(
                    "colmrothm",
                    h5("Column name rotation"),
                    value = 45,
                    min = 0,
                    max = 360,
                    step = 1
                  ),
        numericInput(
                    "colonumlowheat",
                    h5("Select Lowest Value"),
                    value = 0
                  ),
        numericInput(
                    "colonummidheat",
                    h5("Select Mid Value"),
                    value = 0.5
                  ),
        numericInput(
                    "colonumhighheat",
                    h5("Select Highest Value"),
                    value = 1
                  ),
        numericInput(
                    "colmclusternum",
                    h5("Select Number of Clusters for Columns"),
                    value = 1
                  ),
        numericInput(
                    "heatmapsizeheight",
                    h5("Select Height of the HeatMap"),
                    value = 20
                  ),
        numericInput(
                    "heatmapsizewidth",
                    h5("Select Width of the HeatMap"),
                    value = 25
                  ),
        selectInput(
                    "heatmapsizeunit",
                    h5("Select unit for Heatmap Size"),
                    choices = c(
                      "cm", "mm", "in"
                    ),
                    selected = "cm",
                    multiple = FALSE
                  )),
        br(),
        br(),
        actionButton("submitheatmapmainparams", "Submit",
        icon("paper-plane"), 
                        style="color: #000000; background-color: #fcffa4; border-color: #2e6da4")
        
        ),


    tabPanel(
    title = "Parameters Guide",
    h5(HTML("1. <b><i>Column name rotation:</b></i> Select the angle of rotation for the column name, eg. 90"),
    br(),
    br(),
    HTML("2. <b><i>Select Lowest Value:</b></i> Select minimum value for heatmap; default is 0"),
    br(),
    br(),
    HTML("3. <b><i>Select Mid Value:</b></i> Select mid value for heatmap; default is 0.5"),
    br(),
    br(),
    HTML("4. <b><i>Select Highest Value:</b></i> Select maximum value for heatmap; default is 1"),
    br(),
    br(),
    HTML("5. <b><i>Select Number of Clusters for Columns:</b></i> Cluster the heatmap based on column, increase the number; default is 1."),
    br(),
    br(),
    HTML("6. <b><i>Select Height of the Heatmap:</b></i> Select height of the heatmap; default is 20"),
    br(),
    br(),
    HTML("7. <b><i>Select Width of the Heatmap:</b></i> Select width of the heatmap; default is 25"),
    br(),
    br(),
    HTML("8. <b><i>Select unit for Heatmap size:</b></i> Select length unit for heatmap; default is centimeters, cm."),
    ),
    br(),
    br(),
    HTML("<b><i>For more details, please see the manual.</b></i>"),
  
      )),

  easyClose = TRUE,
         footer = NULL
              
    ))
})



observeEvent(input$submitheatmapmainparams, {
  removeModal()
})

# Row Annotation Parameters
observeEvent(input$rowannoparams, {
  showModal(modalDialog(
  title = "Select Parameters",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
  numericInput(
              "rowannonum",
              h4("Number for Row Annotation"),
              value = 10,
              min = 0,
              max = 100,
              step = 10
            ),
  # The annotation will be placed in position if the row cluster removed
  numericInput(
              "rowannonum1",
              h5("Set-1: Add First Row Number"),
              value = 101,
              min = 0
            ),
  numericInput(
              "rowannonum2",
              h5("Set-1: Add Last Row Number"),
              value = 110,
              min = 0
            ),
    numericInput(
              "rowannonum3",
              h5("Set-2: Add First Row Number"),
              value = 301,
              min = 0
            ),
      numericInput(
              "rowannonum4",
              h5("Set-2: Add Last Row Number"),
              value = 310,
              min = 0
            ),
      numericInput(
              "rowannonum5",
              h5("Add Single row line"),
              value = 501,
              min = 0
            ),
  numericInput(
              "rowannosplitnum",
              h4("Select Number to split rows"),
              value = 1,
              min = 1
            ),
  numericInput(
              "rowannosplitgap",
              h4("Select distance between slices"),
              value = 0,
              min = 0
            )
  ),
  
  easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitrowannoparams", "OK")
              )
            )
          )
})

observeEvent(input$submitrowannoparams, {
  removeModal()
})

# Heatmap Colour Settings
observeEvent(input$coloheatmap, {
  showModal(modalDialog(
  title = "Select Parameters",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
  colourInput("highcoloheat",
              h4(" Select higher value color"),
              value = "firebrick3"),
    colourInput("lowcoloheat",
              h4(" Select lower value color"),
              value = "dodgerblue3"),
    colourInput("midcoloheat",
              h4(" Select middle value color"),
              value = "white"),
    colourInput("colmcoloheat",
              h4(" Select Column Font color"),
              value = "black"),
    colourInput("rowcoloheat",
              h4(" Select Row Font color"),
              value = "black"),
    colourInput("colmtitlecoloheat",
              h4(" Select Column Title Font color"),
              value = "black"),
    colourInput("rowtitlecoloheat",
              h4(" Select Row Title Font color"),
              value = "black"),
    colourInput("legetitlecoloheat",
              h4(" Select Legend Title Font color"),
              value = "black"),
    colourInput("legelablecoloheat",
              h4(" Select Legend Labels Font color"),
              value = "black")
  ),
    easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitcolosettparams", "OK")
              )
  )
  )
})

observeEvent(input$submitcolosettparams, {
  removeModal()
})

## Font Settings
observeEvent(input$fontheatmap, {
  showModal(modalDialog(
  title = "Select Parameters",
  tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
  selectInput(
    "colmnamefontfam",
    h5("Select Column Names Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  ), 
  selectInput(
    "rownamefontfam",
    h5("Select Row Names Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  ),
  selectInput(
    "colmtitlfontfam",
    h5("Select Column Title Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  ),
  selectInput(
    "rowtitlfontfam",
    h5("Select Row Title Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  ),
  selectInput(
    "legetitlfontfam",
    h5("Select Legend Title Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  ),
  selectInput(
    "legelabfontfam",
    h5("Select Lagend Labels Font Family"),
    choices = c(
      "Arial" = "sans",
      "Times New Roman" = "Times",
      "Courier" = "mono"
    ),
  multiple = FALSE,
  selected = "Times"
  )
  ),

    easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitcolosettparams", "OK")
              )
  )
  )
})

observeEvent(input$submitcolosettparams, {
  removeModal()
})

# Font Size Settings
observeEvent(input$fontsizeheatmap, {
  showModal(modalDialog(
    title = "Select Font Size Parameters",
    tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
    numericInput(
          "colmfontsizeheatmap",
          h5("Column Name Font Size"),
          value = 12
          ),
        numericInput(
          "rowfontsizeheatmap",
          h5("Row Name Font Size"),
          value = 12
          ),
        numericInput(
          "colmtitlfontsizeheatmap",
          h5("Column Title Font Size"),
          value = 12
          ),
        numericInput(
          "rowtitlfontsizeheatmap",
          h5("Row Title Font Size"),
          value = 12
          ),
        numericInput(
          "legetitlfontsizeheatmap",
          h5("Other Legend Title Font Size"),
          value = 12
          ),
        numericInput(
          "legelabfontsizeheatmap",
          h5("Other Legend Label Font Size"),
          value = 12
              )
  ),
easyClose = TRUE,
              footer = tagList(
                modalButton("Cancel"),
                actionButton("submitfontsizeheatmap", "OK")
              )
  ))
})

observeEvent(input$submitfontsizeheatmap, {
  removeModal()
})

# Font Alpha Settings
observeEvent(input$fontalphaheatmap, {
  showModal(modalDialog(
    title = "Select Legend Parameters",
    tags$head(
      tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                 .inline .form-group{display: table-row;}")
      ),
    
  tags$div(class = "inline",
    numericInput(
          "colonamealphahm",
          label = "Change Column Name Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "rownamealphahm",
          label = "Change Row Name Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "colotitlalphahm",
          label = "Change Column Title Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "rowtitlalphahm",
          label = "Change Row Title Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "legetitlalphahm",
          label = "Change Legend Title Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "legelablalphahm",
          label = "Change Legend Label Color Transparency",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        )
  ),

easyClose = TRUE,
footer = tagList(
  modalButton("Cancel"),
  actionButton("submitfontalphaheatmap", "OK")
              )
  ))
})

observeEvent(input$submitfontalphaheatmap, {
  removeModal()
})


#plot
heatmplot <- eventReactive(input$submitheatmap, {
id_heatmapplot1 <- showNotification("generating plot, please wait...", duration = NULL, closeButton = FALSE)
on.exit(removeNotification(id_heatmapplot1), add = TRUE)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(magick))

col_fun = circlize::colorRamp2(c(input$colonumlowheat,input$colonummidheat,input$colonumhighheat), 
            c(input$lowcoloheat, 
            input$midcoloheat, 
            input$highcoloheat))

ht_opt(heatmap_column_names_gp = gpar(fontfamily = input$colmnamefontfam,
                                      fontsize = input$colmfontsizeheatmap,
                                      col = input$colmcoloheat,
                                      alpha = input$colonamealphahm), 
       heatmap_row_names_gp = gpar(fontfamily = input$rownamefontfam,
                                      fontsize = input$rowfontsizeheatmap,
                                      col = input$rowcoloheat,
                                      alpha = input$rownamealphahm), 
       heatmap_column_title_gp = gpar(fontfamily = input$colmtitlfontfam,
                                      fontsize = input$colmtitlfontsizeheatmap,
                                      col = input$colmtitlecoloheat,
                                      alpha = input$colotitlalphahm),
       heatmap_row_title_gp = gpar(fontfamily = input$rowtitlfontfam,
                                      fontsize = input$rowtitlfontsizeheatmap,
                                      col = input$rowtitlecoloheat,
                                      alpha = input$rowtitlalphahm),
       legend_title_gp = gpar(fontfamily = input$legetitlfontfam,
                                      fontsize = input$legetitlfontsizeheatmap,
                                      col = input$legetitlecoloheat,
                                      alpha = input$legetitlalphahm),
       legend_labels_gp = gpar(fontfamily = input$legelabfontfam,
                                      fontsize = input$legelabfontsizeheatmap,
                                      col = input$legelablecoloheat,
                                      alpha = input$legelablalphahm))

ha = rowAnnotation(foo = anno_mark(at=c(input$rowannonum1:input$rowannonum2,
                                      input$rowannonum3:input$rowannonum4,
                                      input$rowannonum5), 
                   labels = rownames(file1heat())[1:input$rowannonum]))

if(input$annotationheatmap == TRUE & input$annoagain1 == FALSE & input$annoagain2 == FALSE & input$annoagain3 == FALSE){
.annoval1 = anno_lines(file2heat()[,input$annotationval1], 
                      add_points = TRUE,
                      pch = as.numeric(input$annopointtype), 
                      size = unit(input$annopointsize, "mm"),                      
                      pt_gp = gpar(col = input$annopointcolo1)
                      )

ha2 = HeatmapAnnotation(val1 = .annoval1,
                        height = unit(input$annoheight, input$annoheightunit),
                        anngrp1 = file2heat()[,input$groupannoval1], # Group1
                        anngrp2 = file2heat()[,input$groupannoval2], #Group2
                        show_legend = TRUE,
                              col = list(
                                anngrp1 = c("M" = input$groupannocolo1,
                                           "F" = input$groupannocolo2),
                                anngrp2 = c("CON" = input$groupannocolo3, 
                                           "PT1" = input$groupannocolo4,
                                           "PT2" = input$groupannocolo5)
                                         )
                        )

ha2@anno_list$val1@label <- input$annoval1name
ha2@anno_list$anngrp1@label <- input$groupannoname1
ha2@anno_list$anngrp2@label <- input$groupannoname2
ha2@anno_list$anngrp1@color_mapping@name <- input$groupannoname1
ha2@anno_list$anngrp2@color_mapping@name <- input$groupannoname2

} else if(input$annotationheatmap == TRUE & input$annoagain1 == TRUE & input$annoagain2 == FALSE & input$annoagain3 == FALSE){
.annoval1 = anno_lines(file2heat()[,input$annotationval1], 
                      add_points = TRUE,
                      pch = as.numeric(input$annopointtype), 
                      size = unit(input$annopointsize, "mm"),
                      pt_gp = gpar(col = input$annopointcolo1)
                      )
.annoval2 = anno_lines(file2heat()[,input$annotationval2], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype2), 
                      size = unit(input$annopointsize2, "mm"),
                      pt_gp = gpar(col = input$annopointcolo2)
                      )

ha2 = HeatmapAnnotation(val1 = .annoval1, 
                        val2 = .annoval2,
                        height = unit(input$annoheight, input$annoheightunit),
                        anngrp1 = file2heat()[,input$groupannoval1], # Group1
                        anngrp2 = file2heat()[,input$groupannoval2], #Group2
                        show_legend = TRUE,
                              col = list(
                                anngrp1 = c("M" = input$groupannocolo1,
                                           "F" = input$groupannocolo2),
                                anngrp2 = c("HC" = input$groupannocolo3, 
                                           "PT1" = input$groupannocolo4,
                                           "PT2" = input$groupannocolo5)
                                         ))

ha2@anno_list$val1@label <- input$annoval1name
ha2@anno_list$val2@label <- input$annoval2name
ha2@anno_list$anngrp1@label <- input$groupannoname1
ha2@anno_list$anngrp2@label <- input$groupannoname2
ha2@anno_list$anngrp1@color_mapping@name <- input$groupannoname1
ha2@anno_list$anngrp2@color_mapping@name <- input$groupannoname2

} else if(input$annotationheatmap == TRUE & input$annoagain1 == TRUE & input$annoagain2 == TRUE & input$annoagain3 == FALSE) {
.annoval1 = anno_lines(file2heat()[,input$annotationval1], 
                      add_points = TRUE,
                      pch = as.numeric(input$annopointtype), 
                      size = unit(input$annopointsize, "mm"),
                      pt_gp = gpar(col = input$annopointcolo1)
                      )
.annoval2 = anno_lines(file2heat()[,input$annotationval2], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype2), 
                      size = unit(input$annopointsize2, "mm"),
                      pt_gp = gpar(col = input$annopointcolo2)
                      )
.annoval3 = anno_lines(file2heat()[,input$annotationval3], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype3), 
                      size = unit(input$annopointsize3, "mm"),
                      pt_gp = gpar(col = input$annopointcolo3)
                      )

ha2 = HeatmapAnnotation(val1 = .annoval1, 
                        val2 = .annoval2,
                        val3 = .annoval3,
                        height = unit(input$annoheight, input$annoheightunit),
                        anngrp1 = file2heat()[,input$groupannoval1], # Group1
                        anngrp2 = file2heat()[,input$groupannoval2], #Group2
                        show_legend = TRUE,
                              col = list(
                                anngrp1 = c("M" = input$groupannocolo1,
                                           "F" = input$groupannocolo2),
                                anngrp2 = c("PAT" = input$groupannocolo3, 
                                           "CON" = input$groupannocolo4)
                                         ))

ha2@anno_list$val1@label <- input$annoval1name
ha2@anno_list$val2@label <- input$annoval2name
ha2@anno_list$val3@label <- input$annoval3name

ha2@anno_list$anngrp1@label <- input$groupannoname1
ha2@anno_list$anngrp2@label <- input$groupannoname2
ha2@anno_list$anngrp1@color_mapping@name <- input$groupannoname1
ha2@anno_list$anngrp2@color_mapping@name <- input$groupannoname2

} else {
.annoval1 = anno_lines(file2heat()[,input$annotationval1], 
                      add_points = TRUE,
                      pch = as.numeric(input$annopointtype), 
                      size = unit(input$annopointsize, "mm"),
                      pt_gp = gpar(col = input$annopointcolo1)
                      )
.annoval2 = anno_lines(file2heat()[,input$annotationval2], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype2), 
                      size = unit(input$annopointsize2, "mm"),
                      pt_gp = gpar(col = input$annopointcolo2)
                      )
.annoval3 = anno_lines(file2heat()[,input$annotationval3], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype3), 
                      size = unit(input$annopointsize3, "mm"),
                      pt_gp = gpar(col = input$annopointcolo3)
                      )
.annoval4 = anno_lines(file2heat()[,input$annotationval4], 
                      add_points = TRUE, 
                      pch = as.numeric(input$annopointtype4), 
                      size = unit(input$annopointsize4, "mm"),
                      pt_gp = gpar(col = input$annopointcolo4)
                      )


ha2 = HeatmapAnnotation(val1 = .annoval1, 
                        val2 = .annoval2,
                        val3 = .annoval3,
                        val4 = .annoval4,
                        height = unit(input$annoheight, input$annoheightunit),
                        anngrp1 = file2heat()[,input$groupannoval1], # Group1
                        anngrp2 = file2heat()[,input$groupannoval2], #Group2
                        show_legend = TRUE,
                              col = list(
                                anngrp1 = c("M" = input$groupannocolo1,
                                           "F" = input$groupannocolo2),
                                anngrp2 = c("HC" = input$groupannocolo3, 
                                           "PT1" = input$groupannocolo4,
                                           "PT2" = input$groupannocolo5
                                           )
                                         ))

ha2@anno_list$val1@label <- input$annoval1name
ha2@anno_list$val2@label <- input$annoval2name
ha2@anno_list$val3@label <- input$annoval3name
ha2@anno_list$val4@label <- input$annoval4name

ha2@anno_list$anngrp1@label <- input$groupannoname1
ha2@anno_list$anngrp2@label <- input$groupannoname2
ha2@anno_list$anngrp1@color_mapping@name <- input$groupannoname1
ha2@anno_list$anngrp2@color_mapping@name <- input$groupannoname2
}

library(dendextend)
rowdendro = as.dendrogram(hclust(dist(file1heat())))
rowdendro = color_branches(rowdendro, k = input$rowclustnum)

ht_list = Heatmap(file1heat(), 
                        col = col_fun,
                        show_heatmap_legend = FALSE,  
                        column_names_rot = input$colmrothm,
                        top_annotation = ha2,
                        show_column_dend = input$colmdend,
                        column_dend_side = input$colmdend1,
                        cluster_rows = input$rowdend,
                        row_dend_side = input$rowdend1,
                        row_km = input$rowannosplitnum,
                        row_gap = unit(input$rowannosplitgap, "mm"),
                        column_km = input$colmclusternum,
                        width = unit(input$heatmapsizewidth, input$heatmapsizeunit),
                        height = unit(input$heatmapsizeheight, input$heatmapsizeunit)
) + ha 

lgd = Legend(col_fun = col_fun, 
             title = input$legetitle, 
             at = c(input$colonumlowheat,input$colonummidheat/2, input$colonummidheat,(input$colonumhighheat + input$colonummidheat)/2,input$colonumhighheat),
             labels_rot = input$legelabrot,
             direction = input$legedir,
             type = input$legetype,
             labels_gp = gpar(fontsize = input$legefontsize),
             title_position = input$legetitlpos)

hm <- draw(ht_list, annotation_legend_list = lgd, annotation_legend_side = input$legendside)
                     
return(hm)
}
)

output$heatmapplot <- renderPlot(
  heatmplot()
)

output$heatDownload <- downloadHandler(
          filename = function() {
            paste("HeatMap", Sys.Date(), tolower(input$filetype_heat), sep = ".")},
          content = function(file) {
          if(input$filetype_heat == "PNG") 
          {png(file, 
            width = input$heat_pngwidth , 
            height = input$heat_pngheight, 
            units = input$heat_pngunits, 
            res = input$heat_pngresol, 
            type = "cairo")} 
          else if(input$filetype_heat == "PDF")
          {pdf(file,
            width = input$heat_pdfwidth , 
            height = input$heat_pdfheight
            )}
          else if(input$filetype_heat == "SVG")
          {svg(file, 
            width = input$heat_svgwidth , 
            height = input$heat_svgheight
            )}
          else (tiff(file, 
            width = input$heat_tiffwidth , 
            height = input$heat_tiffheight, 
            units = input$heat_tiffunits, 
            res = input$heat_tiffresol,
          type = "cairo"))
            

          draw(heatmplot())

          dev.off()
          }
          )
#####============== Heatmap module end==============#######

  }

shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))