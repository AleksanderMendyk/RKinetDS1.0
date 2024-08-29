##RKinetDS v1.2 - software for modeling dissolution profiles
##Authors: 
#Natalia Obajtek: nat.obajtek@gmail.com
#Aleksander Mendyk: mfmendyk@cyf-kr.edu.pl; aleksander.mendyk@uj.edu.pl
#Jakub Szlęk: j.szlek@uj.edu.pl
#Adam Pacławski: adam.paclawski@uj.edu.pl
##License: GPLv3

library(shiny)
library(shinyjs)
library(shinyFiles)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(yaml)
library(plotly)
library(stringr)
library(DT)

########################################### INITIAL ########################################### 

#Read main engine function
source('RKinetDS_computational_core.R')

cfg <- yaml::yaml.load_file('user_config.yml')
cfg_name <- "user_config.yml"

#CSS background color
backgr_color_css <-"
/* Select the sidebar by its class and set its background color */
.well {
background-color: #FDFDFD;
  }
"

#CSS clicked button color
check_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button to blue. */
background: DodgerBlue;
/* Change the text size to 15 pixels. */
font-size: 15px;
}"


########################################### USER INTERFACE ###########################################

ui <- shiny::fluidPage(
  useShinyjs(), # useShinyjs is necessary to enable/disable elements of UI
  theme = shinytheme("lumen"),
  tags$style(backgr_color_css),
  tags$style(check_color_css),

  #Navbar structure for UI
  navbarPage("RKinetDS",
             tabPanel("Main panel", fluid = TRUE, icon = icon("glyphicon glyphicon-pencil", lib="glyphicon"),
                      titlePanel("Main panel"),
                      sidebarLayout(
                        sidebarPanel(
                          width = 4,

                          #graph
                          p(strong("Upload your data:")),
                          shinyFilesButton(id = "file_data",
                                           label = "Choose" ,
                                           title = "Choose a data file (tab limited):",
                                           multiple = FALSE),
                          br(),
                          p("Name of a selected file:"),
                          verbatimTextOutput("chosen_file"),
                          br(),
                          p(strong("Click to start computation:")),
                          actionButton("run", "Run app")
                        ),

                        mainPanel(
                          plotlyOutput("plot")
                        )
                      ),
             ),

             tabPanel("Settings", fluid = TRUE, icon = icon("glyphicon glyphicon-play", lib = "glyphicon"),
                      titlePanel("Settings"),
                      sidebarLayout(
                        sidebarPanel(
                          width = 4,
                          radioButtons("selection_options", strong("Selection options"),
                                       choices = c("Select all models", "Deselect all models","Restore standard settings"),
                                       selected = 0
                          ),

                          actionButton("save", "Save settings")
                        ),

                        mainPanel(

                          tabsetPanel(type = "tabs",
                                      tabPanel("Models",
                                               br(),

                                               checkboxGroupInput(inputId = "selected_models",
                                                                  label = (strong("Models list:")),
                                                                  choices = cfg$models,
                                                                  selected = cfg$models[unlist(cfg$selected_modelsx$selected_models)])
                                      ),

                                      tabPanel("Optimization methods",
                                               br(),
                                               checkboxGroupInput(inputId = "selected_optimmethods",
                                                                  label = (strong("Optimization method:")),
                                                                  choices = names(cfg$optim_method),
                                                                  selected = ifelse((cfg$optim_method)==TRUE, (names(cfg$optim_method==TRUE)), 0),
                                               ),

                                               br(),
                                               p(strong("Optimization parameters:
                                   ")),

                                               selectInput(inputId = "selected_optitrace",
                                                           label = "Tracing of optimization function evaluations:",
                                                           choices = c("Yes", "No"),
                                                           selected = (if (TRUE %in% cfg$running_params$opti_trace) "Yes" else "No")

                                               ),

                                               textInput(inputId = "BFGS_params",
                                                         label = "Maximum number of iterations in BFGS method:",
                                                         value = ifelse(length(cfg$running_params$maxit_BFGS)==0, 5000, cfg$running_params$maxit_BFGS)),

                                               textInput(inputId = "rel_tol_params",
                                                         label = "Value of stop criterion for optimizing:",
                                                         value = ifelse(length(cfg$running_params$optim_rel_tol)==0, 1e-20, cfg$running_params$optim_rel_tol)),

                                               uiOutput("optim_methods")


                                      ),
                                      tabPanel("Data format",
                                               br(),

                                               selectInput(inputId = "selected_headers",
                                                           label = (strong("Headers in input data:")),
                                                           choices = c("Yes", "No"),
                                                           selected = (if (TRUE %in% cfg$data$headers) "Yes" else "No")

                                               ),
                                      )),
                        )
                      ),),
#
             navbarMenu("Results", icon = icon("glyphicon glyphicon-stats", lib = "glyphicon", id = "overall_tab"),
                        tabPanel("Used settings", fluid = TRUE,
                                 p(strong("Used models:")), verbatimTextOutput("used_models_rep"),
                                 p(strong("Used optimization methods:")), verbatimTextOutput("used_optim_rep"),
                                 p(strong("Used optimization parameters:")), verbatimTextOutput("used_param_rep"),
                                 p(strong("Used data format:")), verbatimTextOutput("used_data_rep")
                        ),
                        tabPanel("Overall results", fluid = TRUE, verbatimTextOutput("display_report")),
                        tabPanel("Model results", uiOutput("tabs_with_models_results")),
                        tabPanel("Error ranking",  
                                 fluidRow(
                                   column(6, dataTableOutput("table_error_RMSE")),
                                 ),
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(6, dataTableOutput("table_error_R2")),
                                 ),
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(6, dataTableOutput("table_error_R2_ad")),
                                 ),
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(6, dataTableOutput("table_error_AIC")),
                                 ),
                                 
                                 
                        
                        ),
                        ),


             navbarMenu("More", icon = icon("glyphicon glyphicon-info-sign", lib="glyphicon"),
                        # tabPanel("About", fluid = TRUE, verbatimTextOutput("help")
                        # ),
                        tabPanel("Manual", fluid = TRUE, htmlOutput("manual")
                        ),

                        tabPanel("License", fluid = TRUE, htmlOutput("license")
                        ),
             ),

  )

)

########################################### SERVER ###########################################

server <- function(input, output, session) {

  # Yaml config file reader (reactive)
  user_config_reader <- reactiveFileReader(intervalMillis = 5000, session, filePath =
                                 "user_config.yml", readFunc = yaml::yaml.load_file)
  
  # Functions to check logical values from UI - getters
  #
  #Opti trace
  get_optitrace <- function(selected_optitrace, optitrace, session) {
    input <- session$input
    if (input[[selected_optitrace]] == "Yes") {
      as.logical(1)
    } else {as.logical(0)}
  }
  
  #Headers
  get_headers <- function(selected_headers, header, session) {
    input <- session$input
    if (input[[selected_headers]] == "Yes") {
      as.logical(1)
    } else {as.logical(0)}
  }
  
  #Choice of optimization methods
  get_optim <- function(selected_optimmethods, choice_name, session) {
    input <- session$input
    if (choice_name %in% input[[selected_optimmethods]]) {
      as.logical(1)
    } else {as.logical(0)}
  }
  

  # save_config function - it does not take any arguments
  save_config <- function(session){

    # # Yaml user_config file ### TODO - make tryCatch work properly
    # tryCatch(
    #   current_user_config <- yaml::yaml.load_file('./user_config.yml'),
    #   finally = {
    #     current_user_config <- NULL
    #     }
    # )
    current_user_config <- yaml::yaml.load_file('user_config.yml')

    # read session input
    input <- session$input
    
    # Yaml factory config file
    factory_cfg <- yaml::yaml.load_file('factory_config.yml')

    # Choice of models
    all_used_models <- cfg$models
    chosen_models <- input$selected_models
    used_models <- which(all_used_models %in% chosen_models)
    models_to_use <- unname(used_models)
    
    #File name
    fileName <- uploadedFileName()
    
    # Filling yaml_content by getters + other objects/variables
    config$yaml_content$optim_method$SANN <- get_optim("selected_optimmethods", "SANN", session)
    config$yaml_content$optim_method$nloptr <- get_optim("selected_optimmethods", "nloptr", session)
    config$yaml_content$optim_method$NelderMead <- get_optim("selected_optimmethods", "NelderMead", session)
    config$yaml_content$optim_method$genSA <- get_optim("selected_optimmethods", "genSA", session)
    config$yaml_content$optim_method$rgenoud <- get_optim("selected_optimmethods", "rgenoud", session)
    config$yaml_content$running_params$opti_trace <- get_optitrace("selected_optitrace", "Yes", session)
    config$yaml_content$data$headers <- get_headers("selected_headers", "Yes", session)
    
    config$yaml_content$selected_modelsx$selected_models <- models_to_use
    config$yaml_content$data$input_data <- fileName
    
    # Optimization parameters
    ### TODO - test it !
    # IMHO:
    # when the code is invalidated, it triggers somehow save_config() function
    # even if the user_config.yml file is present
    # then following lines are causing problems - in case the input$*params are not already created (because they are
    # created dynamically, when a user clicks on a checkbox) at the UI side - the value is NULL
    # I propose temporary fix - these values will be read either from input$[some_method_param] or 'user_config.yml' 
    # or from factory settings
      
    # SANN
    ifelse(length(input$SANN_params)!=0,
           config$yaml_content$running_params$maxit_SANN <- as.integer(input$SANN_params),
           ifelse(length(current_user_config$running_params$maxit_SANN)!=0 || config$yaml_content$optim_method$SANN==TRUE,
                  config$yaml_content$running_params$maxit_SANN <- current_user_config$running_params$maxit_SANN, # as.integer(input$SANN_params)
                  config$yaml_content$running_params$maxit_SANN <- list() # as.integer(input$SANN_params)
                  )
           
    )
    
    # nloptr
    ifelse(length(input$nloptr_params)!=0,
           config$yaml_content$running_params$maxit_nloptr <- as.integer(input$nloptr_params),
           ifelse(length(current_user_config$running_params$maxit_nloptr)!=0 || config$yaml_content$optim_method$nloptr==TRUE,
                  config$yaml_content$running_params$maxit_nloptr <- current_user_config$running_params$maxit_nloptr, 
                  config$yaml_content$running_params$maxit_nloptr <- list()
                  )
           )
    
    # NelderMead
    ifelse(length(input$gensa_params)!=0,
           config$yaml_content$running_params$max_iter_gensa <- as.integer(input$gensa_params),
           ifelse(length(current_user_config$running_params$max_iter_gensa)!=0 || config$yaml_content$optim_method$max_iter_gensa==TRUE,
                  config$yaml_content$running_params$max_iter_gensa <- current_user_config$running_params$max_iter_gensa, 
                  config$yaml_content$running_params$max_iter_gensa <- list()
           )
    )
    
    # rgenround
    ifelse(length(input$rgenoud_params)!=0,
           config$yaml_content$running_params$max_iter_rgenoud <- as.integer(input$rgenoud_params),
           ifelse(length(current_user_config$running_params$max_iter_rgenoud)!=0 || config$yaml_content$optim_method$max_iter_rgenoud==TRUE,
                  config$yaml_content$running_params$max_iter_rgenoud <- current_user_config$running_params$max_iter_rgenoud, 
                  config$yaml_content$running_params$max_iter_rgenoud <- list()
           )
    )
    
    # Nelder-Mead
    ifelse(length(input$NM_params)!=0,
           config$yaml_content$running_params$maxit_NM <- as.integer(input$NM_params),
           ifelse(length(current_user_config$running_params$maxit_NM)!=0 || config$yaml_content$optim_method$maxit_NM==TRUE,
                  config$yaml_content$running_params$maxit_NM <- current_user_config$running_params$maxit_NM, 
                  config$yaml_content$running_params$maxit_NM <- list()
                  )
           )
    
    # BFGS
    ifelse(length(input$BFGS_params)!=0,
           config$yaml_content$running_params$maxit_BFGS <- as.integer(input$BFGS_params),
           ifelse(length(current_user_config$running_params$maxit_BFGS)!=0 || config$yaml_content$optim_method$maxit_BFGS==TRUE,
                  config$yaml_content$running_params$maxit_BFGS <- current_user_config$running_params$maxit_BFGS, 
                  config$yaml_content$running_params$maxit_BFGS <- list()
                  )
           )
    
    # rel_tol
    ifelse(length(input$rel_tol_params)!=0,
           config$yaml_content$running_params$optim_rel_tol <- as.double(input$rel_tol_params),
           ifelse(length(current_user_config$running_params$optim_rel_tol)!=0 || config$yaml_content$optim_method$optim_rel_tol==TRUE,
                  config$yaml_content$running_params$optim_rel_tol <- current_user_config$running_params$optim_rel_tol, 
                  config$yaml_content$running_params$optim_rel_tol <- list()
                  )
           )
    
    # Write modified YAML content to file
    write_yaml(config$yaml_content, "user_config.yml")

  } # save_config() ENDS

  # Upload file button
  my_current_dir <- getwd()

  if (Sys.info()[[1]]=="Linux"){
    volumes <- getVolumes()
    shinyFileChoose(input, "file_data", roots = volumes,defaultPath = my_current_dir, filetypes = c("*", "txt"))
    }else{
      volumes<-c(wd=".")
      shinyFileChoose(input, "file_data", roots = volumes,filetypes = c("*", "txt"))
    }

  uploadedFileName <- reactiveVal(NULL)

  # Observe if file_data was provided
  observeEvent(input$file_data, {
    fileInfo <- parseFilePaths(roots = volumes, input$file_data)
    
   
    if (!is.null(fileInfo$datapath)) {
      fileName <- fileInfo$datapath
      uploadedFileName(fileName)
    }
  })

  # print out chosen file - TODO make it reactive
  output$chosen_file <- renderText({
    paste("No file selected")
  })

  # Create example data
  time <- c(2,8,15,25,38,45,60)
  dissolved <- c(3,12,35,48,55,58,66)
  
  # Plot example graph
  output$plot <- renderPlotly({plot_ly(x=time) %>%
      add_trace(y=dissolved,type="scatter",mode="lines+markers",name="dissolution") %>%
      layout(title="Example dissolution profile",xaxis=list(title="Time",linewidth=2),yaxis=list(title="Q [%]",linewidth=2))})

  # Create a graph - based on the loaded data - and make it reactive - each time file load - plot is generated
  plot_data <- reactive({

    # get file path
    file_info <- parseFilePaths(roots=volumes,input$file_data)
    filepath <- file_info$datapath
    filename <- file_info$name
    
    output$chosen_file <- renderText({
      filename
    })

    # put data into uploaded_data object via read.table
    uploaded_data <- read.table(filepath, header = cfg$data$headers, sep = "\t", colClasses="numeric")

    # plot
    plot_ly(x=~uploaded_data[,1]) %>%
      add_trace(y=~uploaded_data[,2],type="scatter",mode="lines+markers",name="dissolution") %>%
      layout(title="Dissolution profile",xaxis=list(title="Time",linewidth=2),yaxis=list(title="Q [%]",linewidth=2))

  })

  observeEvent(uploadedFileName(), {
    if (length(uploadedFileName())==0) {
      # Create example graph
      output$plot <- renderPlotly({
          plot_ly(x=time) %>%
          add_trace(y=dissolved,type="scatter",mode="lines+markers",name="dissolution") %>%
          layout(title="Example dissolution profile",xaxis=list(title="Time",linewidth=2),yaxis=list(title="Q [%]",linewidth=2))})
    } else {
    output$plot <- renderPlotly({
      plot_data()
    })
    }
  })


  # Display optim methods
  output$optim_methods <- renderUI({

    my_method_lst <- list()

    if ("SANN" %in% input$selected_optimmethods) {
      tmp_lst <- list(textInput(inputId = "SANN_params",
                                label = "Maximum number of iterations in SANN method:",
                                value = ifelse(length(cfg$running_params$maxit_SANN)==0, 10000, cfg$running_params$maxit_SANN)
      )
      )
      my_method_lst <- append(my_method_lst, tmp_lst)
    }

    if ("nloptr" %in% input$selected_optimmethods) {
      tmp_lst <- list(textInput(inputId = "nloptr_params",
                                label = "Maximum number of iterations in nloptr method:",
                                value = ifelse(length(cfg$running_params$maxit_nloptr)==0, 10000, cfg$running_params$maxit_nloptr)))
      my_method_lst <- append(my_method_lst, tmp_lst)
    }

    if ("NelderMead" %in% input$selected_optimmethods) {
      tmp_lst <- list(textInput(inputId = "NM_params",
                                label = "Maximum number of iterations in NelderMead method:",
                                value = ifelse(length(cfg$running_params$maxit_NM)==0, (as.integer(500000)), cfg$running_params$maxit_NM)))
      my_method_lst <- append(my_method_lst, tmp_lst)
    }

    if ("genSA" %in% input$selected_optimmethods) {
      tmp_lst <- list(textInput(inputId = "gensa_params",
                                label = "Maximum number of iterations in genSA method:",
                                value = ifelse(length(cfg$running_params$max_iter_gensa)==0, 1500, cfg$running_params$max_iter_gensa)))
      my_method_lst <- append(my_method_lst, tmp_lst)
    }

    if ("rgenoud" %in% input$selected_optimmethods) {
      tmp_lst <- list(textInput(inputId = "rgenoud_params",
                                label = "Maximum number of iterations in rgenoud method:",
                                value = ifelse(length(cfg$running_params$max_iter_rgenoud)==0, 500, cfg$running_params$max_iter_rgenoud)))
      my_method_lst <- append(my_method_lst, tmp_lst)
    }

    if(!is.null(my_method_lst)){
      do.call("tagList", my_method_lst)
    }

  })

  # Restore config settings
  reset_app <- function(working_conf){
    updateSelectInput(session, "selected_models", selected = working_conf$models[working_conf$selected_modelsx$selected_models])
    updateSelectInput(session, "selected_optimmethods", selected = names(working_conf$optim_method)[2])
    updateTextInput(session, "BFGS_params", value = working_conf$running_params$maxit_BFGS)
    updateTextInput(session, "rel_tol_params", value = working_conf$running_params$optim_rel_tol)

    updateSelectInput(session, "selected_optitrace", selected = (if (TRUE %in% working_conf$running_params$opti_trace) "Yes" else "No"))

    updateTextInput(session, "SANN_params", value = working_conf$running_params$maxit_SANN)
    updateTextInput(session, "nloptr_params", value = working_conf$running_params$maxit_nloptr)
    updateTextInput(session, "NM_params", value = working_conf$running_params$maxit_NM)
    updateTextInput(session, "gensa_params", value = working_conf$running_params$max_iter_gensa)
    updateTextInput(session, "rgenoud_params", value = working_conf$running_params$max_iter_rgenoud)

    updateSelectInput(session, "selected_headers", selected = (if (TRUE %in% working_conf$data$headers) "Yes" else "No"))
  }

  observeEvent(input$restore, {
    factory_cfg <- yaml::yaml.load_file('factory_config.yml')
    reset_app(factory_cfg)
  })

  # Tick all models button
  observeEvent(input$selection_options, {
    if (input$selection_options == "Select all models") {
      updateSelectInput(session, "selected_models", selected = cfg$models)}
    if (input$selection_options == "Deselect all models") {
      updateSelectInput(session, "selected_models", selected = 0)}
    if (input$selection_options == "Restore standard settings") {
      factory_cfg <- yaml::yaml.load_file('factory_config.yml')
      reset_app(factory_cfg)}
  })


  # Save yml file
  config <- reactiveValues(yaml_content = cfg) # list(yaml_content = cfg)
  save_button_to_listen <- reactive({
    input$save
  })

  # observe if any save button was pressed
  observeEvent(save_button_to_listen(), {
    # all logic was moved to save_config() - need to test if everything works here
    # this function writes configuration to 'user_config.yml'
    save_config(session = session)

  })

  observeEvent(input$save, {
    # Output config status message
    showModal(modalDialog(
      title = "Status","Settings saved",
      easyClose = TRUE,
      footer = tagList(
        modalButton("OK")
      )
    ))
  })
  
  #Observe if "Save" button is clicked
  saved <- reactiveVal(FALSE)
  
  observeEvent(input$save, {
    saved(TRUE)
  })
  
  #Observe if "Choose file" button is clicked
  chosen <- reactiveVal(FALSE)
  
  observeEvent(input$file_data, {
    chosen(TRUE)
  })

  # Function to disable/enable input elements.
  # input_list: List of inputs, reactiveValuesToList(input)
  # enable_inputs: Enable or disable inputs?
  # Only buttons: Toggle all inputs, or only buttons?
  toggle_inputs <- function(input_list, enable_inputs=TRUE, only_buttons=FALSE)  {
    # Subset if only_buttons is TRUE.
    if(only_buttons){
      buttons <- which(sapply(input_list,function(x) {any(grepl('Button',attr(x,"class")))}))
      input_list = input_list[buttons]
    }
    
    # Toggle elements
    for(x in names(input_list))
      if(enable_inputs){
        shinyjs::enable(x)} else {
          shinyjs::disable(x) }
  }
  
  #Perform computation
  observeEvent(input$run, {
    
    # First save config to user_config.yml
    # this can be contracted into more global function
    # and used - let's try -> first tests are promising - it seems this works well
    save_config(session = session)

    # Disable buttons/inputs
    input_list <- reactiveValuesToList(input)
    toggle_inputs(input_list, FALSE, FALSE)

    if (!chosen()) {
      showModal(modalDialog(
        title = "Warning",
        "Please upload file data before running!",
        easyClose = TRUE
      ))
    } else {
      
      # ############################# #
      # Main - start computation core #
      # ############################# #
      
      #Show message before computing
      showModal(modalDialog(
        title = "Status", "Wait, model computation in progress",
        easyClose = TRUE,
        footer = tagList(
          modalButton("OK")
        )
      ))
      
      # withProgress() - simple progress bar;
      # In order to enable you need to pass parameter progress=TRUE
      # then in computational core - an incProgress() function from
      # shiny updates the progress bar
      withProgress(message = 'Running... (this may take a while)', # withProgress function START
                   value = 0, {
                     
       #Running main function script
       RKinetDS_comp_core(progress=TRUE)
                     
      }) # withProgress() function END

      output$display_report <- renderText({
        file_report <- file.path("overall_report.txt")
        paste(readLines(file_report), collapse = "\n")
      })
        
      #Show RMSE error ranking
      RMSE_error <- file.path("RMSE_error_ranking.txt")
      
      output$table_error_RMSE <- renderDT({
      read.table(RMSE_error, header=TRUE, fill = TRUE)
        })
      
      #Show R2 error ranking
      R2_error <- file.path("R2_error_ranking.txt")
  
      output$table_error_R2 <- renderDataTable ({
        read.table(R2_error, header=TRUE, fill = TRUE)
      })
      
      #Show R2 adjusted error ranking
      R2_error_ad <- file.path("R2_adjusted_error_ranking.txt")
      
      output$table_error_R2_ad <- renderDataTable ({
        read.table(R2_error_ad, header=TRUE, fill = TRUE)
      })
      
      #Show AIC error ranking
      AIC_error <- file.path("AIC_error_ranking.txt")
      
      output$table_error_AIC <- renderDataTable ({
        read.table(AIC_error, header=TRUE, fill = TRUE)
      })
        
        # Show message after computing
        showModal(modalDialog(
          title = "Status", "Model computation finished",
          easyClose = TRUE,
          footer = tagList(
            modalButton("OK")
          )
        ))

        # Create tabs logic via renderUI
        # First take which models were selected - from config
        cfg <- yaml::yaml.load_file('user_config.yml')
        list_of_selected_model <- cfg$selected_modelsx$selected_models
        #
        # Here we have logic for creating tabs based on results obtained
        # in the main process from RKinetDS_comp_core()
        #
        output$tabs_with_models_results <- renderUI({
            
            results_tabs_list <- lapply(1:length(list_of_selected_model), function(i) {
              # Get models name 
              optimized_eq_path <- paste0("Results_eq_", list_of_selected_model[i], "/optimizedEquation_",
                                         list_of_selected_model[i],".txt", sep="")
              optimized_equation_file <- readLines(optimized_eq_path)
              models_name_line <- optimized_equation_file[grepl("Model name: ", optimized_equation_file)]
              models_name <- str_split(models_name_line, "Model name: ")[[1]][2]
              
              # plot, path and text output indices
              plot_model_index <- paste("plot_model_", i, sep="")
              txt_out_model_index <- paste("info_model_", i, sep="")
              path_info_index <- paste("folderPath_", i, sep="")
              # 
              # Tabs logic and construction
              tabPanel(
                "",
                title = models_name,
                inputId = "model_tabs",
                br(),
                plotlyOutput(plot_model_index),
                br(),
                p(strong("Model results:")),
                verbatimTextOutput(txt_out_model_index),
                p(strong("Results file path:")),
                textOutput(path_info_index),
                br(),
                br()
              )
            
              }) # lapply function ENDS HERE
            
            # Next do.call() need to be called
            do.call(tabsetPanel, results_tabs_list)
          })
          
          # Call renderPlot for each one. Plots are only actually generated when they
          # are visible on the web page.
          for (i in 1:length(list_of_selected_model)) {
            
            # Need local so that each item gets its own number. Without it, the value
            # of i in the renderPlot() will be the same across all instances, because
            # of when the expression is evaluated.
            local({
              my_i <- i
              # Some parts need to be repeated here also
              # in order to properly create plots and txt
              # Get path for txt files 
              data_path <- paste0("Results_eq_", list_of_selected_model[my_i], "/results_", list_of_selected_model[my_i],".txt", sep="")
              optimized_eq_path <- paste0("Results_eq_", list_of_selected_model[my_i], "/optimizedEquation_", list_of_selected_model[my_i],".txt", sep="")
              
              # plot, path and text output indices
              plot_model_index <- paste0("plot_model_", my_i, sep="")
              txt_out_model_index <- paste0("info_model_", my_i, sep="")
              path_info_index <- paste0("folderPath_", my_i, sep="")
              
              # Read data from files
              models_data <- read.table(data_path, header = TRUE, colClasses="numeric")
              optimized_equation_file <- readLines(optimized_eq_path)
              models_name_line <- optimized_equation_file[grepl("Model name: ", optimized_equation_file)]
              models_name <- str_split(models_name_line, "Model name: ")[[1]][2]
              params_after_optim <- optimized_equation_file[grepl("Parameters for equation after optimization: ", optimized_equation_file)]
              eq_algebraic_form <- optimized_equation_file[grepl("Algebraic form: ", optimized_equation_file)]
              eq_rmse <- optimized_equation_file[grepl("RMSE:  ", optimized_equation_file)]
              eq_rmse_value <- str_split(eq_rmse, "RMSE:  ")[[1]][2]
              
              # renderPlotly logic - put all into output$[plot_model_index]
              output[[plot_model_index]] <- renderPlotly({
                plot_ly(x = ~models_data[, 1]) %>%
                  add_trace(y = ~models_data[, 2], type = "scatter", mode = "markers", name = "OBS") %>%
                  add_trace(y = ~models_data[, 3], type = "scatter", mode = "lines+markers", name = "PRED") %>%
                  layout(title = paste("Model dissolution profile: ", models_name, sep=""),
                         xaxis = list(title = "Time", linewidth = 2), yaxis = list(title = "Q [%]", linewidth = 2))
              })

              output[[txt_out_model_index]] <- renderText({
                if (file.exists(optimized_eq_path)) {
                  paste(optimized_equation_file, collapse = "\n")
                } else {
                  "Results not found."
                }
              })

              #File path
              output[[path_info_index]] <- renderText({
              normalizePath(path.expand(optimized_eq_path))
              })
            
              })
            
          }

        plot_the_best_data <- reactive({
        # Read the best data from RMSE_error_ranking.txt
        error_table_txt <- read.csv("RMSE_error_ranking.txt", header = TRUE, sep="\t")
        the_best_model_index <- error_table_txt$Number[1]
        the_best_model_filename <- paste0("Results_eq_", the_best_model_index,
                                          "/results_", the_best_model_index,".txt", sep="")
        the_best_optimized_eq_filename <- paste0("Results_eq_", the_best_model_index, "/optimizedEquation_", the_best_model_index,".txt", sep="")
        the_best_optimized_eq_file <- readLines(the_best_optimized_eq_filename)
        the_best_models_name_line <- the_best_optimized_eq_file[grepl("Model name: ", the_best_optimized_eq_file)]
        the_best_models_name <- str_split(the_best_models_name_line, "Model name: ")[[1]][2]
        the_best_eq_algebraic_form <- the_best_optimized_eq_file[grepl("Algebraic form: ", the_best_optimized_eq_file)]
        the_best_eq_algebraic_form <- str_split(the_best_eq_algebraic_form, "Algebraic form: ")[[1]][2]
        the_best_eq_rmse <- the_best_optimized_eq_file[grepl("RMSE:  ", the_best_optimized_eq_file)]
        the_best_eq_rmse_value <- str_split(the_best_eq_rmse, "RMSE:  ")[[1]][2]
        the_best_data <- read.csv(the_best_model_filename, header = TRUE, sep="\t")

        # plot
        plot_ly(x=~the_best_data[,1]) %>%
          add_trace(y = ~the_best_data[, 2], type = "scatter", mode = "markers", name = "OBS") %>%
          add_trace(y = ~the_best_data[, 3], type = "scatter", mode = "lines+markers", name = "PRED") %>%
          layout(title = paste0("Model dissolution profile: ", the_best_models_name, sep=""),
                                annotations = list(xref="paper",
                                                   yref="paper",
                                                   x=1,
                                                   y=c(0.2,0.15,0.1),
                                                   text=c(paste("Model: ", the_best_models_name, sep=""),
                                                          paste("Equation: ", the_best_eq_algebraic_form, sep=""),
                                                          paste("RMSE: ", the_best_eq_rmse_value, sep="")
                                                        ),
                                    textposition = "left center",
                                    textfont = list(family= "Times", size= 18, color= "DarkOrange"),
                                    showarrow=FALSE
                                    ),
          xaxis = list(title = "Time", linewidth = 2), yaxis = list(title = "Q [%]", linewidth = 2))
        
        #  add_trace(y=~the_best_data[,2],type="scatter",mode="lines+markers",name="dissolution") %>%
        # layout(title="Dissolution profile",xaxis=list(title="Time",linewidth=2),yaxis=list(title="Q [%]",linewidth=2))
        })

        output$plot <- renderPlotly({
        plot_the_best_data()
        })

    } # # Main - start computation core  ELSE condition END

    # Enable inputs once again
    toggle_inputs(input_list,TRUE,FALSE)

  }) # observeEvent(input$run) FUNCTION END
 

  #Display report
  output$display_report <- renderText({
    "Results not found"})

  # Display used settings
   output$used_models_rep <- renderText({
    # Read the config, and make it a consistent with results section
    used_config <- user_config_reader() # user_config_reader() is a function defined at the beginning of the server
    paste(used_config$models[unlist(used_config$selected_modelsx$selected_models)], collapse = "\n")
    })

  output$used_optim_rep <- renderText({
    # Read the config, and make it a consistent with results section
    used_config <- user_config_reader() # user_config_reader() is a function defined at the begining of the server
    values_optim <- unlist(unname(used_config$optim_method))
    if (any(values_optim)) {
      paste(names(used_config$optim_method), ": ", used_config$optim_method, sep = "", collapse = "\n")
    }
    })

  output$used_param_rep <- renderText({
    # Read the config, and make it a consistent with results section
    used_config <- user_config_reader() # user_config_reader() is a function defined at the begining of the server
    paste("Maximum number of iterations in BFGS method: ", used_config$running_params$maxit_BFGS, sep="", "\n",
          "Maximum number of iterations in SANN method: ", used_config$running_params$maxit_SANN, "\n",
          "Maximum number of iterations in nloptr method: ", used_config$running_params$maxit_nloptr, "\n",
          "Maximum number of iterations in genSA method: ", used_config$running_params$max_iter_gensa, "\n",
          "Maximum number of iterations in rgenoud method: ", used_config$running_params$max_iter_rgenoud, "\n",
          "Maximum number of iterations in NelderMead method: ", used_config$running_params$maxit_NM, "\n",
          "Value of stop criterion for optimizing: ", used_config$running_params$optim_rel_tol, "\n",
          "Tracing of optimization function evaluations: ", used_config$running_params$opti_trace)
  })
  output$used_data_rep <- renderText({
    # Read the config, and make it a consistent with results section
    used_config <- user_config_reader() # user_config_reader() is a function defined at the begining of the server
    paste("Uploaded data included headers: ", used_config$data$headers, sep="")
  })
  
  #Display manual HTML
  output$manual <- renderUI({
    get_manual_page()
  })

  get_manual_page <- function(){
      return(includeHTML("Documentation_RKinetDS.html"))
  }

  #Display license HTML
  output$license <- renderUI({
    get_license_page()
  })
  
  get_license_page <- function() {
    return(includeHTML("license.html"))
  }
  
} # shiny server function END

shiny::shinyApp(ui, server,options=list(launch.browser=TRUE))
