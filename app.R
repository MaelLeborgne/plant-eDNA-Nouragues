
# ===== app.R =====
# Contain the app for dashbord with shinny package

# ===== Libraries =====
library(shiny)           # Dashboard package
library(shinydashboard)  # Dashboard package
library(fresh)           # Personnalise the theme (colors)
library(patchwork)       # Arrange plots
library(DT)              # Tables

# ===== Import data =====
load(
  "~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.DataImport.RData",
  envir = .GlobalEnv
  ) # import gene, bota...


# ===== Source data and plots =====
source("data_prep.R")
source("plots.R")

# ===== User Interface =====
ui <- dashboardPage(
  # Title
  dashboardHeader(title="ADNe Nourague"),
  # Tabs
  dashboardSidebar(
    sidebarMenu(
      menuItem("Vue globale", tabName="global", icon=icon("chart-line"))
    ),
    # Filters
    selectInput("marker","Marker", # Variable # Label
                choices = c("GH","18S","ITS")),
    selectInput("level","Level",
                choices = c("Family","Genus","Species")),
    selectInput("name","Name",
                choices = c("All")), # à update dans le serveur
    selectInput("var_bota","Variable bota",
                choices = c("Count","NCI", "Gaussian","OTU","Shannon")),
    selectInput("norm","Normalisation", 
                choices = c(TRUE,FALSE)),
    selectInput("var_genet","Variable genet",
                choices = c("Weight","Weight_density", "OTU","Shannon"))
    ),
  
  dashboardBody(
    # ===== ITEMS =====
    tabItems(
      # ===== ONGLET 1 : VUE GLOBALE =====
      tabItem(tabName="global",
              
              fluidRow(
                box(plotOutput("renderPlot_compare"), width=9), # 100% = 12
                valueBoxOutput("Name", width=3),
                valueBoxOutput("R2", width=3),
                valueBoxOutput("Corr_pears", width=3),
                DTOutput("R2_tab"),
                box(plotOutput("renderPlot_x_y"), width=6)
              )
      )
    )
  )
)


# ===== SERVER =====
server <- function(input, output, session){
  
  # ===== FILTER ===== 
  # # Update 'name' filter
  observeEvent(
    list(input$marker, input$level),
    {
      
      req(input$marker, input$level)
      
      noms <- tryCatch(
        give_level_names(input$marker, input$level),
        error = function(e) character(0)
      )
      
      # Sécurisation du format
      noms <- as.character(noms)
      noms <- noms[!is.na(noms)]
      noms <- unique(noms)
      
      # Filtrage pour ne garder que ceux présents dans bota[[input$level]]
      if(!is.null(bota[[input$level]])) {
        noms <- noms[noms %in% bota[[input$level]]]
      } else {
        noms <- character(0)
      }
      
      updateSelectInput(
        session,
        "name",
        choices = c("All", noms)
      )
    }
  )
  # Bota grid 
  bota_grid <- reactive({
    
    bota_grid_compute(
      bota_df    = bota,
      level      = input$level,
      name       = input$name,
      var_bota   = input$var_bota,
      user_norm  = input$norm,
      radius = 10
    )
    
  })
  
  # Genet grid 
  genet_grid <- reactive({
    
    genet_grid_read( PROJECT_PATH,
                     input$marker, # Marker
                     input$level,  # Level
                     input$name,    # Name
                     from_csv = TRUE,
                     total = FALSE)
    
  })
  
  # Inventaire botanique
  bota_f <- reactive({
    
    bota %>% select_tax(input$level, input$name)
    
  })
  
  # Dynamic compare
  df_compare <- reactive({
    
    df_compare_join(
      grid_pred = genet_grid(),
      grid_obs  = bota_grid(),
      var_pred  = input$var_genet,
      var_obs   = input$var_bota)
    
  })
  
  # Dynamic compare
  df_R2 <- reactive({
    
    R2_taxon(Marker     = input$marker,
             Level      = input$level,
             var_bota   = input$var_bota,
             var_genet  = input$var_genet,
             user_norm  = input$norm)
    
  })
  
  # ===== VUE GLOBALE =====
  # Plot compare
  output$renderPlot_compare <- renderPlot({
    stems = FALSE
    if(input$name != "All"){
      stems = TRUE
    }
    patch_compare(var_bota = input$var_bota,
                  var_genet = input$var_genet,
                  b_grid = bota_grid(),
                  g_grid = genet_grid(),
                  b_f    = bota_f(),
                  stems  = stems)
  })
  
  # R-squared
  output$R2 <- renderValueBox({
    r2 <- R_squared(
      df_compare_join( # data prep
        grid_pred = genet_grid(),
        grid_obs  = bota_grid(),
        var_pred  = input$var_genet,
        var_obs   = input$var_bota
      )
    ) %>% round(digits = 3)
    valueBox(paste0("R2 = ",r2),"", color = "blue")
  })
  # Corrélation pears
  output$Corr_pears <- renderValueBox({
    corr <- corr_b_g(
      df_compare_join( # data prep
        grid_pred = genet_grid(),
        grid_obs  = bota_grid(),
        var_pred  = input$var_genet,
        var_obs   = input$var_bota
      )
    ) %>% round(digits = 2)
    valueBox(paste0("Corr = ",corr),"", color = "blue")
  })
  # Name box
  output$Name <- renderValueBox({
    valueBox(input$name,paste0("Level : ",input$level," / Marker : ", input$marker), color = "blue")
  })
  # Barplot R_suaqred
  output$renderBarplot_R2_taxon <- renderPlot({ 
      barplot_R2_taxon( df_R2()) 
      },
    height = function() {nrow(df_R2()) * 40}
  )
  
  # Tableau de R2
  output$R2_tab <- renderDT({
    datatable(
      df_R2(),        # dataframe
      options = list(pageLength = 10), # 10 lignes par page
      rownames = TRUE
    )
  })
  
  # Plot 2D Bota / Genet
  output$renderPlot_x_y <- renderPlot({
    plot_x_y(
      df_compare_join( # data prep
        grid_pred = genet_grid(),
        grid_obs  = bota_grid(),
        var_pred  = input$var_genet,
        var_obs   = input$var_bota
      ), 
      x_param = paste0("Bota :",input$var_bota), 
      y_param = paste0("Genet :",input$var_genet)
      )
  })
  
  # Tableau de Richness
  # output$tab_rich <- renderDT({
  #   datatable(
  #     tab_OTU_stat(
  #       Marker = input$marker,
  #       Name = 
  #     ),        # dataframe
  #     options = list(pageLength = 10), # 10 lignes par page
  #     rownames = TRUE
  #   )
  # })
  
}

# ===== DAHBORD LAUNCH =====
shinyApp(ui, server)



