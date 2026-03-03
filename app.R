
# ===== app.R =====
# Contain the app for dashbord with shinny package

# ===== Libraries =====
library(shiny)           # Dashboard package
library(shinydashboard)  # Dashboard package
library(fresh)           # Personnalise the theme (colors)
library(patchwork)       # Arrange plots
library(DT)              # Tables

# ===== Import data =====
#load("~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.DataImport") # import gene


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
    selectInput("marker","Marker", # choix utilisateur # Label
                choices = c("GH","18S","ITS")),
    selectInput("level","Level", # choix utilisateur # Label
                choices = c("Family","Genus","Species")),
    selectInput("name","Name", # choix utilisateur # Label
                choices = c("All")), # à update dans le serveur
    selectInput("var_bota","Variable bota", # choix utilisateur # Label
                choices = c("Count","NCI", "Gaussian","OTU","Div_alpha")),
    selectInput("var_genet","Variable genet", # choix utilisateur # Label
                choices = c("Weight","Weight_density", "OTU","Div_alpha"))
    ),
  
  dashboardBody(
    # ===== ITEMS =====
    tabItems(
      # ===== ONGLET 1 : VUE GLOBALE =====
      tabItem(tabName="global",
              
              fluidRow(
                box(plotOutput("renderPlot_compare"), width=10), # 100% = 12
                valueBoxOutput("Name", width=2),
                valueBoxOutput("R2", width=2),
                DTOutput("R2_tab")
                #box(plotOutput("renderBarplot_R2_taxon"),width=12, height = "5000px")
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
      bota       = bota,
      level      = input$level,
      name       = input$name,
      var_bota = input$var_bota
    )
    
  })
  
  # Genet grid 
  genet_grid <- reactive({
    
    genet_grid_read(
      PROJECT_PATH,
      input$marker, # Marker
      input$level,  # Level
      input$name    # Name
      )
    
  })
  
  # Inventaire botanique
  bota_f <- reactive({
    df <- bota %>% select_tax(input$level, input$name)
    return(df)
  })
  
  # Dynamic compare
  df_compare <- reactive({
    df <- R2_taxon(Marker = input$marker,
                   Level = input$level,
                   var_bota = input$var_bota,
                   var_genet   = input$var_genet)
    return(df)
  })
  
  # ===== VUE GLOBALE =====
  # Plot compare
  output$renderPlot_compare <- renderPlot({
    p1 <- plot_grid_stem(bota_grid(), bota_f(),input$var_genet) +
      labs(title = "Bota density")
    p2 <- plot_grid_stem(genet_grid(), bota_f(),input$var_genet)+
      labs(title = "Genet density")
    p3 <- plot_compare(
      df_compare_join( # data prep
        genet_grid(),
        bota_grid(),
        var_bota = input$var_bota,
        var_genet = input$var_genet
      ),
      bota_f() # filtered bota
    )
    # Plot patchwork
    (p1 | p2 | p3) 
  })
  # R-squared
  output$R2 <- renderValueBox({
    r2 <- R_squared(
      df_compare_join(
        genet_grid(),
        bota_grid(),
        var_bota = input$var_bota,
        var_genet = input$var_genet
      ), 
      var_obs = input$var_genet,
      var_pred = input$var_bota
    ) %>% round(digits = 3)
    valueBox(paste0("R2 = ",r2),"", color = "blue")
  })
  # Name box
  output$Name <- renderValueBox({
    valueBox(input$name,paste0("Level : ",input$level," / Marker : ", input$marker), color = "blue")
  })
  # Barplot R_suaqred
  output$renderBarplot_R2_taxon <- renderPlot({ 
      barplot_R2_taxon( df_compare()) 
      },
    height = function() {nrow(df_compare()) * 40}
  )
  # Tableau de R2
  output$R2_tab <- renderDT({
    datatable(
      df_compare(),        # dataframe
      options = list(pageLength = 10), # 10 lignes par page
      rownames = TRUE
    )
  })
  
}

# ===== DAHBORD LAUNCH =====
shinyApp(ui, server)



