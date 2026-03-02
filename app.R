
# ===== app.R =====
# Contain the app for dashbord with shinny package

# ===== Libraries =====
library(shiny)           # Dashboard package
library(shinydashboard)  # Dashboard package
library(fresh)           # Personnalise the theme (colors)
library(patchwork) # Arrange plots

# ===== Import data =====
#load("~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.ModelEnv") # import gene


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
                choices = c("All")) # à update dans le serveur
    ),
  
  dashboardBody(
    # ===== ITEMS =====
    tabItems(
      # ===== ONGLET 1 : VUE GLOBALE =====
      tabItem(tabName="global",
              
              fluidRow(
                box(plotOutput("renderPlot_compare"), width=12) # 100% = 12
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
    # Compute grid
    df <- select_tax(
      bota,            # inventaire botanique
      input$level,     # Level
      input$name) %>%  # Name
      bota_grid_NCI()
    # Return
    return(df)
  })
  # Genet grid 
  genet_grid <- reactive({
    # Read grid
    df <- genet_grid_read(
      PROJECT_PATH,
      input$marker, # Marker
      input$level,  # Level
      input$name    # Name
      )
    # Return
    return(df)
  })
  # Inventaire botanique
  bota_f <- reactive({
    df <- bota %>% select_tax(input$level, input$name)
    return(df)
  })
  
  # ===== VUE GLOBALE =====
  # Plot compare
  output$renderPlot_compare <- renderPlot({
    p1 <- plot_grid_stem(bota_grid(), bota_f(),"density") +
      labs(title = "Bota density")
    p2 <- plot_grid_stem(genet_grid(), bota_f(),"density")+
      labs(title = "Genet density")
    p3 <- plot_compare(
      df_compare_join( # data prep
        genet_grid(), # reactive data
        bota_grid()   # reactive data
      ),
      bota_f() # filtered bota
    )
    # Plot patchwork
    (p1 | p2 | p3) 
  })
  
}

# ===== DAHBORD LAUNCH =====
shinyApp(ui, server)



