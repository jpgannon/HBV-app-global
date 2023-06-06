library(shiny)
library(shinythemes)
library(tidyverse)
library(lubridate)
library(plotly)
library(leaflet)
library(sf)
library(rtop)

theme_set(theme_classic())
source("HBV.R")
source("hbvnse.R")
source("NPE-KlingGupta.R")
source("hbvNPEKG.R")

end <- ymd("2014-10-01")
start <-  ymd("2009-10-01")

template <- tribble(
    ~"DATE",~"Q_mmd",~"prcp_mmd",~"T_mean_C",~"swe_mm",
    2009-10-02,0,0,0,0
)

NSEpars <- data.frame(matrix(ncol = 14, nrow = 100)) %>%
  mutate_all(~replace(., is.na(.), 0))
colnames(NSEpars) <- c("FC", "beta", "LP", "SFCF", "TT", 
                       "CFMAX", "k0", "k1", "k2", "UZL", 
                       "PERC", "MAXBAS", "NSE", "NPEKG")
write_csv(NSEpars, "NSEpars.csv")

#Read shapefile
watersheds <- st_read("gages700.shp")

gages <- read_csv("gage_information_caravan700.csv")

rando <- function(){
    FC    <- runif(1, min = 40   , max = 400)  #Max soil moisture storage, field capacity
    beta  <- runif(1, min = 1    , max = 6)    #Shape coefficient governing fate of water input to soil moisture storage
    LP    <- runif(1, min = 0.3  , max = 1)    #Threshold for reduction of evap
    SFCF  <- runif(1, min = 0.4  , max = 1.2)  #Snowfall correction factor
    TT    <- runif(1, min = -1.5 , max = 1.2)  #Threshold temperature
    CFMAX <- runif(1, min = 1    , max = 8)    #Degree-day factor
    k0    <- runif(1, min = 0.05 , max = 0.5)  #Recession constant (upper storage, near surface)
    k1    <- runif(1, min = 0.01 , max = 0.3)  #Recession constant (upper storage)
    k2    <- runif(1, min = 0.001, max = 0.15) #Recession constant (lower storage)
    UZL   <- runif(1, min = 0    , max = 70)   #Threshold for shallow storage
    PERC  <- runif(1, min = 0    , max = 4)    #Percolation, max flow from upper to lower storage
    MAXBAS<- rep(1)
    
    c(FC, beta, LP, SFCF, TT, CFMAX, k0, k1, k2, UZL, PERC, MAXBAS)
}

startpars <- rando()

# Define UI for application 
ui <- fluidPage(
    theme = shinytheme("lumen"),
    # Application title
    titlePanel("Run the HBV Model for one of 700 CARAVAN Watersheds"),

    # Sidebar with a sliders for parameters
    sidebarLayout(
        sidebarPanel(width = 4,
            selectInput("gage", "Choose a watershed here or using the map", 
                        choices = c(unique(gages$GAGE_NAME), "User Data"), 
                        selected = gages$GAGE_NAME[5]),
            dateRangeInput("dates", "Run the model for these dates:",
                           start = start, end = end, min = start, max = end),
            hr(style = "border-top: 1px solid #000000;"),
            h5(HTML("<b>Parameters for water movement through storages</b>")),
            sliderInput("FC"     ,"Soil field capacity", min = 40   , max = 400 , value = startpars[1]),  #Max soil moisture storage, field capacity
            sliderInput("beta"   ,"Shape coeff for water delivery to soil", min = 1    , max = 6   , value = startpars[2]),    #Shape coefficient governing fate of water input to soil moisture storage
            sliderInput("UZL"    ,"Threshold for shallow storage", min = 0    , max = 70  , value = startpars[10]),   #Threshold for shallow storage
            sliderInput("k0"     ,"Recession constant: near surface storage", min = 0.05 , max = 0.5 , value = startpars[7]),  #Recession constant (upper storage, near surface)
            sliderInput("k1"     ,"Recession constant: upper storage", min = 0.01 , max = 0.3 , value = startpars[8]),  #Recession constant (upper storage)
            sliderInput("k2"     ,"Recession constant: lower storage", min = 0.001, max = 0.15, value = startpars[9]), #Recession constant (lower storage)
            sliderInput("PERC"   ,"Percolation: max flow from upper-lower", min = 0    , max = 4   , value = startpars[11]),    #Percolation, max flow from upper to lower storage
            sliderInput("LP"     ,"Threshold for Evap reduction", min = .3   , max = 1   , value = startpars[3]),    #Threshold for reduction of evap
            h5(HTML("<b>Parameters for snow routine</b>")),
            sliderInput("SFCF"   ,"Snowfall correction factor", min = 0.4  , max = 1.2 , value = startpars[4]),  #Snowfall correction factor
            sliderInput("TT"     ,"Threshold temperature", min = -1.5 , max = 1.2 , value = startpars[5]),  #Threshold temperature
            sliderInput("CFMAX"  ,"Degree-day factor", min = 1    , max = 8   , value = startpars[6]),    #Degree-day factor
            hr(style = "border-top: 1px solid #000000;"),
            #sliderInput("MAXBAS" ,"Base of Triangular routing function", min = 1    , max = 3   , value = startpars[12]),    #base of the triangular routing function, days
            numericInput("lat"   ,"Latitude: (for PET calculation)", value = gages$LAT[gages$GAGE_NAME == gages$GAGE_NAME[1]]),
            #numericInput("elev"  ,"Elevation: ", value = gages$Elevation_m[gages$GAGE_NAME == gages$GAGE_NAME[1]]),
            hr(style = "border-top: 1px solid #000000;"),
            h5(HTML("<b>Parameterization tools</b>")),
            actionButton("genrando", "Generate Random Parameter Set"),
            radioButtons("objfxn", label = "Objective Function to optimize", choices = c("NSE", "NPKGE"), selected = "NSE"),
            sliderInput("snowWT" , "Weight of snow in objective function calculation", min = 0, max = 1, value = 0),
            numericInput("runs"  , "# of Runs for Monte Carlo", min = 100, max = 10000, value = 100),
            actionButton("runMC" , "Run Monte Carlo"),
            h5(HTML("Click below to run a SCEUA optimization. This will take a few minutes.")),
            actionButton("runSCEUA", "Run SCEUA Optimization")
        ),

        # Show output plots
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("HBV Output", plotlyOutput("QPlot", height = 1000)),
            
            tabPanel("Select on Map",
                     leafletOutput("Map")
                     ),
            tabPanel("Monte Carlo Plots",
                     h4("Click button to create parameter vs. nse plots after running a Monte Carlo. 
                        If plots show a single point, the run hasn't been completed."),
                     actionButton("Refresh", "Create Monte Carlo Plots"),
                     plotOutput("NSEplot", height = "800px")),
           tabPanel("Upload Data",
                    h4("Click the button below to download a template input file.
                        Add your DAILY data to this sheet and then upload as a CSV
                        using the upload button below."),
                    downloadButton("downloadTemplate","Download Input Data Template"),
                    radioButtons("datetype", "Select the date format of your input file",
                                 choices = c("mdy", "ymd", "dmy"), selected = "ymd"),
                    fileInput("UserData", "Upload filled-in template file by clicking Browse below",
                              accept = c("text/csv",
                                         "text/comma-separated-values/text/plain",
                                         ".csv")),
                    h5("Select User Data in the Choose your watershed dropdown to the left.
                       If everything went well your data will be shown below."),
                    actionButton("clear", "Clear User Input"),
                    tableOutput("input_contents")
           ),
           tabPanel("Model Diagram",
                    plotOutput("diagram"))
           
            )
           
        )
    )
)

# Define server logic 
server <- function(input, output, session) {
    
    output$downloadTemplate <- downloadHandler(
        filename = "InputTemplate.csv",
        content = function(file) {
            write.csv(template, file, row.names = FALSE)
        }
    )
    
    inFile <- NULL
    
    #if no user input file, read default
    inputs_all <- reactive({
    
    #get user input if it exists
    inFile <- input$UserData
        
    if(input$gage != "User Data"){ #is.null(inFile)){
        inputs_all <- read_csv(paste0("Caravan700/", 
                                      gages$GAGE_ID[gages$GAGE_NAME == input$gage],
                                      ".csv")) %>%
                    rename(Precip = prcp_mmd, Qobs = Q_mmd, 
                           Temp = T_mean_C, Snow = swe_mm) 
    }else{
        inputs_all <- read_csv(inFile$datapath) %>%
            rename(Precip = prcp_mmd, Qobs = Q_mmd, 
                   Temp = T_mean_C, Snow = swe_mm) %>%
            mutate(DATE = case_when(
                input$datetype == "mdy" ~ mdy(DATE),
                input$datetype == "ymd" ~ ymd(DATE),
                input$datetype == "dmy" ~ dmy(DATE)))
    }
    return(inputs_all)
    })
    
    observeEvent(input$clear,{
       # values$upload_state <- 'reset'
    })
    
    observeEvent(input$gage,{
    #ELEV <- gages$Elevation_m[gages$GAGE_NAME == input$gage]
    LAT <- gages$LAT[gages$GAGE_NAME == input$gage]
    
    #updateNumericInput(inputId = "elev", 
      #                 value = ELEV)
    
    updateNumericInput(inputId = "lat", 
                       value = LAT)
    
    updateDateRangeInput(inputId = "dates", 
                         start = inputs_all()$DATE[1], 
                         end = inputs_all()$DATE[length(inputs_all()$DATE)],
                         min = inputs_all()$DATE[1], 
                         max = inputs_all()$DATE[length(inputs_all()$DATE)])
    })
    
    input_data <- reactive({
        input_data <- filter(inputs_all(), 
                        DATE >= input$dates[1] & DATE <= input$dates[2])
        input_data
    })
    
    #generate output table for input file tab to show user what data is read in
    output$input_contents <- renderTable({
         if(input$gage == "User Data") mutate(input_data(), DATE = as.character(DATE))
        else mutate(template, DATE = as.character(DATE))
    })
    
    #Calculate PET
    
    PET <- reactive({
        
        latrad <- (input$lat/360) * 2 * pi #convert to radians
        
        PET1 <- input_data() %>% select(DATE) %>%
        mutate(DOY = yday(DATE)) %>% #DOY for dates
        mutate(tempvar = (2 * pi / 365) * DOY) %>%
        mutate(delta_h = 0.4093 * sin(tempvar - 1.405)) %>% # declination of the sun above the celestial equator in radians on day JulDay of the year#
        mutate(daylen = (2 * acos(-tan(delta_h) * tan(latrad)) / 0.2618)) %>% #day length, h
        mutate(PET = 29.8 * daylen * 0.611 * exp(17.3 * input_data()$Temp / (input_data()$Temp + 237.3)) / (input_data()$Temp + 273.2))  #PET Hamon method
    
        PET1$PET
    })
    
    HBVout <- reactive({
        routing <- 0
        pars <- c(input$FC, input$beta, input$LP, input$SFCF, 
                      input$TT, input$CFMAX, input$k0, input$k1, 
                      input$k2, input$UZL, input$PERC, input$MAXBAS)
        
        modeloutput <- HBV(pars, input_data()$Precip, input_data()$Temp, PET(), routing)
        
       #add observations for plotting
        bind_cols(modeloutput, input_data()$Qobs)
    })
    
    output$QPlot <- renderPlotly({
        output <- bind_cols(HBVout(), 
                            DATE = input_data()$DATE,
                            Qobs = input_data()$Qobs)
        snowout <- bind_cols(output, 
                             Snow = input_data()$Snow)
        
        #make snow record to calc NSE
        justsnow <- snowout %>% drop_na(Snow)
        
        EvalStart <- floor(length(output$q) * 0.4)
        EvalEnd <- length(output$q)
        
        outputNSE <- output[EvalStart:EvalEnd,]
        justsnow <- justsnow %>%
            filter(DATE >= min(outputNSE$DATE), DATE <= max(outputNSE$DATE))
        
        #Calculate NSE and add to parameter set
        NSE <- 1 - ((sum((outputNSE$Qobs - outputNSE$q) ^ 2, na.rm = TRUE)) / 
        sum((outputNSE$Qobs - mean(outputNSE$Qobs, na.rm = TRUE)) ^ 2, na.rm = TRUE))
        
        snowNSE <- 1 - ((sum((justsnow$Snow - justsnow$SWE) ^ 2, na.rm = TRUE)) / 
        sum((justsnow$Snow - mean(justsnow$Snow, na.rm = TRUE)) ^ 2, na.rm = TRUE))
        
        if(is.na(snowNSE) == TRUE) snowNSE <- 0
        if(snowNSE == -Inf) snowNSE <- 0
        
        NSE <- (NSE * (1 - input$snowWT)) + (snowNSE * input$snowWT)
        
        NPEKG <- RNP(sim = outputNSE$q,
                     obs = outputNSE$Qobs)
        
        snowNPEKG <- RNP(sim = justsnow$SWE,
                         obs = justsnow$Snow)
        
        if(snowNPEKG == -Inf) snowNPEKG <- 0
        if(is.na(snowNPEKG) == TRUE) snowNSE <- 0
        
        NPEKG <- (NPEKG * (1 - input$snowWT)) + (snowNPEKG * input$snowWT)
        
        vline <- function(x = 0, color = "red") {
            list(
                type = "line", 
                y0 = 0, 
                y1 = 1, 
                yref = "paper",
                x0 = x, 
                x1 = x, 
                line = list(color = color, dash = "dot")
            )
        }
        
        table1 <- plot_ly(
          type = 'table',
          header = list(
            values = c('<b>Parameter</b>', '<b>Value</b>'),
            line = list(color = '#506784'),
            fill = list(color = '#119DFF'),
            align = c('left','center'),
            font = list(color = 'white', size = 12)
          ),
          cells = list(
            values = rbind(
              c('NSE', 'NPE-KG', 'Observed mean Q', 'Modeled mean Q'),
              c(round(NSE, 2), round(NPEKG, 2), 
                paste(round(mean(output$Qobs, na.rm = TRUE),2), "mm/day"), 
                paste(round(mean(output$q, na.rm = TRUE),2), "mm/day"))),
            line = list(color = '#506784'),
            fill = list(color = c('white', 'white')),
            align = c('left', 'center'),
            font = list(color = c('#506784'), size = 12)
          ))
        
        plot1 <- output %>% plot_ly(x = ~DATE)
        plot1 <- plot1 %>% add_trace(y = ~q, name = 'Q Modeled', mode = 'lines')
        plot1 <- plot1 %>% add_trace(y = ~Qobs, name = 'Q Measured', mode = 'lines')
    
        plot1 <- plot1 %>% layout(yaxis = list(title = "Discharge(mm)")) %>%
                        layout(shapes = list(vline(output$DATE[EvalStart]))) 

        plot2 <- output %>% plot_ly(x = ~DATE, y = ~P, 
                                    type = 'bar', name = "Precip") %>%
            layout(yaxis = list(title = "Precip (mm)"))%>%
            layout(shapes = list(vline(output$DATE[EvalStart]))) 

        plot3 <- output %>% plot_ly(x = ~DATE) %>%
            add_trace(y = ~Storage, mode = 'lines', name = " Total Storage") %>%
            add_trace(y = ~S1, name = 'Upper Storage', mode = 'lines')%>%
            add_trace(y = ~S2, name = 'Lower Storage', mode = 'lines') %>%
            layout(yaxis = list(title = "Storage (mm)"))%>%
            layout(shapes = list(vline(output$DATE[EvalStart]))) 

        plot4 <- output %>% plot_ly(x = ~DATE) %>%
            add_trace(y = ~AET, mode = 'lines', name = "AET") %>%
            add_trace(y = ~PET, name = 'PET', mode = 'lines') %>%
            layout(yaxis = list(title = "ET (mm)"))%>%
            layout(shapes = list(vline(output$DATE[EvalStart]))) 

        plot5 <- snowout %>% plot_ly(x = ~DATE) %>%
            add_trace(y = ~SWE, mode = 'lines', name = "Snow (modeled)") 
        
        plot5 <- plot5 %>% add_trace(y = ~Snow, name = 'Snow (measured)', 
                                       mode = 'lines' )%>%
            layout(yaxis = list(title = "Snowpack (mm)")) %>%
            layout(shapes = list(vline(output$DATE[EvalStart]))) 
        
        plot6 <- output %>% plot_ly(x = ~DATE) %>%
          add_trace(y = ~Temp, mode = 'lines', name = "Temp (C) measured") %>%
          layout(yaxis = list(title = "Temp (deg C)"))
        
        subplot(table1, plot1, plot2, plot3, plot4, plot5, plot6,
                shareX = TRUE, titleY = TRUE, 
                nrows = 7, heights = c(.10, .2, .14, .14, .14, .14, .14)) #.25, .15, .15, .15, .15, .15
    })
    
    
      
      observeEvent(input$Refresh,
      { if(input$objfxn == "NSE") output$NSEplot <- renderPlot({
                  NSEparsdat <- read_csv("NSEpars.csv") %>% 
                  select(-MAXBAS) %>%
                  top_n(100, NSE) %>%
                  pivot_longer(cols = -c(NSE, NPEKG))
                NSEparsdat %>%
                  ggplot(aes(x = value, y = NSE))+
                  geom_point()+
                  facet_wrap(facets = vars(name), scales = "free_x", ncol = 3)+
                  ggtitle("Parameter values and NSE of top 100 runs")+
                  theme(text = element_text(size = 20))})
        if(input$objfxn == "NPKGE")output$NSEplot <- renderPlot({
                  NSEparsdat <- read_csv("NSEpars.csv") %>% 
                  select(-MAXBAS) %>%
                  top_n(100, NPEKG) %>%
                  pivot_longer(cols = -c(NSE, NPEKG))
                NSEparsdat %>%
                  ggplot(aes(x = value, y = NPEKG))+
                  geom_point()+
                  facet_wrap(facets = vars(name), scales = "free_x", ncol = 3)+
                  ggtitle("Parameter values and NPEKG of top 100 runs")+
                  theme(text = element_text(size = 20))
      })
    })
    
    #run monte carlo simulation when button is clicked
    HBVoutMC <- reactive({
        NSEmax <- 0
        bestpars <- 0
        routing <- 0
        dummy <- input$FC
       
        NSEpars <- data.frame(matrix(ncol = 14, nrow = input$runs))
        colnames(NSEpars) <- c("FC", "beta", "LP", "SFCF", "TT", 
                               "CFMAX", "k0", "k1", "k2", "UZL", 
                               "PERC", "MAXBAS", "NSE", "NPEKG")
        
        for (i in 1:input$runs){
            pars <- rando()
            
            modeloutput <- HBV(pars, input_data()$Precip, input_data()$Temp, PET(), routing)
        
            #add observations for plotting
            results <- bind_cols(modeloutput, 
                                 Qobs = input_data()$Qobs,
                                 Snow = input_data()$Snow)
            
            EvalStart <- floor(length(input_data()$Qobs) * 0.4)
            EvalEnd <- length(input_data()$Qobs)
            
            #trim the first 40% of the record so it isn't included in the NSE calculation
            results <- results[EvalStart:EvalEnd,]
            
            #make snow record to calc NSE
            
            justsnow <- select(results, Snow, SWE)
            
            #Calculate NSE and add to parameter set
            NSE <-  1 - ((sum((results$Qobs - results$q) ^ 2, na.rm = TRUE)) / 
                             sum((results$Qobs - mean(results$Qobs, na.rm = TRUE)) ^ 2, na.rm = TRUE))
            
            snowNSE <- 1 - ((sum((justsnow$Snow - justsnow$SWE) ^ 2, na.rm = TRUE)) / 
                                sum((justsnow$Snow - mean(justsnow$Snow, na.rm = TRUE)) ^ 2, na.rm = TRUE))
            
            NPEKG <- RNP(sim = results$q,
                         obs = results$Qobs)
            
            snowNPEKG <- RNP(sim = justsnow$SWE,
                             obs = justsnow$Snow)
            
            if(snowNSE == -Inf) snowNSE <- 0
            if(snowNPEKG == -Inf) snowNPEKG <- 0
            
            NSE <- (NSE * (1 - input$snowWT)) + (snowNSE * input$snowWT)
            NPEKG <- (NPEKG * (1 - input$snowWT)) + (snowNPEKG * input$snowWT)
            #add NSE and pars df output
            NSEpars[i,] <- c(pars, NSE, NPEKG)
         
        }
        #bestpars
        if(input$objfxn == "NSE"){
        NSEpars <- arrange(NSEpars, desc(NSE))
        write_csv(NSEpars, "NSEpars.csv")
        NSEpars
        }
        
        if(input$objfxn == "NPKGE"){
          NSEpars <- arrange(NSEpars, desc(NPEKG))
          write_csv(NSEpars, "NSEpars.csv")
          NSEpars
        }
    })
    
    HBVoutSCEUA <- reactive({
        dummy <- input$FC
        
        lowpar <- c(40, 1, 0.3, 0.4, -1.5, 1, 0.05, 0.01, 0.001, 0, 0, 1)   
        highpar <- c(400, 6, 1, 1.2, 1.2, 8, 0.5, 0.3, 0.15, 70, 4, 3)
        
        pars <- rando()
        
        dat <- bind_cols(Precip = input_data()$Precip,
                         Temp =  input_data()$Temp,
                         PET = PET())
        
        obs <- input_data()$Qobs
        
        if(input$objfxn == "NSE"){
        sceuaout <- sceua(hbvnse, pars = pars, 
                          lower = lowpar, upper = highpar, 
                          dat = dat, obs = obs, routing = 0)
        }
        if(input$objfxn == "NPKGE"){
          sceuaout <- sceua(hbvnpe, pars = pars, 
                            lower = lowpar, upper = highpar, 
                            dat = dat, obs = obs, routing = 0)
        }
        sceuaout$par

        })
    
    observeEvent(input$runMC, {
        
        updateNumericInput(inputId = "runs", value = input$runs)
        
        allpars <- HBVoutMC()
        mcpars <- as.numeric(allpars[1,])
        
        updateSliderInput(inputId = "FC", value =    mcpars[1])
        updateSliderInput(inputId = "beta" , value = mcpars[2])
        updateSliderInput(inputId = "LP"   , value = mcpars[3])
        updateSliderInput(inputId = "SFCF" , value = mcpars[4])
        updateSliderInput(inputId = "TT"   , value = mcpars[5])
        updateSliderInput(inputId = "CFMAX", value = mcpars[6])
        updateSliderInput(inputId = "k0"   , value = mcpars[7])
        updateSliderInput(inputId = "k1"   , value = mcpars[8])
        updateSliderInput(inputId = "k2"   , value = mcpars[9])
        updateSliderInput(inputId = "UZL"  , value = mcpars[10])
        updateSliderInput(inputId = "PERC" , value = mcpars[11])
        }, once = FALSE)
    
    observeEvent(input$runSCEUA, {

        scpars <- HBVoutSCEUA()
        
        updateSliderInput(inputId = "FC",    value = scpars[1])
        updateSliderInput(inputId = "beta" , value = scpars[2])
        updateSliderInput(inputId = "LP"   , value = scpars[3])
        updateSliderInput(inputId = "SFCF" , value = scpars[4])
        updateSliderInput(inputId = "TT"   , value = scpars[5])
        updateSliderInput(inputId = "CFMAX", value = scpars[6])
        updateSliderInput(inputId = "k0"   , value = scpars[7])
        updateSliderInput(inputId = "k1"   , value = scpars[8])
        updateSliderInput(inputId = "k2"   , value = scpars[9])
        updateSliderInput(inputId = "UZL"  , value = scpars[10])
        updateSliderInput(inputId = "PERC" , value = scpars[11])
    }, once = FALSE)
    
    
    #generate new parameter set if button is clicked 
    observeEvent(input$genrando, {
        
        newpars <- rando()
        
        updateSliderInput(inputId = "FC", value = newpars[1])
        updateSliderInput(inputId = "beta" , value = newpars[2])
        updateSliderInput(inputId = "LP"   , value = newpars[3])
        updateSliderInput(inputId = "SFCF" , value = newpars[4])
        updateSliderInput(inputId = "TT"   , value = newpars[5])
        updateSliderInput(inputId = "CFMAX", value = newpars[6])
        updateSliderInput(inputId = "k0"   , value = newpars[7])
        updateSliderInput(inputId = "k1"   , value = newpars[8])
        updateSliderInput(inputId = "k2"   , value = newpars[9])
        updateSliderInput(inputId = "UZL"  , value = newpars[10])
        updateSliderInput(inputId = "PERC" , value = newpars[11])
    })
    
    output$Map <- renderLeaflet({
        leaflet() %>% #addTiles() %>%
            addProviderTiles("OpenTopoMap",
                             options = providerTileOptions(noWrap = TRUE)) %>%
            addPolygons(data = watersheds) %>%
            addCircleMarkers(layerId = gages$GAGE_NAME, 
                             lat = gages$LAT, #input$lat_input, 
                             lng = gages$LONG,
                             radius = 4,
                             fillOpacity = 0.65) %>%
            addCircleMarkers(lat = gages$LAT[gages$GAGE_NAME == input$gage],
                             lng = gages$LONG[gages$GAGE_NAME == input$gage],
                             #radius = 4, 
                             color = "red") %>%
            setView(lat = gages$LAT[gages$GAGE_NAME == input$gage], 
                    lng = gages$LONG[gages$GAGE_NAME == input$gage],
                    zoom = 5)
    })
    #select sites by clicking on them
    observeEvent(input$Map_marker_click, {
        site <- input$Map_marker_click
        siteID <- site$id
        updateSelectInput(session, "gage", selected = siteID)
    })
    
    output$diagram <- renderImage({
      filename <- normalizePath(file.path('./www/HBV_diagram.png'))
      list(src = filename)
    }, deleteFile = FALSE)
}
# Run the application 
shinyApp(ui = ui, server = server)
