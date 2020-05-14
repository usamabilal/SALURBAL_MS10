library(shiny)
library(tidyverse)
library(DT)
library(leaflet)
library(htmltools)
library(plotly)
library(sf)
load("data.rdata")
load("l1ad_shp.rdata")
shp<-inner_join(shp, dta)
ui <- fluidPage(
    navbarPage("Life Expectancy in Latin American Cities",
               tabPanel("Introduction",
                        h1("Introduction"),
                        hr(),
                        p("This app provides extra information and data on the manuscript entitled", em("Variation and Predictors of Life Expectancy across 363 Cities in 9 Latin American Countries: the SALURBAL study")),
                        p("The following tabs are available:"),
                        tags$ol(tags$li(strong("Distribution: "), "shows a boxplot with the distribution of life expectancy at different ages and sexes for all cities"), 
                                tags$li(strong("Uncertainty: "), "shows a linerange plot with the distribution of life expectancies at different ages and sexes for all cities [or restricted to a specific country] and the corresponding 95% Credible Intervals"),
                                tags$li(strong("Table: "), "shows a data table with the life expectancies at different ages and sexes for all cities [or restricted to a specific country] and the corresponding 95% Credible Intervals"),
                                tags$li(strong("Map: "), "maps, for a specific country, the life expectancies at different ages and sexes")),
                        p("Code for the app and analysis is available here: "),
                        a(href="https://github.com/usamabilal/SALURBAL_MS10",
                          "https://github.com/usamabilal/SALURBAL_MS10",  target="_blank")),
               tabPanel("Distribution",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("age", "Life Expectancy at Age:", 
                                             choices=c("birth", "20", "40", "60"), 
                                             selected="birth"),
                                radioButtons("sex", "Sex:", 
                                             choices=c("men", "women"), 
                                             selected="women")
                            ),
                            mainPanel(
                                plotOutput("plot1")
                            )
                        )
               ),
               tabPanel("Uncertainty",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("age5", "Life Expectancy at Age:", 
                                             choices=c("birth", "20", "40", "60"), 
                                             selected="birth"),
                                radioButtons("sex5", "Sex:", 
                                             choices=c("men", "women"), 
                                             selected="women"),
                                selectInput("country5", label="Country", 
                                            choices=c("All", "AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "AR")
                            ),
                            mainPanel(
                                plotlyOutput("plot5", height = "200%"),
                                textOutput("text5")
                            )
                        )
               ),
               tabPanel("Table",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("age2", "Life Expectancy at Age:", 
                                             choices=c("birth", "20", "40", "60"), 
                                             selected="birth"),
                                radioButtons("sex2", "Sex:", 
                                             choices=c("men", "women"), 
                                             selected="women"),
                                selectInput("country2", label="Country", 
                                            choices=c("All", "AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "All")
                            ),
                            mainPanel(
                                DT::dataTableOutput("table2")
                            )
                        )
               ),
               tabPanel("Map",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("age4", "Life Expectancy at Age:", 
                                             choices=c("birth", "20", "40", "60"), 
                                             selected="birth"),
                                radioButtons("sex4", "Sex:", 
                                             choices=c("men", "women"), 
                                             selected="women"),
                                selectInput("country4", label="Country", 
                                            choices=c("All","AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "AR"),
                                p("Note: Be patient...map takes some time to load. You can interact with it (click, zoom, pan, etc.)")
                            ),
                            mainPanel(
                                leafletOutput("plot4")
                            )
                        )
               ))
)


server <- function(input, output) {
    output$plot1 <- renderPlot({
        tage<-as.numeric(ifelse(input$age=='birth', 0, input$age))
        tsex<-ifelse(input$sex=="men", "M", "F")
        temp<-dta %>% filter(age==tage&sex==tsex)
        
        title<-paste0("Life Expectancy at ",
                      ifelse(tage==0, "birth", paste0("age ", tage)),
                      " for ", ifelse(tsex=="M", "Men", "Women")," in 363 Latin American cities")
        ggplot(temp, aes(x=country_name, y=le, group=country_name)) +
            geom_boxplot(aes(group=as.factor(country_name)), fill=NA, outlier.color = NA, width=0.5)+
            geom_jitter(aes(fill=as.factor(country_name)), width=0.1, height=0, alpha=1, size=2, 
                        color="black", pch=21) +
            guides(color=F, fill=F, size=F)+
            labs(x="",
                 y="Years")+
            scale_y_continuous(sec.axis=dup_axis())+
            theme_bw() +
            theme(legend.position = "bottom",
                  legend.key.width = unit(50, "points"),
                  panel.grid.major.x = element_blank(),
                  axis.text.x=element_text(face="bold", size=14, angle=90, hjust=1, 
                                           color="black", vjust=.5),
                  axis.text.y=element_text(size=16, color="black"),
                  axis.title.y=element_text(face="bold", size=20),
                  plot.title=element_text(face="bold", size=25))
    })
    output$text5 <- renderText({ "Note: You can interact with the plot, zoom, pan, etc." })
    output$plot5 <- renderPlotly({
        tage<-as.numeric(ifelse(input$age5=='birth', 0, input$age5))
        tsex<-ifelse(input$sex5=="men", "M", "F")
        temp<-dta %>% filter(age==tage&sex==tsex)
        if (input$country5=="CR/PA/SV"){
            temp<-temp %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country5=="All"){
            temp<-temp
        } else {
            temp<-temp %>% filter(iso2==input$country5)    
        }
        labelsize<-ifelse(input$country5%in%c("MX", "BR"), 6,
                          ifelse(input$country5=="All", 4, 10))
        temp<-temp[order(temp$le),]
        temp$id<-1:nrow(temp)
        title<-paste0("Life Expectancy at ",
                      ifelse(tage==0, "birth", paste0("age ", tage)),
                      " for ", ifelse(tsex=="M", "Men", "Women"))
        temp$label<-paste0(temp$city_link, ". LE = ", 
                           round(temp$le, digits=1), " (", 
                           round(temp$lci,digits=1), "-", 
                           round(temp$uci, digits=1), ")")
        range<-c(min(temp$lci), max(temp$uci))
        range<-c(floor(range[1]/5)*5,
                 ceiling(range[2]/5)*5)
        
        ggplotly(ggplot(temp, aes(x=id, y=le, group=city_link)) +
                     geom_linerange(aes(ymin=lci, ymax=uci), size=.5, alpha=1) +
                     geom_point(size=1, aes(text=label)) +
                     #geom_jitter(aes(color=as.factor(iso2)), width=0.1, height=0, alpha=0.5, size=2) +
                     guides(size=F, color=guide_legend(override.aes = list(alpha = 1, size=1)))+
                     coord_flip()+
                     xlab("") + ylab("Years") +
                     scale_color_discrete(name="")+
                     scale_y_continuous(limits = c(NA, NA), sec.axis = dup_axis(name = ""),
                                        breaks=seq(range[1],range[2], by=5))+
                     scale_x_continuous(expand=c(0.01, 0.01),
                                        breaks=temp$id,
                                        labels=paste0(temp$city_link, ", ", temp$iso2))+
                     #ggtitle(title)+
                     theme_classic() +
                     theme(legend.position = "bottom",
                           legend.key.width = unit(50, "points"),
                           panel.grid.major.x = element_line(color=scales::alpha("gray", 1), linetype=2),
                           panel.grid.minor.x = element_line(color=scales::alpha("gray", 0.5), linetype=2),
                           axis.text.x=element_text(size=14, color="black"),
                           axis.text.y=element_text(size=labelsize, color="black"),
                           axis.ticks.y=element_blank(),
                           axis.title.y=element_text(face="bold", size=20),
                           plot.title=element_text(face="bold", size=25)),
                 tooltip = c("text")
        ) %>% config(displayModeBar=T, displayLogo=F,
                     modeBarButtonsToRemove = list(
                         'sendDataToCloud',
                         'toImage',
                         'autoScale2d',
                         'hoverClosestCartesian',
                         'hoverCompareCartesian',
                         'select2d','lasso2d'
                     ), collaborate = F) %>% 
            layout(dragmode="pan")
    })
    output$table2 <- DT::renderDataTable({
        tage<-as.numeric(ifelse(input$age2=='birth', 0, input$age2))
        tsex<-ifelse(input$sex2=="men", "M", "F")
        temp<-dta %>% filter(age==tage&sex==tsex)
        title<-paste0("Life Expectancy at ",
                      ifelse(tage==0, "birth", paste0("age ", tage)),
                      " for ", ifelse(tsex=="M", "Men", "Women")," in 363 Latin American cities")
        if (input$country2=="CR/PA/SV"){
            temp<-temp %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country2=="All"){
            temp<-temp
        } else {
            temp<-temp %>% filter(iso2==input$country2)    
        }
        temp<-temp %>% arrange(desc(le))
        temp<-temp %>% ungroup() %>% select(country_name, city_link, le, lci, uci)
        datatable(temp, colnames = c("Country", "City", "Estimate", 
                                     "Lower 95% CrI", "Upper 95% CrI"),
                  options=list(pageLength=363)) %>% 
            formatRound(columns=c('le', 'lci', 'uci'), digits=1)
    })
    output$plot3 <- renderLeaflet({
        tage<-as.numeric(ifelse(input$age3=='birth', 0, input$age3))
        tsex<-ifelse(input$sex3=="men", "M", "F")
        temp<-dta %>% filter(age==tage&sex==tsex)
        shp<-inner_join(shp, temp)
        gradient_pal =  colorNumeric(palette="RdYlGn",  # colors we want to use
                                     domain=shp$le, reverse = F) #supply range of possible values 
        geo=c(-66, -8)
        shp %>% 
            mutate(tag = str_c("<b>",city_link,"</b>","<br/>",   #Create label
                               "<i>LE: ",paste0(round(le, digits=1), " (", round(lci, digits=1), ";", round(uci, digits=1), ")"),"</i>") %>%
                       map(htmltools::HTML)) %>%
            leaflet(options= leafletOptions(dragging = T, #remove panning
                                            minZoom = 1,  #set zoom limits
                                            maxZoom = 10)) %>% 
            addTiles() %>% 
            setView(lng=geo[1],
                    lat=geo[2],
                    zoom = 2) %>% 
            addPolygons(weight = 3, # change boundary pixel to 1 (default is 5)
                        fillOpacity = 0.7,
                        color= ~gradient_pal(le),
                        label = ~tag,
                        highlight = highlightOptions(weight=3,
                                                     color = 'red',
                                                     bringToFront = T)) %>% 
            addLegend(pal = gradient_pal, values = ~le,
                      title ='Life Expectancy',opacity = 1,group='le') 
        
    })
    output$plot4 <- renderLeaflet({
        tage<-as.numeric(ifelse(input$age4=='birth', 0, input$age4))
        tsex<-ifelse(input$sex4=="men", "M", "F")
        shp2<-shp %>% filter(age==tage&sex==tsex)
        if (input$country4=="CR/PA/SV"){
            shp2<-shp2 %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country4=="All"){
            shp2<-shp2
        } else {
            shp2<-shp2 %>% filter(iso2==input$country4)    
        }
        gradient_pal =  colorNumeric(palette="RdYlGn",  # colors we want to use
                                     domain=shp2$le, reverse = F) #supply range of possible values 
        
        shp2 %>% 
            mutate(tag = str_c("<b>",city_link,"</b>","<br/>",   #Create label
                               "<i>LE: ",paste0(round(le, digits=1), " (", round(lci, digits=1), ";", round(uci, digits=1), ")"),"</i>") %>%
                       map(htmltools::HTML)) %>%
            leaflet(options= leafletOptions(dragging = T, #remove panning
                                            minZoom = 1,  #set zoom limits
                                            maxZoom = 10)) %>% 
            addTiles() %>% 
            addPolygons(weight = 3, # change boundary pixel to 1 (default is 5)
                        fillOpacity = 1,
                        color= ~gradient_pal(le),
                        label = ~tag,
                        highlight = highlightOptions(weight=3,
                                                     color = 'red',
                                                     bringToFront = T)) %>% 
            addLegend(pal = gradient_pal, values = ~le,
                      title ='Life Expectancy',opacity = 1,group='le') 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
