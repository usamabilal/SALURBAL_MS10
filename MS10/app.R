library(shiny)
library(tidyverse)
library(DT)
library(leaflet)
library(htmltools)
library(plotly)
library(sf)
library(scales)
load("MS10_exposure_data.RData")
load("MS9_outcome_data.rdata")
load("data_LE.rdata")
load("l1s.rdata")
load("l1ad_shp.rdata")
shp_le<-inner_join(shp, dta_le)

dta<-mortality_cause %>% select(SALID1, all_of(cause_titles_coll), total) %>% 
    gather(type, deaths, -SALID1, -total) %>% 
    mutate(prop=deaths/total) %>% 
    left_join(l1s) %>% 
    # shortening tucuman
    mutate(city_link=ifelse(grepl("Tucuman", city_link), "Tucuman", city_link))

shp_mortality<-inner_join(shp, dta %>% mutate(prop=prop*100))

# testing stuff
# input<-list()
# input$age<-input$age5<-input$age4<-input$age3<-input$age2<-"birth"
# input$sex<-input$sex5<-input$sex4<-input$sex3<-input$sex2<-"men"
# input$country4<-input$country2<-input$country3<-input$country5<-"All"
ui <- fluidPage(
    navbarPage("Life Expectancy and Mortality Profiles in Latin American Cities",
               tabPanel("Introduction",
                        h1("Introduction"),
                        hr(),
                        p("This app provides extra information and data on the manuscript entitled", 
                          em("Life expectancy and mortality in 363 cities of Latin America: the SALURBAL project")),
                        p("The following tabs are available:"),
                        tags$ol(tags$li(strong("Life Expectancy Distribution: "), "shows a boxplot with the distribution of life expectancy at different ages and sexes for all cities"), 
                                tags$li(strong("Life Expectancy Uncertainty: "), "shows a linerange plot with the distribution of life expectancies at different ages and sexes for all cities [or restricted to a specific country] and the corresponding 95% Credible Intervals"),
                                tags$li(strong("Life Expectancy Table: "), "shows a data table with the life expectancies at different ages and sexes for all cities [or restricted to a specific country] and the corresponding 95% Credible Intervals"),
                                tags$li(strong("Life Expectancy Map: "), "maps, for a specific country, the life expectancies at different ages and sexes"),
                                tags$li(strong("PM Distribution: "), "shows an interactive stacked bar plot with the proportion of deaths due to each cause"), 
                                tags$li(strong("PM Table: "), "shows a data table with proportionate mortality by cause"),
                                tags$li(strong("PM Map: "), "maps, for a specific country, the proportionate mortality by cause")),
                        p("Code for the app and analysis is available here: ",
                          a(href="https://github.com/usamabilal/SALURBAL_MS10",
                            "https://github.com/usamabilal/SALURBAL_MS10",  target="_blank")),
                        p("The manuscript and derived data from this analysis is licensed under a ", a(href="https://creativecommons.org/licenses/by/4.0/", "CC-BY 4.0 license.", target="_blank"))),
               tabPanel("Life Expectancy Distribution",
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
               tabPanel("Life Expectancy Uncertainty",
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
               tabPanel("Life Expectancy Table",
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
                                            selected = "All"),
                                downloadButton("downloadData_le", "Download")
                            ),
                            mainPanel(
                                DT::dataTableOutput("table2")
                            )
                        )
               ),
               tabPanel("Life Expectancy Map",
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
               ),
               tabPanel("PM Distribution",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("cause1", "Sort by:", 
                                             choices=c("CMNN", "Cancer", "CVD/NCDs","Unintentional","Violence"), 
                                             selected="Violence"),
                                selectInput("country1", label="Country", 
                                            choices=c("All", "AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "All")
                            ),
                            mainPanel(
                                plotlyOutput("plot1_pm", height = "200%"),
                                textOutput("text1_pm")
                            )
                        )
               ),
               tabPanel("PM Table",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("country2pm", label="Country", 
                                            choices=c("All", "AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "All"),
                                downloadButton("downloadData_pm", "Download")
                            ),
                            mainPanel(
                                DT::dataTableOutput("table2_pm")
                            )
                        )
               ),
               tabPanel("PM Map",
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("cause3", "Show %:", 
                                             choices=c("CMNN", "Cancer", "CVD/NCDs","Unintentional","Violence"), 
                                             selected="Violence"),
                                selectInput("country3", label="Country", 
                                            choices=c("All", "AR", "BR", "CL", "CO", "CR/PA/SV",
                                                      "MX", "PE"), 
                                            selected = "AR"),
                                p("Note: Be patient...map takes some time to load. You can interact with it (click, zoom, pan, etc.)")
                            ),
                            mainPanel(
                                leafletOutput("plot3_pm")
                            )
                        )
               ))
)


server <- function(input, output) {
    output$plot1 <- renderPlot({
        tage<-as.numeric(ifelse(input$age=='birth', 0, input$age))
        tsex<-ifelse(input$sex=="men", "M", "F")
        temp<-dta_le %>% filter(age==tage&sex==tsex)
        
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
        temp<-dta_le %>% filter(age==tage&sex==tsex)
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
        
        plot<-ggplot(temp, aes(x=id, y=le, group=city_link)) +
                     geom_linerange(aes(ymin=lci, ymax=uci), size=.5, alpha=1) +
                     geom_point(size=1, aes(text=label)) +
                     #geom_jitter(aes(color=as.factor(iso2)), width=0.1, height=0, alpha=0.5, size=2) +
                     guides(size=F, color=guide_legend(override.aes = list(alpha = 1, size=1)))+
                     coord_flip()+
                     xlab("") + ylab("Years") +
                     scale_color_discrete(name="")+
                     scale_y_continuous(limits = c(NA, NA), sec.axis = dup_axis(name = ""),
                                        breaks=seq(range[1],range[2], by=2.5))+
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
                           plot.title=element_text(face="bold", size=25))
                 
        ggplotly(plot, tooltip = c("text")) %>% 
            config(displayModeBar=T, displaylogo=F,
                     modeBarButtonsToRemove = list(
                         'sendDataToCloud',
                         'toImage',
                         'autoScale2d',
                         'hoverClosestCartesian',
                         'hoverCompareCartesian',
                         'select2d','lasso2d'
                     ), collaborate = F) %>% 
            layout(dragmode="zoom", 
                   xaxis= list(fixedrange=T))
    })
    output$downloadData_le <- downloadHandler(
        filename = function() {
            paste("salurbal_le_data.csv")
        },
        content = function(file) {
            tage<-as.numeric(ifelse(input$age2=='birth', 0, input$age2))
            tsex<-ifelse(input$sex2=="men", "M", "F")
            temp<-dta_le %>% filter(age==tage&sex==tsex)
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
            temp<-temp %>% arrange(desc(le)) %>% 
                ungroup() %>% 
                select(country_name, city_link, age, sex,le, lci, uci) %>% 
                mutate(le=round(le, digits=1),
                       lci=round(lci, digits=1),
                       uci=round(uci, digits=1),
                       sex=ifelse(sex=="F", "Women", "Men")) %>% 
                rename(country=country_name, city=city_link)
            write.csv(temp, file, row.names = FALSE)
        }
    )
    output$table2 <- DT::renderDataTable({
        tage<-as.numeric(ifelse(input$age2=='birth', 0, input$age2))
        tsex<-ifelse(input$sex2=="men", "M", "F")
        temp<-dta_le %>% filter(age==tage&sex==tsex)
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
    output$plot4 <- renderLeaflet({
        tage<-as.numeric(ifelse(input$age4=='birth', 0, input$age4))
        tsex<-ifelse(input$sex4=="men", "M", "F")
        shp2<-shp_le %>% filter(age==tage&sex==tsex)
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
    output$text1_pm <- renderText({ "Note: You can interact with the plot, zoom, pan, etc." })
    output$plot1_pm <- renderPlotly({
        type_var<-case_when(
            input$cause1=="CMNN" ~ "cmnn",
            input$cause1=="Cancer" ~ "cancer",
            input$cause1=="CVD/NCDs" ~ "ncd",
            input$cause1=="Unintentional" ~ "accident",
            input$cause1=="Violence" ~ "violent",
        )
        if (input$country1=="CR/PA/SV"){
            temp<-dta %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country1=="All"){
            temp<-dta
        } else {
            temp<-dta %>% filter(iso2==input$country1)    
        }
        # get ordering of countries based on selected variable
        country_order<-temp %>% 
            group_by(iso2, type) %>% 
            summarise(deaths=sum(deaths),
                      total=sum(total)) %>% 
            filter(type==type_var) %>% 
            mutate(prop=deaths/total) %>% 
            arrange(prop) %>% select(iso2) %>% 
            ungroup() %>% 
            mutate(country_id=row_number()) 
        temp<-full_join(temp, country_order)
        # get ordering of cities overall
        city_order<-temp %>% 
            filter(type==type_var) %>% 
            arrange(country_id, prop) %>% select(SALID1) %>% 
            mutate(new_id=row_number())
        temp<-temp %>% full_join(city_order) %>% 
            mutate(id=case_when(
                type=="cmnn"~"1cmnn",
                type=="cancer"~"2cancer",
                type=="ncd"~"3ncd",
                type=="accident"~"4accident",
                type=="violent"~"5violent"),
                id2=case_when(
                    type=="cmnn"~"Communicable/Maternal/Neonatal/Nutritional",
                    type=="cancer"~"Cancer",
                    type=="ncd"~"CVD and other NCD",
                    type=="accident"~"Unintentional Injuries",
                    type=="violent"~"Violent Injuries"),
                id2=factor(id2, levels=c("Communicable/Maternal/Neonatal/Nutritional", 
                                         "Cancer", "CVD and other NCD",
                                         "Unintentional Injuries", "Violent Injuries")),
                # label for x axis
                city_label=paste0(city_link, ", ", iso2),
                # label for bar
                city_label2=paste0(city_label, "; % ", id2, " ", 
                                   round(prop*100, digits=1), "%")) %>% 
            arrange(country_id, new_id)
        # separating countries
        limits<-temp %>% group_by(country_id) %>% 
            summarise(limit=max(new_id)+0.5) %>% 
            arrange(limit) %>% slice(-n()) %>% pull(limit)
        # if length(limits) is 0, just make it NA
        if (length(limits)==0) limits<-NA
        labelsize<-ifelse(input$country1%in%c("MX", "BR"), 8,
                          ifelse(input$country1=="All", 6, 12))
        plot<-ggplot(temp, aes(x=as.factor(new_id), y=prop, 
                               fill=as.factor(id2),text=city_label2))+
            geom_bar(width = 1, stat = "identity", color=NA, size=0)+
            geom_vline(xintercept = limits, lty=2, color="black")+
            scale_fill_discrete(labels=c("Communicable/Maternal/Neonatal/Nutritional", 
                                         "Cancer", "CVD and other NCD",
                                         "Unintentional Injuries", "Violent Injuries"),
                                name="Cause")+
            scale_alpha_manual(values=c(0,1))+
            coord_cartesian(clip = "off", ylim=c(0, 1))+
            scale_x_discrete(labels=unique(temp$city_label))+
            scale_y_continuous(expand=c(0,0),
                               labels=percent,
                               limits=c(0, 1.03),
                               breaks=seq(0, 1, by=.2),
                               sec.axis = sec_axis(~.*1, name = "Proportionate Mortality",
                                                   labels=percent,
                                                   breaks=seq(0, 1, by=.2)))+
            labs(x="", y="Proportionate Mortality", title="")+
            guides(fill=guide_legend(nrow=2))+
            theme_classic() +
            theme(legend.position = "bottom",
                  legend.text=element_text(size=14, color="black"),
                  legend.title=element_text(size=20, color="black", face="bold"),
                  legend.background=element_blank(),
                  plot.title=element_text(face="bold", size=20),
                  axis.title.x=element_text(face="bold", color="black"),
                  axis.text.x=element_text(color="black", size=labelsize,angle=90, 
                                           hjust=1, vjust=.5),
                  axis.title.y=element_text(face="bold", color="black", size=14),
                  axis.text.y=element_text(color="black", size=14))
        plot_plotly<-ggplotly(plot,tooltip = c("text")) %>% 
            config(displayModeBar=T, displaylogo=F,
                   modeBarButtonsToRemove = list(
                       'sendDataToCloud',
                       'toImage',
                       'autoScale2d',
                       'hoverClosestCartesian',
                       'hoverCompareCartesian',
                       'select2d','lasso2d')) %>% 
            layout(dragmode="zoom", showlegend=T,
                   yaxis= list(fixedrange=T),
                   legend = list(orientation = "h",
                                 x=.5, y=1,
                                 yanchor="bottom",
                                 xanchor="left",
                                 itemclick=F,
                                 font=list(size=10)))
        # cleaning legend, from: https://stackoverflow.com/questions/49133395/strange-formatting-of-legend-in-ggplotly-in-r
        for (i in 1:length(plot_plotly$x$data)){
            if (!is.null(plot_plotly$x$data[[i]]$name)){
                plot_plotly$x$data[[i]]$name =  gsub("\\(","",
                                                     str_split(plot_plotly$x$data[[i]]$name,",")[[1]][1])
            }
        }
        plot_plotly
    })
    output$downloadData_pm <- downloadHandler(
        filename = function() {
            paste("salurbal_pm_data.csv")
        },
        content = function(file) {
            if (input$country2pm=="CR/PA/SV"){
                temp<-dta %>% filter(iso2%in%c("CR", "PA", "SV"))
            } else if (input$country2pm=="All"){
                temp<-dta
            } else {
                temp<-dta %>% filter(iso2==input$country2pm)    
            }
            temp<-temp %>% 
                arrange(desc(total)) %>% 
                ungroup() %>% 
                select(country_name, city_link, type, prop) %>% 
                mutate(prop=prop*100,
                       prop=round(prop, digits=1)) %>% 
                spread(type, prop) %>% 
                select(country_name, city_link, cmnn, cancer, ncd, accident, violent) %>% 
                rename(country=country_name, city=city_link)
            write.csv(temp, file, row.names = FALSE)
        }
    )
    output$table2_pm <- DT::renderDataTable({
        if (input$country2pm=="CR/PA/SV"){
            temp<-dta %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country2pm=="All"){
            temp<-dta
        } else {
            temp<-dta %>% filter(iso2==input$country2pm)    
        }
        temp<-temp %>% 
            arrange(desc(total)) %>% 
            ungroup() %>% 
            select(country_name, city_link, type, prop) %>% 
            mutate(prop=prop*100) %>% 
            spread(type, prop) %>% 
            select(country_name, city_link, cmnn, cancer, ncd, accident, violent)
        datatable(temp, colnames = c("Country", "City", 
                                     "% CMNN", "% Cancer", "% CVD/NCDs", "% Unintentional Injuries", "% Violent Injuries"),
                  options=list(pageLength=363)) %>% 
            formatRound(columns=c('cmnn', 'cancer', 'ncd', 'accident', 'violent'), digits=1)
    })
    output$plot3_pm <- renderLeaflet({
        type_var<-case_when(
            input$cause3=="CMNN" ~ "cmnn",
            input$cause3=="Cancer" ~ "cancer",
            input$cause3=="CVD/NCDs" ~ "ncd",
            input$cause3=="Unintentional" ~ "accident",
            input$cause3=="Violence" ~ "violent",
        )
        shp2<-shp_mortality %>% filter(type==type_var)
        if (input$country3=="CR/PA/SV"){
            shp2<-shp2 %>% filter(iso2%in%c("CR", "PA", "SV"))
        } else if (input$country3=="All"){
            shp2<-shp2
        } else {
            shp2<-shp2 %>% filter(iso2==input$country3)    
        }
        gradient_pal =  colorNumeric(palette="RdYlGn",  # colors we want to use
                                     domain=shp2$prop, reverse = T) #supply range of possible values 
        
        #geo=as.numeric(colMeans(st_coordinates(st_centroid(shp2))))
        shp2 %>% 
            mutate(tag = str_c("<b>",city_link,"</b>","<br/>",   #Create label
                               "<i>LE: ",paste0(round(prop, digits=1), "%"),"</i>") %>%
                       map(htmltools::HTML)) %>%
            leaflet(options= leafletOptions(dragging = T, #remove panning
                                            minZoom = 1,  #set zoom limits
                                            maxZoom = 10)) %>% 
            addTiles() %>% 
            # setView(lng=geo[1],
            #         lat=geo[2],
            #         zoom = 4) %>% 
            addPolygons(weight = 3, 
                        fillOpacity = 1,
                        color= ~gradient_pal(prop),
                        label = ~tag,
                        highlight = highlightOptions(weight=3,
                                                     color = 'red',
                                                     bringToFront = T)) %>% 
            addLegend(pal = gradient_pal, values = ~prop,
                      title =paste0('%',input$cause3),opacity = 1,group='prop') 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
