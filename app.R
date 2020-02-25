library(raster)
library(ggplot2)
library(deldir)
library(vegan)
library(plyr)
library(gridExtra)
library(reshape2)
library(FactoMineR)
library(factoextra)
#library(chngpt)
library(stringr)
library(gstat)
library(sp)
library(scales)
library(lme4)
library(nlme)
#library(glmmTMB)
#library(indicspecies)
#library(umap)
#library(UpSetR)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(shinyjs)
library(mdatools)

# shiny::runGitHub("jberry47/Shiny-SorgSys_Fields","jberry47")

f16soil <- read.table("data/5ff4c4c88319c6a96e14c5b061e3c91f-schachtman_2016_soil.tsv",header=T,stringsAsFactors=F,sep="\t")
f16biom <- read.table("data/5ff4c4c88319c6a96e14c5b061e3c91f-schachtman_2016_biom.tsv",header=T,stringsAsFactors=F,sep="\t")
f17soil <- read.table("data/5ff4c4c88319c6a96e14c5b061e3c91f-schachtman_2017_soil.tsv",header=T,stringsAsFactors=F,sep="\t")
f17biom <- read.table("data/5ff4c4c88319c6a96e14c5b061e3c91f-schachtman_2017_biom.tsv",header=T,stringsAsFactors=F,sep="\t")
f17iso <- read.table("data/4dd9c063253b200b2ca408cc25e346c9-cousins_2017_field.tsv",header=T,stringsAsFactors=F,sep="\t")
metab <- read.csv("data/SpecAbundGROUPS_SM_Harvest_LCMS_4db.csv",header = T, stringsAsFactors = T)

ui <- dashboardPage(skin="black", title="DOE Sorghum Systems",
                    dashboardHeader(
                      title = tagList(
                        tags$span(
                          class = "logo-mini", "DOE SS"
                        ),
                        tags$span(
                          class = "logo-lg", "DOE Sorghum Systems"
                        )
                      ),
                      titleWidth = 450
                    ),
                    dashboardSidebar(
                      tags$script(HTML("$('body').addClass('sidebar-mini');")),
                      width = 150,
                      sidebarMenu(
                        shinyjs::useShinyjs(),
                        id="tabs",
                        menuItem(text="2017 Fields",tabName = "2017fields",icon = icon("sitemap")),
                        menuItem(text="2016 Fields", tabName = "2016fields")
                      )
                    ),
                    dashboardBody(
                      useShinyjs(),
                      tags$head(tags$style("#container * {display: inline;}")),
                      tags$style(HTML("
                                        .shiny-progress-notification .progress-bar {background-color: darkgreen;}
                                        .shiny-progress-notification .progress-text {
                                        font-size: 17pt;
                                        }
                                        .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
                                        .multicol{
                                        -webkit-column-count: 4; /* Chrome, Safari, Opera */
                                        -moz-column-count: 4; /* Firefox */
                                        column-count: 4;
                                        }
                                        .twocol{
                                        -webkit-column-count: 2; /* Chrome, Safari, Opera */
                                        -moz-column-count: 2; /* Firefox */
                                        column-count: 2;
                                        }
                                        .warning { 
                                        color: red;
                                        }"
                      )),
                      tabItems(
                        tabItem(tabName = "2016fields",
                                fluidRow(
                                  box(style = "overflow-y:scroll",width=12,title = "Harvest Traits",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=T,
                                      column(width=3,
                                             selectInput("which_trait_16","Which Trait",c("fresh_weight_stalk_g","fresh_weight_panicle_g","total_fresh_weight_kg_hectare","dry_weight_stalk_g","dry_weight_panicle_g","total_dry_weight_kg_hectare","number_of_stalks","plant_height_cm"),"plant_height_cm",width=180)
                                      ),
                                      column(width=3,
                                             selectInput("which_harvest_16","Which Harvest",c("small 1","large","small 2"),"large",width=180)
                                      ),
                                      column(width=12,
                                             uiOutput("harvest_plot_16")
                                      )
                                  )
                                )
                        ),
                        tabItem(tabName = "2017fields",
                                           fluidRow(
                                             box(style = "overflow-y:scroll",width=12,title = "Harvest Traits",solidHeader = T,status = 'success',collapsible = T,collapsed=T,
                                                 tabsetPanel(
                                                   tabPanel(title="Boxplots",
                                                            column(width=3,
                                                              selectInput("which_harvest_trait_17","Which Trait",c("fresh_weight_stalk_g","fresh_weight_panicle_g","total_fresh_weight_kg_hectare","dry_weight_stalk_g","dry_weight_panicle_g","total_dry_weight_kg_hectare","number_of_stalks","plant_height_cm"),"plant_height_cm",width=180),
                                                              uiOutput("which_harvest_17_date"),
                                                              tags$b("Highlight Genos"),
                                                              uiOutput("which_harvest_17_genos"),
                                                            ),
                                                            column(width=9,
                                                                   plotOutput("harvest_plot_17")
                                                            )
                                                   ),
                                                   tabPanel(title="Spatial",
                                                            p("test")
                                                   )
                                                 )
                                             ),
                                             box(style = "overflow-y:scroll",width=12,title = "Soil Nutrients",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=T,
                                                 tabsetPanel(
                                                   tabPanel(title="Boxplots",
                                                            column(width=3,
                                                                   selectInput("which_soil_trait_17","Which Trait",c("Soil_pH","WDRF_Buffer_pH","S_Salts_mmho.cm","Texture_No","Organic_Matter_LOI_.","Nitrate.N_ppm_N","lbs_N.A","Potassium._ppm_K","Sulfate.S_ppm_S","Calcium_ppm_Ca","Magnesium_ppm_Mg","Sodium_ppm_Na","CEC.Sum_of_Cations_me.100g","Mehlich_P.III_ppm_P"),"Soil_pH",width=180),
                                                                   uiOutput("which_soil_17_date"),
                                                                   tags$b("Highlight Genos"),
                                                                   uiOutput("which_soil_17_genos"),
                                                            ),
                                                            column(width=9,
                                                                   plotOutput("soil_plot_17")
                                                            )
                                                   ),
                                                   tabPanel(title="Spatial",
                                                            p("test")
                                                   )
                                                 )
                                             ),
                                             box(style = "overflow-y:scroll",width=12,title = "Isotopes",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=T,
                                                 tabsetPanel(
                                                   tabPanel(title="Boxplots",
                                                            column(width=3,
                                                                   selectInput("which_iso_trait_17","Which Trait",c("lma_gm2","d13c_per_il","d15n_per_mil","cent_n","cent_c","sla_m2g"),"lma_gm2",width=180),
                                                                   uiOutput("which_iso_17_date"),
                                                                   tags$b("Highlight Genos"),
                                                                   uiOutput("which_iso_17_genos"),
                                                            ),
                                                            column(width=9,
                                                                   plotOutput("iso_plot_17")
                                                            )
                                                   ),
                                                   tabPanel(title="Spatial",
                                                            p("test")
                                                   )
                                                 )
                                             ),
                                             box(style = "overflow-y:scroll",width=12,title = "Microbiome",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=T,
                                                 column(width=12,
                                                        tags$b("Overview"),
                                                        p("2017 Fields")
                                                 )
                                             ),
                                             box(style = "overflow-y:scroll",width=12,title = "Metabolomics",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=T,
                                                 tabsetPanel(
                                                   tabPanel(title="ANOSIM",
                                                            uiOutput("meta_parse"),
                                                            hr(),
                                                            actionButton("use_config","Use subset"),
                                                            uiOutput("pca_output")  
                                                   )
                                                 )
                                             )
                                           )
                            )
                        )
                    )
)

server = function(input, output, session){
  options(shiny.maxRequestSize=5000*1024^2)
  
  #**********************************************************************************************************************
  # 2017 Harvest Data
  #**********************************************************************************************************************
  output$which_harvest_17_date <- renderUI({
    sels <- as.character(na.omit(unique(c(f17biom$biomass_date[!is.na(f17biom[,input$which_harvest_trait_17])],f17biom$plant_height_date[!is.na(f17biom[,input$which_harvest_trait_17])]))))
    selectInput("which_harvest_17_date_sel","Which Harvest",sels,sels[1],width=180)
  })
  
  output$which_harvest_17_genos <- renderUI({
    sels <- unique(f17biom$geno[f17biom$biomass_date == input$which_harvest_17_date_sel | f17biom$plant_height_date == input$which_harvest_17_date_sel])
    names(sels) <- sels
    tags$div(class = "multicol",checkboxGroupInput("which_harvest_17_geno_sel", NULL,sels,inline = F))
  })

  output$harvest_plot_17 <- renderPlot({
    sub <- f17biom[f17biom$biomass_date == input$which_harvest_17_date_sel | f17biom$plant_height_date == input$which_harvest_17_date_sel,]
    sub <- sub[!is.na(sub[,input$which_harvest_trait_17]),]
    if(nrow(sub) > 0){
      highlights <- sub[sub$geno %in% input$which_harvest_17_geno_sel,]
      ggplot(data=sub,aes_string("geno",input$which_harvest_trait_17))+
        facet_wrap(~treatment)+
        geom_boxplot()+
        geom_boxplot(data=highlights,aes(fill=geno))+
        theme_light()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.title.y= element_text(size = 18),
              axis.title.x= element_blank())+
        theme(axis.text.y = element_text(size = 14),
              axis.text.x = element_text(size = 10))+
        theme(axis.ticks.length=unit(0.2,"cm"))
    }else{
      ggplot()
    }
  })
  
  
  #**********************************************************************************************************************
  # 2017 Soil Data
  #**********************************************************************************************************************
  output$which_soil_17_date <- renderUI({
    sels <- na.omit(unique(f17soil$sample_date[!is.na(f17soil[,input$which_soil_trait_17])]))
    selectInput("which_soil_17_date_sel","Which Date",sels,sels[1],width=180)
  })
  
  output$which_soil_17_genos <- renderUI({
    sels <- unique(f17soil$geno[f17soil$sample_date == input$which_soil_17_date_sel])
    names(sels) <- sels
    tags$div(class = "multicol",checkboxGroupInput("which_soil_17_geno_sel", NULL,sels,inline = F))
  })
  
  output$soil_plot_17 <- renderPlot({
    sub <- f17soil[f17soil$sample_date == input$which_soil_17_date_sel,]
    sub <- sub[!is.na(sub[,input$which_soil_trait_17]),]
    if(nrow(sub) > 0){
      highlights <- sub[sub$geno %in% input$which_soil_17_geno_sel,]
      ggplot(data=sub,aes_string("geno",input$which_soil_trait_17))+
        facet_wrap(~treatment)+
        geom_boxplot()+
        geom_boxplot(data=highlights,aes(fill=geno))+
        theme_light()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.title.y= element_text(size = 18),
              axis.title.x= element_blank())+
        theme(axis.text.y = element_text(size = 14),
              axis.text.x = element_text(size = 10))+
        theme(axis.ticks.length=unit(0.2,"cm"))
    }else{
      ggplot()
    }
  })
  
  
  #**********************************************************************************************************************
  # 2017 Isotope Data
  #**********************************************************************************************************************
  output$which_iso_17_date <- renderUI({
    sels <- na.omit(unique(f17iso$harvest.date[!is.na(f17iso[,input$which_iso_trait_17])]))
    selectInput("which_iso_17_date_sel","Which Date",sels,sels[1],width=180)
  })
  
  output$which_iso_17_genos <- renderUI({
    sels <- unique(f17iso$genotype[f17iso$harvest.date == input$which_iso_17_date_sel])
    names(sels) <- sels
    tags$div(class = "multicol",checkboxGroupInput("which_iso_17_geno_sel", NULL,sels,inline = F))
  })
  
  output$iso_plot_17 <- renderPlot({
    sub <- f17iso[f17iso$harvest.date == input$which_iso_17_date_sel,]
    sub <- sub[!is.na(sub[,input$which_iso_trait_17]),]
    if(nrow(sub) > 0){
      highlights <- sub[sub$genotype %in% input$which_iso_17_geno_sel,]
      ggplot(data=sub,aes_string("genotype",input$which_iso_trait_17))+
        facet_wrap(~treatment)+
        geom_boxplot()+
        geom_boxplot(data=highlights,aes(fill=genotype))+
        theme_light()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.title.y= element_text(size = 18),
              axis.title.x= element_blank())+
        theme(axis.text.y = element_text(size = 14),
              axis.text.x = element_text(size = 10))+
        theme(axis.ticks.length=unit(0.2,"cm"))
    }else{
      ggplot()
    }
  })
  
  
  #**********************************************************************************************************************
  # 2017 Metabolomics Data
  #**********************************************************************************************************************
  output$meta_parse <- renderUI({
    if(!is.null(metab)){
      dataTableOutput("meta_config")
    }
  })
  
  output$meta_config <- renderDataTable({
    datatable(metab[,1:11],rownames = F,selection = "none",filter = 'top',options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 5, autoWidth = TRUE))
  })
  
  meta_sub <- reactiveValues(data = NULL)
  observeEvent(input$use_config,{
    rows <- input$meta_config_rows_all
    meta_sub$data <- metab[rows,]
    
  })
  
  output$pca_output <- renderUI({
    if(!is.null(meta_sub$data)){
       tagList(
         hr(),
         column(width=4,
            column(width = 6,
                   selectInput("simca_grouping_1",width = 280,label = "Grouping 1",choices = colnames(meta_sub$data[1:11]),selected = "genotype"),
                   selectInput("simca_x_axis",width = 280, label = "X-Axis",choices = c("Component 1", "Component 2", "Component 3", "Component 4", "Component 5", "Component 6", "Component 7", "Component 8", "Component 9", "Component 10", "Component 11", "Component 12", "Component 13", "Component 14", "Component 15"), selected = "Component 1")
            ),
            column(width = 6,
                   selectInput("simca_grouping_2",width = 280,label = "Grouping 2",choices = colnames(meta_sub$data[1:11]),selected = "treatment"),
                   selectInput("simca_y_axis",width = 280, label = "Y-Axis",choices = c("Component 1", "Component 2", "Component 3", "Component 4", "Component 5", "Component 6", "Component 7", "Component 8", "Component 9", "Component 10", "Component 11", "Component 12", "Component 13", "Component 14", "Component 15"), selected = "Component 2")
            ),
            textInput("simca_plot_name","Plot Name (ex: E17_leafOnly_cbTreatment.png)",width=400),
            downloadButton("simca_download","Download Plot"),
            br(),
            downloadButton("data_download","Download Data (tsv)")
         ),
         column(width = 8,
            plotOutput("simca_plot_ind")
         )
      )
    }
  })
  
  x <- reactive({
    as.numeric(strsplit(input$simca_x_axis," ")[[1]][2])
  })
  
  y <- reactive({
    as.numeric(strsplit(input$simca_y_axis," ")[[1]][2])
  })
  
  simca_df <- reactive({
    x <- x()
    y <- y()
    comps <- meta_sub$data
    simca <- simca(comps[,12:ncol(comps)],"blah",scale = T)
    df <- data.frame(predict(simca, comps[,12:ncol(comps)])$scores[,1:15],stringsAsFactors = F)
    df$Grouping1 <- comps[,input$simca_grouping_1]
    df$Grouping2 <- comps[,input$simca_grouping_2]
    varexp <- signif(100*simca$eigenvals[c(x,y)]/sum(simca$eigenvals),4)
    list(df,varexp)
  })
  
  simca_data <- reactive({
    cbind(meta_sub$data[,1:11],simca_df()[[1]])
  })
  
  
  simca_plot <- reactive({
    x <- x()
    y <- y()
    simca_df <- simca_df()
    simca_df[[1]]$Grouping1 <- as.factor(simca_df[[1]]$Grouping1)
    simca_df[[1]]$Grouping2 <- as.factor(simca_df[[1]]$Grouping2)
    ggplot(data=simca_df[[1]], aes_string(colnames(simca_df[[1]])[x],colnames(simca_df[[1]])[y]))+
      stat_ellipse(color="gray40",fill="white",geom = "polygon",type = "norm")+
      geom_point(aes(color=Grouping2,shape=Grouping1),alpha=0.6,size=4)+
      xlab(paste(gsub("Comp.","Component ",colnames(simca_df[[1]])[x])," (",simca_df[[2]][1],"%)",sep = ""))+
      ylab(paste(gsub("Comp.","Component ",colnames(simca_df[[1]])[y])," (",simca_df[[2]][2],"%)",sep = ""))+
      geom_vline(xintercept = 0,linetype="dashed")+
      geom_hline(yintercept = 0,linetype="dashed")+
      theme_bw()+
      theme(panel.background = element_rect(fill = "gray95"))+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(axis.ticks.length=unit(0.2,"cm"))+
      theme(aspect.ratio=1)
  })
  
  output$simca_plot_ind <- renderPlot({
    simca_plot()
  })
  
  output$simca_download <- downloadHandler(
    filename = function() {input$simca_plot_name},
    content=function(file){
      ggsave(file,simca_plot(),device = "png",width = 6,height = 5,dpi = 300)
    }
  )
  
  output$data_download <- downloadHandler(
    filename = function() {"metabolomics_simca_data.tsv"},
    content = function(file){
      head(simca_data())
      write.table(simca_data(),file,row.names = FALSE, quote = FALSE,sep = "\t")
    }
  )
}

shinyApp(ui, server)
                    