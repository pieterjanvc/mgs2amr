library(tidyverse)
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(RSQLite)
library(shinyFiles)
library(RColorBrewer)

ui <- dashboardPage(
  

  dashboardHeader(title = "MGS2AMR explorer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Detected ARG", tabName = "detectedARG", icon = icon("dna")),
      menuItem("Detected Bacteria", tabName = "detectedBact", icon = icon("bacterium"))
    )
  ),
  
  dashboardBody(
    
    #Add CSS to hover color sankey diagram
    tags$style(HTML("
      .sankey-link:hover {
      fill: gold !important;
      fill-opacity: 1 !important;
      }
    ")),
    
    #Tabs
    tabItems(
      
      tabItem(tabName = "dashboard", 
              h2("Import Data"),
              fluidRow(
              column(12,wellPanel(
                     shinyFilesButton('linkDB', 
                                      label='Connect to local MGS2AMR database or upload an output tar.gz file', 
                                      title='Please select a .db or tar.gz file', 
                                      multiple=FALSE),br(),br(),
                     dataTableOutput("pipelineTable")
              ))),
              fluidRow(
                infoBoxOutput("fileInfo", 12)
              ),
              fluidRow(column(12,
                plotlyOutput("ARGoverview")
              ))
              ),
      
      tabItem(tabName = "detectedARG", 
              h2("Detected ARG"),
              fluidRow(column(12,dataTableOutput("ARGtable"))),
              uiOutput("ARGtoBact_title"),
              fluidRow(column(12,dataTableOutput("ARGtoBact"))),
              uiOutput("ARGstrains_title"),
              fluidRow(column(12,dataTableOutput("ARGstrains")))
              ),
      
      tabItem(tabName = "detectedBact", 
              h2("Top 25 High Scoring Bacteria"),
              fluidRow(column(12,plotlyOutput("bactHeatmap", height = 600))),
              fluidRow(column(12,dataTableOutput("bactTable"))),
              uiOutput("bactARG_title"),
              fluidRow(column(12,plotlyOutput("bactPlot"))),
              fluidRow(column(12,dataTableOutput("bactARG")))
              )
      
    )
  )
  
)

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=500*1024^2)
  
  # ---- VARIABLES ---- 
  #*********************
  
  ARGtoBact_title = reactiveVal()
  ARGstrains_title = reactiveVal()
  tempDir = reactiveVal()
  
  # ---- FUNCTIONS ---- 
  #*********************
  
  createLink <- function(link, value, btn = F) {
    sprintf('<a href="%s" target="_blank" %s>%s</a>',
            link, ifelse(btn, '"btn btn-primary"', ""), value)
  }
  
  alphaToHex = function(alpha){
    hex = as.integer(alpha*255) %>% as.hexmode() %>% as.character()
    return(ifelse(nchar(hex) == 1, paste0("0", hex), hex))
  }

  # ---- INPUT DATA ---- 
  #*********************
  
  #Connect to a local MGS2AMR database - fs::path_home()
  volumes <- c(Home = "~", getVolumes()())
  shinyFileChoose(input, "linkDB", roots = volumes, session = session, filetypes=c("gz","db"))
  
  # Show all completed pipeline runs in a table
  pipelineTable = reactive({
    req(!is.integer(input$linkDB))
    
    filePath = parseFilePaths(volumes, input$linkDB)
    
    if(str_detect(filePath$name, "db$")){
      myConn = dbConnect(SQLite(), filePath$datapath)
      myTables = dbListTables(myConn)
      check = all(c("pipeline", "detectedARG", "annotation") %in% myTables)
      if(!check){
        showModal(modalDialog("The database is not in valid MGS2AMR format"))
      }
      req(check)
      
      pipelineTable = tbl(myConn, "pipeline") %>% filter(statusCode > 3) %>%
        select(pipelineId, name, start = startTimestamp, modified = modifiedTimestamp) %>% 
        collect()
      dbDisconnect(myConn)
    } else {
      
      fileName = str_match(filePath$datapath, "([^/]+)\\.tar\\.gz$")[2]
      td = tempdir()
      untar(filePath$datapath, exdir = td)
      td = sprintf("%s/%s", td, fileName)
      myFiles = list.files(td)
      
      if(length(myFiles) == 0 | !all(str_detect(myFiles, "info|detectedARG|annotation"))){
        showModal(modalDialog("The upoaded file is not an MGS2AMR output file"))
      }
      
      req(length(myFiles) != 0 & all(str_detect(myFiles, "info|detectedARG|annotation")))
      
      isolate(tempDir(list(dir = td, files = myFiles, name = fileName)))
      
      pipelineTable = read_csv(sprintf("%s/%s", td, myFiles[str_detect(myFiles, "info")]),
                               show_col_types = FALSE) %>% 
        select(pipelineId, name, start = startTimestamp, modified = modifiedTimestamp)
    }
    
    pipelineTable
    
  })
  
  output$pipelineTable = renderDataTable({
    
    datatable(pipelineTable(), rownames = F, selection = 'single')
    
  })


  #Get all data for a specific pipeline run
  myData = reactive({
    
    req(!is.null(input$pipelineTable_rows_selected))

    if(!is.null(tempDir())){

      myInfo = read_csv(sprintf(
        "%s/%s", tempDir()$dir, tempDir()$files[str_detect(tempDir()$files, "info")]),
        show_col_types = FALSE)
      detectedARG = read_csv(sprintf(
        "%s/%s", tempDir()$dir, tempDir()$files[str_detect(tempDir()$files, "detectedARG")]),
        show_col_types = FALSE)
      output = read_csv(sprintf(
        "%s/%s", tempDir()$dir, tempDir()$files[str_detect(tempDir()$files, "annotation")]),
        show_col_types = FALSE)

    } else {

      DBpath = parseFilePaths(volumes, input$linkDB)

      myConn = dbConnect(SQLite(), DBpath$datapath)
      
      myInfo = tbl(myConn, "pipeline") %>%
        filter(pipelineId == local(
          pipelineTable()$pipelineId[input$pipelineTable_rows_selected])) %>% 
        collect()
      
      detectedARG = tbl(myConn, "detectedARG") %>%
        filter(pipelineId == local(
          pipelineTable()$pipelineId[input$pipelineTable_rows_selected])) %>% 
        # select(-gene, -subtype) %>% 
        collect() %>% 
        left_join(tbl(myConn, "ARG"), by = "geneId") %>% 
        mutate(
          cover = ifelse(type == "noFragments", cover1, cover2),
          across(c(cover, startPerc), function(x) round(x * 100, 2)),
          startDepth = round(startDepth, 2)) %>% 
        arrange(desc(cover), desc(startDepth))
      
      output = tbl(myConn, "annotation") %>%
        filter(pipelineId == local(
          pipelineTable()$pipelineId[input$pipelineTable_rows_selected])) %>%
        left_join(tbl(myConn, "bactStrains"), by = "accession") %>%
        left_join(tbl(myConn, "bactTaxa"), by = "taxid") %>% 
        collect() %>% 
        group_by(geneId) %>% mutate(
          top = fullPath / max(fullPath),
          depth = KC / LN
        ) %>% ungroup() %>% left_join(
          detectedARG %>% select(geneId, gene, subtype, cover, type), 
          by = "geneId") %>% 
        mutate(plasmid = as.logical(plasmid))
      
      dbDisconnect(myConn)

    }

    req(all(c('pipelineId','name','outputFolder', #change to output folder later
              'inputfileBP','startTimestamp') %in% colnames(myInfo)))
    
    req(all(c('geneId','accession','taxid','genus','species','plasmid',
              'LN','KC','extension','fullPath','extraBits','path0',
              'maxOrder0','path1','maxOrder1','gene','subtype','cover','type',
              'depth','top','pipelineId') %in% colnames(output)))
    
    removeDuplicates = output %>% left_join(
      detectedARG %>% group_by(ARGgroup) %>% 
        filter(n() > 1) %>% ungroup() %>% 
        select(geneId, ARGgroup), by = "geneId") %>% 
      filter(!is.na(ARGgroup))
    
    removeDuplicates = removeDuplicates %>% 
      group_by(accession, geneId, ARGgroup) %>% 
      summarise(total = sum(fullPath), .groups = "drop") %>% 
      group_by(geneId, ARGgroup) %>% 
      summarise(total = max(total), .groups = "drop") %>% 
      group_by(ARGgroup) %>% 
      filter(geneId != geneId[total == max(total)][1]) %>% ungroup()
    
    detectedARG = detectedARG %>% filter(!geneId %in% removeDuplicates$geneId)
    output = output %>% filter(!geneId %in% removeDuplicates$geneId)
    
    #Only keep genera with at least one top hit across all ARG
    toKeep = output %>% filter(top == 1) %>% 
      group_by(genus) %>% summarise(.groups = "drop")
    
    output = output %>% filter(genus %in% toKeep$genus)
    
    return(list(myInfo = myInfo, detectedARG = detectedARG, output = output))
    
  })
  
  #Populate metadata no run in an infobox
  output$fileInfo = renderInfoBox({
    req(myData())

    infoBox(
      sprintf("Pipeline ID %s (%s)", myData()$myInfo$pipelineId, 
              myData()$myInfo$startTimestamp),
      HTML(sprintf("%s<h5 style='word-wrap: break-word';>%s</h5>",
                   myData()$myInfo$name, myData()$myInfo$outputFolder)),
      icon = icon("microscope"), color = "blue", fill = T
    )
  })
  
  
  # ---- DETECTED ARG ---- 
  #***********************
 
  #Output all detected ARG
  output$ARGtable = renderDataTable({
    datatable(myData()$detectedARG %>% mutate(
      gene = createLink(paste0("https://www.ncbi.nlm.nih.gov/nucleotide/", nucl), gene)) %>% 
        select(gene, allele = subtype, class, name = info, `length (bp)` = nBases,
               `coverage (%)` = cover, `identity (%)` = startPerc, depth = startDepth),
      rownames = F, selection = 'single', escape = F)
  })
  
  #Show an overview of detected ARG as treemap
  output$ARGoverview = renderPlotly({
    
    myPlot = myData()$detectedARG %>% select(id = gene, parent = class) %>% distinct()
    myPlot = bind_rows(data.frame(id = unique(myPlot$parent)), myPlot) %>% 
      mutate(parent = replace_na(parent, "ARG"))
    myPlot = bind_rows(data.frame(id = "ARG", parent = ""), myPlot)
    myPlot = myPlot %>% mutate(
      label = sprintf("<b>%s</b>", id),
      label = ifelse(label == "<b>ARG</b>", "<b>Antimicrobial Resistance Genes</b>", label))
    
    plot_ly(
      type='treemap',
      ids=myPlot$id,
      labels=myPlot$label,
      parents=myPlot$parent,
      maxdepth=2,
      textposition = "middle center")
    
  })
  
  #Bact per ARG table
  ARGtoBact = reactive({
    myData()$output %>% 
      filter(gene %in% myData()$detectedARG$gene[input$ARGtable_rows_selected]) %>% 
      group_by(taxid, plasmid) %>% 
      filter(fullPath == max(fullPath)) %>% slice(1) %>% ungroup() %>% 
      mutate(across(c(depth, extension, fullPath), round, digits = 2)) %>% 
      arrange(desc(top), desc(fullPath))
  })
  
  output$ARGtoBact = renderDataTable({
    if(is.null(input$ARGtable_rows_selected)){
      
      ARGtoBact_title(NULL)
      NULL
    } else {
      
      ARGtoBact_title(sprintf(
        "Bacterial matches for %s",
        myData()$detectedARG$gene[input$ARGtable_rows_selected]
      ))
      
      datatable(ARGtoBact() %>% mutate(
        taxid = createLink(paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", taxid), taxid)
      ) %>% 
        select(genus, species, taxid, plasmid, `path score` = fullPath, 
               depth, type, top), 
      escape = F, rownames = F, selection = 'single')
    }

  })
  
  output$ARGtoBact_title = renderUI(h2(ARGtoBact_title()))
  
  #Strains for ARG of one bact
  ARGstrains = reactive({
    myData()$output %>% 
      filter(taxid %in% ARGtoBact()$taxid[input$ARGtoBact_rows_selected],
             plasmid == ARGtoBact()$plasmid[input$ARGtoBact_rows_selected]) %>% 
      mutate(across(c(depth, extension, fullPath), round, digits = 2)) %>% 
      arrange(desc(top), desc(fullPath))
       
  })
  
  output$ARGstrains = renderDataTable({
    
    if(is.null(input$ARGtoBact_rows_selected)){
      ARGstrains_title(NULL)
      NULL
    } else {
      
      ARGstrains_title(sprintf(
        "All matching strains for %s in %s %s",
        myData()$detectedARG$gene[input$ARGtable_rows_selected], 
        ARGtoBact()$genus[input$ARGtoBact_rows_selected],
        ARGtoBact()$species[input$ARGtoBact_rows_selected]
      ))
      
      datatable(ARGstrains() %>% mutate(
        accession = createLink(paste0("https://www.ncbi.nlm.nih.gov/nuccore/", accession), accession)) %>%
          select(genus, species, accession, plasmid, 
                 `path score` = fullPath, depth, type, top) ,
        rownames = F, selection = 'single', escape = F)
    }
    
  })
  
  output$ARGstrains_title = renderUI(h2(ARGstrains_title()))
  
  
  # ---- DETECTED BACTERIA ---- 
  #****************************
  
  #Heatmap of top 25 species
  output$bactHeatmap = renderPlotly({
    
    heatmapData = myData()$output %>% group_by(genus, gene) %>% 
      filter(fullPath == max(fullPath)) %>% 
      slice(1) %>% 
      group_by(genus) %>% summarise(totalScore = sum(fullPath), .groups = "drop") %>% 
      arrange(desc(totalScore)) %>% slice(1:25)
    
    heatmapData = myData()$output %>% 
      filter(genus %in% heatmapData$genus) %>% 
      group_by(genus, gene) %>% 
      filter(fullPath == max(fullPath)) %>% 
      slice(1) %>% 
      group_by(geneId) %>% mutate(fullPath = fullPath / max(fullPath)) %>% 
      ungroup() %>% select(genus, gene, fullPath) %>% 
      pivot_wider(id_cols = genus, names_from = "gene", values_from = fullPath,
                  values_fill = 0) %>% as.data.frame()
    
    rownames(heatmapData) = heatmapData$genus
    
    if(nrow(heatmapData) > 1){
      hm = heatmap(heatmapData %>% select(-genus) %>% as.matrix())
      heatmapData = heatmapData %>% select(-genus) %>% as.matrix()
      heatmapData = heatmapData[hm$rowInd, hm$colInd]
    } else {
      heatmapData = heatmapData %>% select(-genus) %>% as.matrix()
    }
    
    if(all(heatmapData[1,1] == heatmapData)){
      colors = "#8080ff"
    } else {
      colors = colorRamp(c("white", "blue"))
    }
    
    plot_ly(y = rownames(heatmapData), x = colnames(heatmapData), 
            z = heatmapData, type = "heatmap", xgap = 1, ygap = 1,
            colors = colors, zmin = 0, zmax = 1)
    
  })
  
  #Most likely bacteria (at genus level, ignoring plasmids)
  bactTable = reactive({

    bactTable = data.frame()
    tempVar = myData()$output %>%
      mutate(plasmid = F) %>% 
      filter(!plasmid, extension > 1) %>%
      group_by(accession, gene) %>%
      filter(geneId == geneId[fullPath == max(fullPath)][1]) %>%
      ungroup()

    while(nrow(tempVar) > 0){
      bact = tempVar %>% group_by(accession, taxid, genus, species, plasmid) %>%
        summarise(
          totalScore = sum(fullPath), nGenes = n_distinct(geneId),
          gDepth = mean(depth[!plasmid],0, na.rm = T),
          pDepth = mean(depth[plasmid],0, na.rm = T),
          nExt = sum(extension > 1), .groups = "drop") %>%
        arrange(desc(totalScore)) %>%
        group_by(genus) %>% filter(totalScore == max(totalScore)) %>% ungroup()

      tempVar = tempVar %>%
        filter(!geneId %in% geneId[accession == bact$accession[1]],
               genus != bact$genus[1])

      #Check if there are muliple top genera, if so do not report
      if(bact %>% filter(totalScore == max(totalScore)) %>%
         pull(genus) %>% n_distinct() > 1){
        next
      }

      #If multiple species are at the top, label as mixed
      if(bact %>% filter(totalScore == max(totalScore)) %>%
         pull(species) %>% n_distinct() > 1){
        bact[1,"species"] = "mixed"
      }

      bactTable = bind_rows(bactTable, bact[1,])
    }

    bactTable %>%
      mutate(across(c(totalScore, gDepth), function(x){round(x, digits = 0)}))
    
    # bactTable = testData$output %>%
    #   group_by(genus, species) %>% 
    #   mutate(spScore = sum(fullPath)) %>% 
    #   group_by(genus) %>% 
    #   mutate(species = species[spScore == max(spScore)][1]) %>% 
    #   group_by(genus, gene) %>% 
    #   filter(fullPath == max(fullPath)) %>% 
    #   slice(1) %>% 
    #   group_by(genus, species) %>% 
    #   summarise(totalScore = sum(fullPath), 
    #             genomeDepth = max(mean(depth[!plasmid],0, na.rm = T), 0),
    #             plasmidDepth = max(mean(depth[plasmid],0, na.rm = T), 0),
    #             .groups = "drop") %>% 
    #   arrange(desc(totalScore)) 
    
  })
  
  output$bactTable = renderDataTable({

    datatable(bactTable() %>%
                select(genus, `closest species` = species, score = totalScore,
                       `# ARG` = nGenes, `estimated depth` = gDepth),
              rownames = F, selection = 'single')
    
    # datatable(bactTable() %>% 
    #             select(genus, `closest species` = species, score = totalScore,
    #                    `estimated  genome depth` = genomeDepth,
    #                    `estimated  plasmid depth` = plasmidDepth),
    #           rownames = F, selection = 'single')
  })
  
  
  #List of ARG for selected bact
  bactARG = reactive({
    
    myData()$output %>% 
    filter(genus == bactTable()$genus[input$bactTable_rows_selected]) %>%
    group_by(geneId) %>% 
    mutate(species = ifelse(
      n_distinct(species[top == max(top)]) > 1,
      paste(n_distinct(species[top == max(top)]), " top species"),
      species)) %>% 
    filter(top == max(top)) %>%  
    slice(1) %>% 
    ungroup() %>% 
    mutate(
      perc = fullPath / max(fullPath),
      across(c(top, perc, fullPath, extension), round, digits = 2)
    ) %>% 
    distinct() %>% 
    arrange(desc(top), desc(perc))
    
  })
  
  output$bactARG = renderDataTable({
    
    if(is.null(input$bactTable_rows_selected)){
      ARGstrains_title(NULL)
      NULL
    } else {
      
      ARGstrains_title(sprintf(
        "All matching ARG for <i>%s (%s)</i>",
        bactTable()$genus[input$bactTable_rows_selected],
        bactTable()$species[input$bactTable_rows_selected]
      ) %>% HTML())
      
      datatable(bactARG() %>%
                  select(genus, species, gene, subtype, plasmid,
                         taxid, fullPath, extension, type, top, perc), 
                rownames = F, selection = 'single')
    }
    
  })
  
  output$bactARG_title = renderUI(h2(ARGstrains_title()))
  
  #Graphically show the ARG associated with a bacterium 
  # (and other bacteria containing this ARG as a top hit)
  output$bactPlot = renderPlotly({
    
    if(is.null(input$bactTable_rows_selected)){
      NULL
    } else {
      
      myGenus = bactTable()$genus[input$bactTable_rows_selected]
    
      links = myData()$output %>% 
        filter(gene %in% bactARG()$gene, genus %in% bactTable()$genus) %>% 
        group_by(genus, gene) %>% filter(fullPath == max(fullPath)) %>% 
        slice(1) %>% ungroup() %>% select(genus, gene, plasmid, top) %>% 
        distinct() %>% 
        left_join(
          myData()$detectedARG %>% select(gene, class) %>% distinct(),  
        by = "gene") %>% 
        filter(genus == myGenus | top == 1)
        
      nodes = data.frame(
        name = sort(unique(links$genus)),
        id = as.factor(sort(unique(links$genus))) %>% as.integer() - 1
      ) %>% mutate(
        color = ifelse(name == myGenus, "#825cdb", "#6363614D"),
        border = ifelse(name == myGenus, "black", "white"))
      
      nodes = bind_rows(
        nodes,
        
        links %>% select(name = gene, class) %>% distinct() %>% arrange(name) %>% 
          mutate(
            id = as.factor(name) %>% as.integer() + max(nodes$id),
            color = brewer.pal(
              max(n_distinct(class),3), 
              "BrBG")[1:n_distinct(class)][as.factor(class) %>% as.integer()],
            border = "white"),

        data.frame(
          name = as.character(sort(unique(links$class))),
          id = as.factor(sort(unique(links$class))) %>% 
            as.integer() + n_distinct(links$genus) + n_distinct(links$gene) -1
        ) %>% mutate(color = brewer.pal(max(n(),3), "BrBG")[1:n()], border = "black")
      )
      
      links = bind_rows(
        
        links %>%
          left_join(nodes %>% select(genus = name, source = id), by = "genus") %>%
          left_join(nodes %>% select(gene = name, target = id), by = "gene") %>%
          mutate(
            color = case_when(
              genus != myGenus ~ "#d1d0c24D",
              plasmid & top == 1 ~ "#F24EA6",
              plasmid ~ "#f07ab9",
              !plasmid & top == 1 ~ "#08B3E7",
              TRUE ~ "#4fc1e3"),
            color = ifelse(
              genus != myGenus, color,
               paste0(color, alphaToHex((exp(1*top) - 1) / (exp(1) - 1))))
          ) %>% 
          select(source, target, color),
        
        links %>% 
          left_join(nodes %>% select(gene = name, source = id), by = "gene") %>% 
          left_join(nodes %>% select(class = name, target = id), by = "class") %>% 
          mutate(color = "#d1d0c24D") %>% select(source, target, color)
        
      )

      myPlot = plot_ly(
        type = "sankey",
        orientation = "h",
        
        node = list(
          label = nodes$name,
          color = nodes$color,
          pad = 10,
          thickness = 10,
          line = list(
            color = "white",
            width = 0.5
          )
        ),
        
        link = data.frame(
          source = links$source,
          target = links$target,
          color = links$color,
          value =  2
        )
      )
      
      myPlot
      
    }
  })
  
}

shinyApp(ui, server)
