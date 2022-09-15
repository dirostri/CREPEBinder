require(shiny)
require(shinyFiles)
require(shinydashboard)
require(shinyjs)
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(phytools))
suppressPackageStartupMessages(require(phylotools))
suppressPackageStartupMessages(require(ape))
options(shiny.maxRequestSize = 100*1024^2)
ui <- ui <- dashboardPage(skin = "black",
                          dashboardHeader(title = 'CREPE TF Cataloguing and Extraction', titleWidth = 450),
                          dashboardSidebar(width = 350,
                                           selectInput(inputId = "select", 
                                                       label = "Choose Analysis To Run", 
                                                       choices = c('TF Cataloguing', 'TF Annotation', 'Custom')),
                                           conditionalPanel(
                                             condition = "input.select == 'TF Cataloguing'",
                                             fileInput(inputId = "file1",
                                                       label = "Protein Fasta",
                                                       accept = c(".fa", ".fasta")
                                             ),
                                             #fileInput(inputId = "file2",
                                             #label = "PFAM HMM path",
                                             #accept = c(".hmm")
                                             #),
                                             textInput(inputId = 'species', label = 'Species Name'),
                                             h4("Citations for Tools and Dependencies"),
                                             h6("Eddy,S (2011) Accelerated Profile HMM Searches. PLOS Computational Biol-ogy, 7(10), e1002195."),
                                             h6("Jinlong Zhang (2017). phylotools: Phylogenetic Tools for Eco-Phylogenetics. R package version 0.2.2"),
                                             h6("Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686"),
                                             h6("Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223"),
                                             h6("Mistry,J. et al. (2020) Pfam: The protein families database in 2021. Nucleic Acids Research, Vol. 49 Issue D1, D412-D419."),
                                             downloadButton("cat.download", "Download TF Catalogue")
                                           ),
                                           conditionalPanel(
                                             condition = "input.select == 'TF Annotation'",
                                             h2("Phylogenetic Trees"),
                                             h6("Please Note: Tree tips must contain the string 'CREPE' in order to work. Please visit our manual for more information."),
                                             shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE),
                                             fileInput(inputId = "metadata",
                                                       label = "Metadata",
                                                       accept = c(".csv", ".tsv") #cat.download
                                             ),
                                             h6("Metadata must reference sequence IDs with a Gene Symbol"),
                                             selectInput("dataset", "Pick Mapping to Download", c('human', 'fly', 'nearest')),
                                             tableOutput("preview"),
                                             h4("Citations for Tools and Dependencies"),
                                             h6("Eddy,S (2011) Accelerated Profile HMM Searches. PLOS Computational Biol-ogy, 7(10), e1002195."),
                                             h6("Jinlong Zhang (2017). phylotools: Phylogenetic Tools for Eco-Phylogenetics. R package version 0.2.2"),
                                             h6("Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686"),
                                             h6("Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223"),
                                             h6("Mistry,J. et al. (2020) Pfam: The protein families database in 2021. Nucleic Acids Research, Vol. 49 Issue D1, D412-D419."),
                                             downloadButton("download", "Download TF Annotation")
                                           ),
                                           conditionalPanel(
                                             condition = "input.select == 'Custom'",
                                             fileInput(inputId = "CustomFa",
                                                       label = "Protein Fasta",
                                                       accept = c(".fa", ".fasta", '.faa')
                                             ),
                                             fileInput(inputId = "file2",
                                                       label = "Tabular Domain Outfile",
                                                       accept = c(".tsv", ".csv")
                                             ),
                                             textInput(inputId = 'Custom', label = 'Custom Accension'),
                                             textInput(inputId = 'speciesCus', label = 'Species Name')
                                           )
                          ),
                          dashboardBody(
                            fluidRow(column(align = "center",
                                            "Cis-Regulatory Element-binding Protein Elucidator | CREPE",
                                            width = 4),
                                     infoBox('Visit our GitHub for the tutorial and much more! github.com/dirostri', width = 10,icon = icon("github"), fill = TRUE),
                                     tableOutput("anno"),
                                     #tableOutput("fasta"),
                                     plotOutput("distribution", width = "1000px", height="900px"),
                                     plotOutput("distributioncus", width = "1000px", height="900px")
                            )
                            #tableOutput("fasta"),
                            #tableOutput("table"),
                            #tableOutput("test")
                            #verbatimTextOutput("text"),
                            #tableOutput("meta"),
                            #tableOutput("anno")
                            #plotOutput("distribution")
                          )
)

server <- function(input, output) {
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  progress$set(message = "Making plot", value = 0)
  
  #Get the PFAM TAB Outputs
  file2 <- reactive({
    req(input$file2)
    
    ext <- tools::file_ext(input$file2$name)
    switch(ext,
           csv = readr::read_csv(input$file2$datapath),
           tsv = readr::read_tsv(input$file2$datapath),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  #Get the Fasta File
  file1 <- reactive({
    req(input$file1)
    ext <- tools::file_ext(input$file1$name)
    switch(ext,
           fasta = phylotools::read.fasta(input$file1$datapath),
           fa = phylotools::read.fasta(input$file1$datapath),
           validate("Invalid file; Please upload a .fasta or .fa file")
    )
  })
  
  #Get Custom Fa
  CustomFa <- reactive({
    req(input$CustomFa)
    ext <- tools::file_ext(input$CustomFa$name)
    switch(ext,
           fasta = phylotools::read.fasta(input$CustomFa$datapath),
           fa = phylotools::read.fasta(input$CustomFa$datapath),
           validate("Invalid file; Please upload a .fasta or .fa file")
    )
  })
  
  file_meta <- reactive({
    req(input$metadata)
    metadata = readr::read_csv(input$metadata$datapath)
  })
  
  #Get the folder of phylogenetic inferences
  volumes = getVolumes() # this makes the directory at the base of your computer.
  shinyDirChoose(input, 'folder', roots=volumes())
  path1 <- reactive({
    req(input$folder)
    return(list.files(path = parseDirPath(volumes(), input$folder))[1:3])
  })
  
  #Shows the trees to be analyed
  output$text <- renderText({
    path1()
  })
  
  #Shows the PFAM table to be used
  output$table <- renderTable({
    head(file2())[,1:4]
  })
  
  #Shows the FASTA to be used
  output$fasta <- renderTable({
    head(file1())
  })
  
  output$meta <- renderTable({
    head(file_meta())
  })
  
  #output$fasta <- renderTable({
    #head(file1())
  #})
  
  # observeEvent(input$select, {
  #   hide("distribution")
  #   hide("anno")
  # })
  
  
  #TF Distribution Output 
  observe({if(input$select == 'TF Cataloguing'){
    output$distribution <- renderPlot({
      file4()
    })
  } else if(input$select == 'Custom'){
    output$distribution <- renderPlot({
      file5()
    })
  }else {
    output$distribution <- renderText("")
  }
  })
  
  #Output annotations
  observe({if(input$select == 'TF Annotation'){
    output$anno <- renderTable({
      annotate()
    })
  } else {
    output$anno <- renderText("")
  }
  })
  
  #Output Custom
  # observe({if(input$select == 'Custom'){
  #   output$distributioncus <- renderPlot({
  #     file5()
  #   })
  # } else {
  #   output$distributioncus <- renderText("")
  # }
  # })
  ##CREPE Functions
  
  #PFAM TF definitions
  suppressMessages(list_of_DBDs <- tibble(
    Pfam_DBDs = c("TF_AP-2", "HLH", "bZIP_1", "bZIP_2", "bZIP_Maf", "GTF2I", "SART-1", "DUF573", 
                  "TCP", "B3", "CG-1", "CSD","LAG1-DNAbind", "CP2", "SRF-TF", "MBD", "NDT80_PhoG", 
                  "P53","RHD", "Runt", "STAT_bind", "TBP", "T-box", "Whirly","KilA-N", "ARID", "BrkDBD", 
                  "CENP-B_N", "CUT", "DP","Ets", "E2F_TDP", "Fork_head", "HSF_DNA-bind", "Homeobox", "IBD",
                  "IRF", "FLO_LFY", "MADF_DNA_bdg", "MAT_Alpha1", "Myb_DNA-binding", "PAX" ,"HTH_psq", 
                  "Pou", "Prox1", "Rap1-DNA-bind", "RFX_DNA_binding", "Sigma70_r2","TEA", "Vhr1", "AFT", 
                  "AP2", "AT_hook", "CBFB_NFYA","EIN3", "GCR1_C", "GRAS", "NAM", "SAND", "HMG_box","S1FA", 
                  "BAF1_ABF1", "GAGA_bind","Copper-fist", "CXC", "zf-CXXC","zf-C2HC", "zf-C2H2", "DM", 
                  "zf-Dof", "FAR1", "GATA","GCM", "mTERF", "zf-NF-X1" ,"zf-C4", "PLATZ", "SBP","MH1", 
                  "THAP" , "WRKY", "Zn_clus"),
    accession =c("PF03299", "PF00010", "PF00170", "PF07716", "PF03131", "PF02946", "PF03343", "PF04504", "PF03634", 
                 "PF02362", "PF03859", "PF00313", "PF09271", "PF04516", "PF00319", "PF01429", "PF05224", "PF00870", 
                 "PF00554", "PF00853", "PF02864", "PF00352", "PF00907", "PF08536", "PF04383", "PF01388", "PF09607",
                 "PF04218", "PF02376", "PF08781", "PF00178", "PF02319", "PF00250", "PF00447", "PF00046", "PF10416", 
                 "PF00605", "PF01698", "PF10545", "PF04769", "PF00249", "PF00292", "PF05225", "PF00157", "PF05044",
                 "PF09197", "PF02257", "PF04542", "PF01285", "PF04001", "PF08731", "PF00847", "PF02178", "PF02045", 
                 "PF04873", "PF12550", "PF03514", "PF02365", "PF01342", "PF00505", "PF04689", "PF04684", "PF06217", 
                 "PF00649", "PF03638", "PF02008", "PF01530", "PF00096", "PF00751", "PF02701", "PF03101", "PF00320", 
                 "PF03615", "PF02536", "PF01422", "PF00105", "PF04640", "PF03110", "PF03165", "PF05485", "PF03106",
                 "PF00172")
  )
  )
  
  #HMMSCAN Colnames
  hmmscan_column_names <- c("target_name", "accession", 
                            "tlen", "query_name", 
                            "accession_sequence", "qlen",
                            "E_value", "bit_score_overall", 
                            "bias_composition", "domain #", 
                            "of",   "c_Evalue",  
                            "i_Evalue",  "bit_score_domain",  
                            "bias_composition_correction",  "hmm_from",
                            "hmm_to", "ali_from",
                            "ali_to",  "env_from", 
                            "env_to", "acc",  
                            "description_of_target")
  
  custom_column_names <-c("target_name", "accession", "query_name", "env_from", "env_to")
    
  file4 <- reactive({
    req(input$file1)
    #req(input$file2) #Commented Out
    req(input$species)
    
    withProgress(message = 'TF Cataloguing', value = 0, {
      Sys.sleep(0.1)
      n <- 7
      incProgress(1/n, detail = 'Running HMMSCAN')
      fasta_tab <- phylotools::read.fasta(input$file1$datapath) %>% 
        separate(seq.name, c("ID", "info"), sep = "([ ])")
      
      #suppressMessages(tabular <- read_tsv(input$file2$datapath, col_names = hmmscan_column_names)) #Commented Out
      
      Sys.sleep(0.1)
      #incProgress(1/n, detail = 'Inputs have been loaded')
      system(paste('hmmscan --cut_nc --domtblout ',input$species,'.tsv TF_DBDs/TF_DBDs.hmm ', input$file1$datapath, sep = '')) #Added
      system(paste("sed '/^#/ d' ", input$species,".tsv | tr -s [:blank:] '\t'  >  tab_from_hmmscan", sep = '')) #Added 
      
      incProgress(1/n, detail = 'Parsing the TFs')
      suppressMessages(tabular <- read_tsv('tab_from_hmmscan', col_names = hmmscan_column_names)) #Added 
      #Processes Output for further analysis
      file.remove(paste(input$species,'.tsv', sep = ''))
      file.remove('tab_from_hmmscan')
      copy_output_temp <-  tabular %>% 
        select(query_name,target_name, env_from, env_to, accession) %>% #From tabular object
        mutate_if(is.character, str_trim) %>% #Subset these columns
        separate(accession, c("accession", "version"), sep = "([.])") %>% #Separate accession into stable acession and version
        filter(accession %in% list_of_DBDs$accession); print("Extracting TF domains") #Filter for all TF DBDs
      
      #Adds FL seqs and DBDs to Output
      copy_output_temp <- merge(x = copy_output_temp, by.x = 'query_name', y = fasta_tab, by.y = 'ID') %>% #Merges fastas for FL seqs
        select(-info) %>% #Removes the info columns
        mutate(DBD_seqs = str_sub(seq.text, start = env_from, end = env_to)) %>% #Creates the DBD_seqs 
        rename(Protein_seqs = seq.text, #Renames seq.txt for Protein seqs
               Pfam_froms = env_from, #Renames env_from for Pfam_froms
               Pfam_tos = env_to)  #Renames env_to for Pfam_tos
      
      incProgress(1/n, detail = 'Processing Inputs')
      
      #print("Preparing Final Output")
      final_copy_v2 <- merge(x = copy_output_temp, y = list_of_DBDs, by = 'accession') %>% 
        mutate(crepe_add = 'CREPE')  %>% 
        tidyr::unite(col = query_name, query_name, crepe_add, sep = '_')#
      
      #remove(copy_output_temp, tabular)
      
      suppressMessages(fams_v2 <- unique(
        final_copy_v2 %>%
          group_by(accession) %>%
          summarise(frequency=n(), 
                    Pfam_DBDs = target_name, 
                    unique = length(
                      unique(query_name)
                    )
          )
      )
      )
      
      #special rule for bzips Added on July 21. These PFs are of bzips
      bzip_pfams <- c("PF00170","PF07716","PF03131")
      temp <- final_copy_v2 %>% filter(final_copy_v2$accession %in% bzip_pfams)
      fams_v2[nrow(fams_v2)+1,'Pfam_DBDs'] <- "bZIP"
      fams_v2[nrow(fams_v2),'unique'] <- length(unique(temp$query_name))
      fams_v2 <- fams_v2 %>% filter(!accession %in% c("PF00170","PF07716","PF03131"))
      
      
      #Prepares df for output writing
      fasta_tab_CREPEdf_dom <- final_copy_v2 %>% 
        unite("fasta_head", c(query_name,Pfam_DBDs, Pfam_froms, Pfam_tos), remove = FALSE, sep = "|") %>%
        select(fasta_head, DBD_seqs)
      fasta_tab_CREPEdf_dom <- fasta_tab_CREPEdf_dom[,c('fasta_head','DBD_seqs')] %>%
        rename('seq.name'=fasta_head,'seq.text'=DBD_seqs)
      fasta_tab_CREPEdf_FL <- unique(final_copy_v2[,c('query_name','Protein_seqs')]) %>%
        rename('seq.name'=query_name,'seq.text'=Protein_seqs)
      
      incProgress(1/n, detail = 'Writing TF Catalogues')
      
      output$cat.download = downloadHandler(
        filename = 'catalagueinfo.zip',
        content = function( file){
          # Set temporary working directory
          owd <- setwd(tempdir())
          on.exit(setwd( owd))
          
          incProgress(1/n, detail = 'Writing TF Distribution Plots') 
          
          # Save Files
          #Writes the outputs
          invisible(phylotools::dat2fasta(fasta_tab_CREPEdf_FL, outfile = paste(input$species, "_CREPE_FL.faa", sep = "")))
          invisible(phylotools::dat2fasta(fasta_tab_CREPEdf_dom, outfile = paste(input$species, "_CREPE_DBD.faa", sep = "")))
          write.csv(final_copy_v2, paste(input$species, "_CREPE_df.csv", sep = ""))
          
          # Zip them up
          zip( file, c( paste(input$species, "_CREPE_FL.faa", sep = ""), 
                        paste(input$species, "_CREPE_DBD.faa", sep = ""), 
                        paste(input$species, "_CREPE_df.csv", sep = ""))
          )
          incProgress(1/n, detail = 'Finishing up!')
        })
    }
    ) 
    
    #output$distribution <- renderPlot({
    fams <- as.vector(ungroup(fams_v2) %>% 
                        select(Pfam_DBDs,unique) %>% 
                        slice_max(unique, n = 20) %>% 
                        pull(Pfam_DBDs))
    
    gg1 <- ggplot(fams_v2 %>% filter(Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      geom_bar(stat = "identity") +
      ylab("Unique Proteins") +
      #xlab("TF Families") +
      ggtitle(paste("Distribution of TF families in", input$species, sep = " ")) +
      geom_text(size = 9, position = position_stack(), aes(label = unique), vjust = -0.3) + 
      theme_classic() +
      scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
      theme(axis.text.x=element_text(angle = 45, size = 20, vjust = 1, hjust=1), axis.title.x = element_blank(), 
            axis.text.y=element_text(size = 15),
            axis.title.y = element_text(size = 15))
    
    gg2 <- ggplot(fams_v2 %>% filter(!Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      geom_bar(stat = "identity") +
      ylab("Unique Proteins") +
      xlab("TF Families") +
      geom_text(size = 9, position = position_stack(), aes(label = unique), vjust = -0.3) + 
      scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
      theme_classic() +
      theme(axis.text.x=element_text(angle = 45,size = 20, vjust = 1, hjust=1), 
            axis.text.y=element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.title.x = element_text(size = 15))
    library(cowplot)
    return( 
      # gg1 <- ggplot(fams_v2 %>% filter(Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      # geom_bar(stat = "identity") +
      # ylab("Unique Proteins") +
      # xlab("TF Families") +
      # ggtitle(paste("Distribution of TF families in", input$species, sep = " ")) +
      # geom_text(size = 4, position = position_stack(), aes(label = unique), vjust = -0.5) + 
      # theme_minimal() +
      # theme(axis.text.x=element_text(angle = 90, hjust = 0, size = 12)),
      # gg2 <- ggplot(fams_v2 %>% filter(!Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      #   geom_bar(stat = "identity") +
      #   ylab("Unique Proteins") +
      #   xlab("TF Families") +
      #   ggtitle(paste("Distribution of TF families in", input$species, sep = " ")) +
      #   geom_text(size = 4, position = position_stack(), aes(label = unique), vjust = -0.5) + 
      #   theme_minimal() +
      #   theme(axis.text.x=element_text(angle = 90, hjust = 0, size = 12))
      plot_grid(gg1, gg2,ncol = 1, align = "hv")
      
    )
    #})
  }) #Catalogue TFs + ggplot 
  
  file5 <- reactive({
    req(input$CustomFa)
    req(input$file2) #Commented Out
    req(input$speciesCus)
    
    withProgress(message = 'TF Cataloguing', value = 0, {
      Sys.sleep(0.1)
      n <- 7
      list_of_DBDs[nrow(list_of_DBDs),'accession'] <- input$Custom
      list_of_DBDs[nrow(list_of_DBDs),'Pfam_DBDs'] <- 'CustomAccession'
      #incProgress(1/n, detail = 'Running HMMSCAN')
      fasta_tab <- phylotools::read.fasta(input$CustomFa$datapath) %>% 
        separate(seq.name, c("ID", "info"), sep = "([ ])")
      
      suppressMessages(tabular <- read_tsv(input$file2$datapath, col_names = custom_column_names))
      
      Sys.sleep(0.1)
      #incProgress(1/n, detail = 'Inputs have been loaded')
      #system(paste('hmmscan --cut_nc --domtblout ',input$species,'.tsv TF_DBDs/TF_DBDs.hmm ', input$file1$datapath, sep = '')) #Added
      #system(paste("sed '/^#/ d' ", input$species,".tsv | tr -s [:blank:] '\t'  >  tab_from_hmmscan", sep = '')) #Added 
      
      incProgress(1/n, detail = 'Parsing the TFs')
      #suppressMessages(tabular <- read_tsv('tab_from_hmmscan', col_names = hmmscan_column_names)) #Added 
      #Processes Output for further analysis
      #file.remove(paste(input$speciesCus,'.tsv', sep = ''))
      #file.remove('tab_from_hmmscan')
      copy_output_temp <-  tabular %>% 
        select(query_name,target_name, env_from, env_to, accession) %>% #From tabular object
        mutate_if(is.character, str_trim) %>% #Subset these columns
        separate(accession, c("accession", "version"), sep = "([.])") %>% #Separate accession into stable acession and version
        filter(accession %in% list_of_DBDs$accession) #Filter for all TF DBDs
      
      #Adds FL seqs and DBDs to Output
      copy_output_temp <- merge(x = copy_output_temp, by.x = 'query_name', y = fasta_tab, by.y = 'ID') %>% #Merges fastas for FL seqs
        select(-info) %>% #Removes the info columns
        mutate(DBD_seqs = str_sub(seq.text, start = env_from, end = env_to)) %>% #Creates the DBD_seqs 
        rename(Protein_seqs = seq.text, #Renames seq.txt for Protein seqs
               Pfam_froms = env_from, #Renames env_from for Pfam_froms
               Pfam_tos = env_to)  #Renames env_to for Pfam_tos
      
      incProgress(1/n, detail = 'Processing Inputs')
      
      #print("Preparing Final Output")
      final_copy_v2 <- merge(x = copy_output_temp, y = list_of_DBDs, by = 'accession') %>% 
        mutate(crepe_add = 'CREPE')  %>% 
        tidyr::unite(col = query_name, query_name, crepe_add, sep = '_')#
      
      #remove(copy_output_temp, tabular)
      
      suppressMessages(fams_v2 <- unique(
        final_copy_v2 %>%
          group_by(accession) %>%
          summarise(frequency=n(), 
                    Pfam_DBDs = target_name, 
                    unique = length(
                      unique(query_name)
                    )
          )
      )
      )
      
      #special rule for bzips Added on July 21. These PFs are of bzips
      bzip_pfams <- c("PF00170","PF07716","PF03131")
      temp <- final_copy_v2 %>% filter(final_copy_v2$accession %in% bzip_pfams)
      fams_v2[nrow(fams_v2)+1,'Pfam_DBDs'] <- "bZIP"
      fams_v2[nrow(fams_v2),'unique'] <- length(unique(temp$query_name))
      fams_v2 <- fams_v2 %>% filter(!accession %in% c("PF00170","PF07716","PF03131"))
      
      
      #Prepares df for output writing
      fasta_tab_CREPEdf_dom <- final_copy_v2 %>% 
        unite("fasta_head", c(query_name,Pfam_DBDs, Pfam_froms, Pfam_tos), remove = FALSE, sep = "|") %>%
        select(fasta_head, DBD_seqs)
      fasta_tab_CREPEdf_dom <- fasta_tab_CREPEdf_dom[,c('fasta_head','DBD_seqs')] %>%
        rename('seq.name'=fasta_head,'seq.text'=DBD_seqs)
      fasta_tab_CREPEdf_FL <- unique(final_copy_v2[,c('query_name','Protein_seqs')]) %>%
        rename('seq.name'=query_name,'seq.text'=Protein_seqs)
      
      incProgress(1/n, detail = 'Writing TF Catalogues')
      
      output$cat.download = downloadHandler(
        filename = 'catalagueinfo.zip',
        content = function( file){
          # Set temporary working directory
          owd <- setwd(tempdir())
          on.exit(setwd( owd))
          
          incProgress(1/n, detail = 'Writing TF Distribution Plots') 
          
          # Save Files
          #Writes the outputs
          invisible(phylotools::dat2fasta(fasta_tab_CREPEdf_FL, outfile = paste(input$speciesCus, "_CREPE_FL.faa", sep = "")))
          invisible(phylotools::dat2fasta(fasta_tab_CREPEdf_dom, outfile = paste(input$speciesCus, "_CREPE_DBD.faa", sep = "")))
          write.csv(final_copy_v2, paste(input$speciesCus, "_CREPE_df.csv", sep = ""))
          
          # Zip them up
          zip( file, c( paste(input$speciesCus, "_CREPE_FL.faa", sep = ""), 
                        paste(input$speciesCus, "_CREPE_DBD.faa", sep = ""), 
                        paste(input$speciesCus, "_CREPE_df.csv", sep = ""))
          )
          incProgress(1/n, detail = 'Finishing up!')
        })
    }
    ) 
    
    #output$distribution <- renderPlot({
    fams <- as.vector(ungroup(fams_v2) %>% 
                        select(Pfam_DBDs,unique) %>% 
                        slice_max(unique, n = 20) %>% 
                        pull(Pfam_DBDs))
    
    gg1 <- ggplot(fams_v2 %>% filter(Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      geom_bar(stat = "identity") +
      ylab("Unique Proteins") +
      #xlab("TF Families") +
      ggtitle(paste("Distribution of TF families in", input$speciesCus, sep = " ")) +
      geom_text(size = 9, position = position_stack(), aes(label = unique), vjust = -0.3) + 
      theme_classic() +
      scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
      theme(axis.text.x=element_text(angle = 45, size = 20, vjust = 1, hjust=1), axis.title.x = element_blank(), 
            axis.text.y=element_text(size = 15),
            axis.title.y = element_text(size = 15))
    
    gg2 <- ggplot(fams_v2 %>% filter(!Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      geom_bar(stat = "identity") +
      ylab("Unique Proteins") +
      xlab("TF Families") +
      geom_text(size = 9, position = position_stack(), aes(label = unique), vjust = -0.3) + 
      scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
      theme_classic() +
      theme(axis.text.x=element_text(angle = 45,size = 20, vjust = 1, hjust=1), 
            axis.text.y=element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.title.x = element_text(size = 15))
    library(cowplot)
    return( 
      # gg1 <- ggplot(fams_v2 %>% filter(Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      # geom_bar(stat = "identity") +
      # ylab("Unique Proteins") +
      # xlab("TF Families") +
      # ggtitle(paste("Distribution of TF families in", input$species, sep = " ")) +
      # geom_text(size = 4, position = position_stack(), aes(label = unique), vjust = -0.5) + 
      # theme_minimal() +
      # theme(axis.text.x=element_text(angle = 90, hjust = 0, size = 12)),
      # gg2 <- ggplot(fams_v2 %>% filter(!Pfam_DBDs %in% fams), aes(reorder(Pfam_DBDs, -unique), unique)) +
      #   geom_bar(stat = "identity") +
      #   ylab("Unique Proteins") +
      #   xlab("TF Families") +
      #   ggtitle(paste("Distribution of TF families in", input$species, sep = " ")) +
      #   geom_text(size = 4, position = position_stack(), aes(label = unique), vjust = -0.5) + 
      #   theme_minimal() +
      #   theme(axis.text.x=element_text(angle = 90, hjust = 0, size = 12))
      plot_grid(gg1, gg2,ncol = 1, align = "hv")
      
    )
    #})
  })
  annotate <- reactive({
    req(input$folder)
    req(input$metadata$datapath)
    
    withProgress(message = 'Begin TF Annotation', value = 0, {
      
      Sys.sleep(0.1)
      n <- 4
      
      #folder <- list.files(path = parseDirPath(volumes(), input$folder), full.names = TRUE)
      folder <- list.files(path = '~/CREPE_shiny/MinimalReproducibleExample/TF_Annotation/', full.names = T)
      #folder <- list.files('~/OrthoFinderR/Results_test/', full.names = T)
      #print(folder)
      
      #metadata <- read_csv(file = input$metadata$datapath) %>% 
      #select("ensembl_gene_id", "external_gene_name") %>% distinct()
      metadata <- read_csv(file = '~/Documents/CREPE_analysis/metadata_updated.csv') %>% 
        select("ensembl_gene_id", "external_gene_name") %>% distinct()
      
      #metadata <- read_csv(file = "~/Documents/CREPE_analysis/metadata_updated.csv") %>% 
      #select("ensembl_gene_id", "external_gene_name") %>% distinct()
      
      head(metadata)
      
      results_L <- data.frame(Distance=NULL, Most_Related=NULL)
      results_F <- data.frame(Distance=NULL, Most_Related=NULL)
      results_D <- data.frame(Distance=NULL, Most_Related=NULL)
      
      incProgress(1/n, detail = 'Commencing Tree Parsing')
      for (tree in folder){
        l <- NULL
        f <- NULL
        d <- NULL
        print(tree)
        if (file.exists(tree) == TRUE){
          gene_tree <- phytools::read.newick(tree)
          #phytools::plotTree(gene_tree)
          #gene_names <- gsub("_gene_ensembl_TFpfam", "", gene_tree$tip.label)
          gene_names <- gene_tree$tip.label
          #gene_names <- gsub("drosophila_melanogaster", "dmelanogaster", gene_names)
          #gene_names <- gsub("_eg_gene_TFpfam", "", gene_names)
          #gene_names <- gsub("20200818_bombyx_moriCREPEdf_FL_", "bombyxCREPE_", gene_names)
          #gene_names <- gsub("20200818_hel_mel_melCREPEdf_FL_", "CREPEHelMel_", gene_names)
          #gene_names <- gsub("Bicyclus_anynana_BaGv2_CREPE_FL_", "CREPEBanyana_", gene_names)
          #gene_names <- gsub("bmori_SILKdb", "SilkdbBmori", gene_names)
          #gene_names <- gsub("dplexippus_eg_gene_TFpfam_", "Danausp_", gene_names)
          #gene_names <- gsub("XP_", "XP", gene_names)
          #gene_names <- gsub("_NP_", "_NP", gene_names)
          #gene_names <- sub("_", "", gene_names)
          #gene_names <- sub("_(.*)pep_all", "", gene_names)
          #gene_names <- sub("\\..*", "", gene_names)
          #gene_names <- sub("SUM149CREPE", "CREPESUM149", gene_names)
          gene_names <- data.frame(names=gene_names) %>% separate(names, into = c("species", "geneids"), sep = "_(?=[^_]+$)", extra = 'merge') %>% mutate(index=1:length(gene_names), ensembl_gene_id = geneids)
          gene_names$geneids <- sub("\\..*", "", gene_names$geneids)
          
          
          
          gene_names2<- merge(x = gene_names, y = metadata, by.x = "geneids",  all.x = T) %>%
            arrange(index)
          #gene_names2 <- inner_join(gene_names,metadata)
          gene_names2$external_gene_name <- replace_na(gene_names2$external_gene_name,"noAnnotationFound")
          gene_names2 <-gene_names2 %>% unite(species,ensembl_gene_id, external_gene_name, col = "forTree", sep = "_")
          gene_names2 <- as.data.frame(apply(gene_names2,2, str_remove_all, " "))
          
          gene_tree[["tip.label"]] <- gene_names2$forTree
          
          homo_ids <- as.character(gene_tree$tip.label[grepl(c("Homo_sapiens"), x = gene_tree$tip.label)])
          crepe_ids <- as.character(gene_tree$tip.label[grepl(c("CREPE"), x = gene_tree$tip.label)])
          dmel_ids <- as.character(gene_tree$tip.label[grepl(c("Drosophila_melanogaster"), x = gene_tree$tip.label)])
          
          f <- as.data.frame(ape::cophenetic.phylo(gene_tree)) 
          l <- f
          d <- f 
          
          f <- f %>% select(all_of(homo_ids), all_of(crepe_ids)) %>% 
            rownames_to_column(var = "species") %>% 
            dplyr::filter(species %in% c(homo_ids, crepe_ids)) %>% 
            select(species,contains("CREPE")) %>%
            filter(!grepl("_noAnnotationFound", species)) %>%
            column_to_rownames(var = "species")
          
          d <- d %>% select(all_of(dmel_ids), all_of(crepe_ids)) %>% 
            rownames_to_column(var = "species") %>% 
            filter(species %in% c(dmel_ids, crepe_ids)) %>% 
            select(species,contains("CREPE")) %>%
            filter(!grepl("_noAnnotationFound", species)) %>%
            column_to_rownames(var = "species")
          
          l <- l %>% select(contains("CREPE")) %>% 
            rownames_to_column(var="species") %>%
            filter(!grepl("noAnnotationFound", species)) %>%
            column_to_rownames(var = "species")
          
          if(nrow(l) >=1){
            results_l <- data.frame(Distance=apply(l[,,drop=F],2,min),
                                    Most_Related=(rownames(x = l)[apply(l[,,drop=F],2,which.min)]))
            results_L <- rbind(results_L, results_l)
            print(head(results_l))
          }
          
          if(nrow(f) >=1){
            results_f <- data.frame(Distance=apply(f[,,drop=F],2,min),
                                    Most_Related=(rownames(x = f)[apply(f[,,drop=F],2,which.min)]))
            results_F <- rbind(results_F, results_f)
            print(head(results_f))
          }
          if(nrow(d) >=1){
            results_d <- data.frame(Distance=apply(d[,,drop=F],2,min),
                                    Most_Related=(rownames(x = d)[apply(d[,,drop=F],2,which.min)]))
            results_D <- rbind(results_D, results_d)
            print(head(results_d))
          }
          
        }else{
          print(paste(tree, sep = "", " does not exist"))
        }
      }
      print('For Loop ended')
      #print(results_L)
      results_L <- results_L %>% 
        separate(Most_Related, sep = "_(?=[^_]+$)", into = c("extra", "Description")) %>% 
        separate(extra, sep = "_(?=[^_]+$)", into = c("MostRelatedSpecies", "ProteinAcc")) %>% rownames_to_column()
      #print('Results L ended')
      
      #print(results_F)
      results_F <- results_F %>%
        separate(Most_Related, sep = "_(?=[^_]+$)", into = c("extra", "Description")) %>% 
        separate(extra, sep = "_(?=[^_]+$)", into = c("MostRelatedSpecies", "ProteinAcc")) %>% rownames_to_column()
      #print('Results F ended')
      
      #print(results_D)
      results_D <- results_D %>%
        separate(Most_Related, sep = "_(?=[^_]+$)", into = c("extra", "Description")) %>% 
        separate(extra, sep = "_(?=[^_]+$)", into = c("MostRelatedSpecies", "ProteinAcc")) %>% rownames_to_column()
      #print('Results D ended')
      results_F$rowname <- gsub("_noAnnotationFound","", results_F$rowname)
      results_L$rowname <- gsub("_noAnnotationFound","", results_L$rowname)
      results_D$rowname <- gsub("_noAnnotationFound","", results_D$rowname)
      
      results_F <- as.data.frame(results_F %>% separate(rowname, sep = "_", into = c("species", "species_prot_acc"), extra = "merge"))
      results_L <- as.data.frame(results_L %>% separate(rowname, sep = "_", into = c("species", "species_prot_acc"), extra = "merge"))
      results_D <- as.data.frame(results_D %>% separate(rowname, sep = "_", into = c("species", "species_prot_acc"), extra = "merge"))
      
      output$anno <- renderTable({
        results_F
      })
      
      incProgress(1/n, detail = 'Preparing outputs for download')
      output$download <- downloadHandler(
        filename = function() {
          paste(input$dataset, 'tsv', sep = '.')
        },
        content = function(file){
          if(input$dataset == 'human'){
            readr::write_tsv(results_F, file = file)
          } else if (input$dataset == 'fly') {
            readr::write_tsv(results_D, file = file)
          }else if (input$dataset == 'nearest'){
            readr::write_tsv(results_L, file = file)
          }
        }
      )
      
      output$download <- downloadHandler(
        filename = function() {
          paste(input$dataset, 'tsv', sep = '.')
        },
        content = function(file){
          if(input$dataset == 'human'){
            readr::write_tsv(results_F, file = file)
          } else if (input$dataset == 'fly') {
            readr::write_tsv(results_D, file = file)
          }else if (input$dataset == 'nearest'){
            readr::write_tsv(results_L, file = file)
          }
        }
      )
      
      incProgress(1/n, detail = 'Finished!')
      
    })
    
  })#
  
  
}
# Run the application 
shinyApp(ui = ui, server = server)

