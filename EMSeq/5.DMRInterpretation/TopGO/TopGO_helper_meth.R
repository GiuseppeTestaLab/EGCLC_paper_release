# AUTHORS: Gaja Matassa and Amaia Tintori

checkGeneVectors <- function(GeneVec, stop=TRUE){
  #Print the number of levels
  message(paste('Gene vector contains levels:',
                paste(levels(GeneVec), collapse = ',')))
  #Specific checks
  if (stop == TRUE) {
    #Check that there are significant terms
    if (length(levels(GeneVec)) == 1){
      stop('There are no DEGs to test for!')
      
      #Check that there are only two levels (T, F)
    } else if (length(levels(GeneVec)) > 2) {
      stop('Only two levels expected')
    }
  }
}

topGOGeneVectors_meth <- function(Annotated_DMR, Cond1, Cond2, genomic_type = "Promoter", gene_type = "protein_coding", Universe){
  
  comparison <- paste0(Cond1, "vs", Cond2)
  DMRs_df <- data.frame(Annotated_DMR[[comparison]])
  Genes_df <- data.frame(Universe)
  colnames(Genes_df) <- "GeneUniv"
  
  if (genomic_type == "all" && gene_type == "all") {
    break
    } 
  else if (genomic_type == "all") {
    DMRs_df <- DMRs_df %>% filter(gene_biotype == gene_type)
    }
  else if (gene_type == "all"){
    DMRs_df <- DMRs_df[grep(genomic_type, DMRs_df$ChipSeekerAnn),]
    }
  else{
    DMRs_df <- DMRs_df[grep(genomic_type, DMRs_df$ChipSeekerAnn),] %>% filter(gene_biotype == gene_type)
    }
  
  if (!all(unique(na.omit(DMRs_df$hgnc_symbol)) %in%  Genes_df$GeneUniv)){
    print("There are some genes associated to DMRs, that are not in the gene universe! They will not be considered...")} 
  else{ print("All the genes associated to DMRs are in the gene universe!") }
  
  Genes_df <- Genes_df %>% mutate(
    DMGenes = case_when(GeneUniv %in% unique(na.omit(DMRs_df$hgnc_symbol)) ~ 1, TRUE ~ 0),
    DMGenesUp = case_when(GeneUniv %in% unique(na.omit(DMRs_df[DMRs_df$diff.Methy<0, "hgnc_symbol"])) ~ 1,
                          TRUE ~ 0),
    DMGenesDown = case_when(GeneUniv %in% unique(na.omit(DMRs_df[DMRs_df$diff.Methy>0, "hgnc_symbol"])) ~ 1,
                            TRUE ~ 0)) %>%
    mutate_at(c('DMGenes', 'DMGenesUp', 'DMGenesDown'), factor) #%>% #set type
  #distinct(.data$GeneUniv, .keep_all = TRUE)
  
  GeneVectors <- Genes_df[,c('DMGenes', 'DMGenesDown', 'DMGenesUp')] %>% as.list()
  
  for (x in names(GeneVectors)){
    names(GeneVectors[[x]]) <- Genes_df$GeneUniv}
  
  return(GeneVectors)
}

topGOResults_new <- function(Genes, gene2GO, desc=NULL, ontology='BP',
                         nodeSize=8, algorithm='weight01', statistic='fisher',
                         EnTh=1.5, PvalTh=0.01, geneTh=4, minTerms=15, maxAnn = 500,
                         saveRes=FALSE, outDir=getwd(), fileName='TopGO', ...) {
  
  #Check that the gene vectors contain significant DEGs
  checkGeneVectors(Genes, stop=TRUE)
  TopGORes <- list()
  
  # 1. Data Preparation
  # Set description
  if (is.null(desc)) {
    desc <- paste(ontology, deparse(substitute(Genes)))
  }
  
  # Create TopGo object
  TopGORes$GOdata <- methods::new('topGOdata',  #object class
                                  description=desc, ontology=ontology,
                                  allGenes=Genes, gene2GO=gene2GO,
                                  nodeSize=nodeSize,
                                  annotat=topGO::annFUN.gene2GO, ...)
  
  # 2. Statistical testing
  TopGORes$Test <- topGO::runTest(TopGORes$GOdata,
                                  algorithm=algorithm, statistic=statistic)
  
  # 3. Generation of result table
  TopGORes$ResAll <- topGO::GenTable(TopGORes$GOdata, Statistics=TopGORes$Test,
                                     topNodes=length(TopGORes$Test@score))
  
  # Save complete results
  if(saveRes == TRUE){  #Write table
    utils::write.table(TopGORes$ResAll, file=paste0(outDir, fileName, '.txt'),
                       sep='\t', quote=TRUE, row.names=FALSE)
  }
  
  # 4. Selection based on enrichment threshold
  ESel <- TopGORes$ResAll %>%
    mutate(ER = round(.data$Significant/.data$Expected, 2)) %>%
    filter(.data$ER > EnTh) %>%
    filter(.data$Significant >= geneTh) %>%
    filter(.data$Annotated < maxAnn) #added by Gaja
  #? filter(if(!is.null(geneTh)) Significant >= geneTh else TRUE)
  
  # Number of significant genes
  nSig <- ESel %>%
    filter(as.numeric(.data$Statistics) <= PvalTh) %>%
    nrow()
  
  TopGORes$ResSel <- ESel[seq(max(nSig, minTerms)),] %>% #max among the two
    stats::na.omit() #Remove NAs created if nrow(ESel) < minTerms
  #n <- max(dim(dplyr::filter(ESel, as.numeric(Statistics) <= PvalTh))[1], minTerms)
  
  return(TopGORes)
}

topGOBarplot_new <- function(TopGORes, terms=15, pvalTh=0.01, plotTitle=NULL,
                             palette=NULL, flip_x=FALSE) {
  
  # 1. Definition of the color palette
  if (is.null(palette)) {
    palette <- grDevices::colorRampPalette(c('gold', 'forestgreen'))(terms)
  }
  
  # 2. Title definition
  if (is.null(plotTitle)) {
    plotTitle = deparse(substitute(TopGORes))
  }
  
  # Add full GO term
  GO <- as.list(GO.db::GOTERM)
  
  # add the complete GO term name
  #TopGORes['extTerm'] <- stringr::str_wrap(sapply(as.character(TopGORes[['GO.ID']]),
  #                                             function(x) X = GO[[x]]@Term), width=50)
  TopGORes['extTerm'] <- stringr::str_wrap(sapply(as.character(TopGORes[['GO.ID']]),
                                                  function(x) X = GO[[x]]@Term))
  
  # 3. Dataframe reorder
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGORes$Statistics <- ifelse(grepl('<', TopGORes$Statistics), 1e-30, TopGORes$Statistics)
  # then I order the results
  ResOrdered <- TopGORes %>%
    dplyr::arrange(desc(as.numeric(Statistics))) %>% #This sort the dataframe but NOT the factor levels
    dplyr::mutate(GO.ID=factor(GO.ID, levels=GO.ID)) %>%   # Updates the factor levels for ggplot2
    dplyr::slice_tail(n=terms)
  
  
  # 4. x-axis limit definition
  MaxVal <- round(max(-log10(as.numeric(ResOrdered$Statistics))), 0) + 1
  
  # 4. BarPlot
  if (flip_x == TRUE){ #Bars going left
    GO_BP <- ggplot(data=ResOrdered, aes(x=.data$GO.ID, y=log10(as.numeric(.data$Statistics)), fill=.data$GO.ID)) +
      geom_bar(stat='identity', aes(alpha=0.75)) +
      geom_text(aes(y=0), label=ResOrdered$extTerm, hjust=1.01) +
      scale_y_continuous(breaks=seq(-MaxVal,0,2), labels=abs(seq(-MaxVal,0,2)),
                         limits=c(-MaxVal,0), expand=c(0.025, 0.025)) +
      geom_hline(yintercept=log10(pvalTh), col='darkred', lty='longdash')
    
  } else {  #Bars going right
    GO_BP <- ggplot(data=ResOrdered, aes(x=.data$GO.ID, y=-log10(as.numeric(.data$Statistics)), fill=.data$GO.ID)) +
      geom_bar(stat='identity', aes(alpha=0.75)) +
      geom_text(aes(y=0), label=ResOrdered$extTerm, hjust=-0.01) +
      scale_y_continuous(breaks=seq(0,MaxVal,2), labels=abs(seq(0, MaxVal, 2)),
                         limits=c(0,MaxVal), expand=c(0.025, 0.025)) +
      geom_hline(yintercept=-log10(pvalTh), col='darkred', lty='longdash')
  }
  #Common parameters
  GO_BP <- GO_BP +
    coord_flip() +
    scale_fill_manual(values=palette) +
    labs(x='', y='-log10 PValue', title=plotTitle) +
    theme_bw() +
    theme(legend.position='none',
          axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12),
          plot.title= element_text(face='bold', colour='darkred', size=12))
  
  return(GO_BP)
}   

topGOBarplotAll_new <- function(TopGOResAll, TopGOResDown, TopGOResUp,
                            terms=15, pvalTh=0.01, plotTitle=NULL, cols=NULL) {
  
  ## 1. Graphical stuff
  if (is.null(plotTitle)) { #Title
    plotTitle = 'TopGO Analysis'
  }
  if (is.null(cols) | is.null(names(cols))){ #Accent colors
    cols <- c(All='forestgreen', Down='blue', Up='red')
    warning("If you did specify the 'cols' argument please make sure to corectly
                set names. See `?topGOBarplotAll`" )
  } else {
    warning("Please make sure you specified names of the cols vector correctly.
                See `?topGOBarplotAll`" )
  }
  #Color palettes
  pals <- lapply(cols, {function(el) grDevices::colorRampPalette(c('gold', el))(terms)})
  
  ## 2. Create specified barplots
  BPs <- list()
  
  if (!missing(TopGOResAll)){ # Barplot for all genes
    BarAll <- topGOBarplot_new(TopGOResAll, terms=terms, pvalTh=pvalTh,
                           plotTitle=paste(plotTitle, ' All'),
                           palette=pals$All)
    BPs[[length(BPs) + 1]] <- BarAll
  }
  if (!missing(TopGOResDown)) { # Barplot for down genes
    BarDown <- topGOBarplot_new(TopGOResDown, terms=terms, pvalTh=pvalTh, flip_x=TRUE,
                            plotTitle=paste(plotTitle, ' Down'),
                            palette=pals$Down)
    BPs[[length(BPs) + 1]] <- BarDown
  }
  if (!missing(TopGOResUp)) { # Barplot for up genes
    BarUp <- topGOBarplot_new(TopGOResUp, terms=terms, pvalTh=pvalTh,
                          plotTitle=paste(plotTitle, ' Up'),
                          palette=pals$Up)
    BPs[[length(BPs) + 1]] <- BarUp
  }
  
  # 3. Plot
  grid.arrange(grobs=BPs, ncol=length(BPs))
}

bubbleplot <- function(Res, terms = 15, Ont = NULL, pvalTh=0.01, plotTitle = NULL, fillCol = "forestgreen") {
  
  # DataFrame containing the previously selected enriched GO terms, the significant genes referring to them, the ER and the PValue.
  GO <- Res[, c("GO.ID", "Term", "Annotated", "Significant", "ER", "Statistics")]
  colnames(GO) <- c("GO.ID", "Term", "Annotated", "Significant", "ER", "PValue") #Changing of dataframe columns names
  
  # List objects and their structure contained in the dataframe 'GO'
  #ls.str(GO)
  
  # Transform the columna 'GeneRatio', 'ER' and 'PValue' into a numeric variable
  GO$GeneRatio <- as.numeric(GO$Significant/GO$Annotated)
  GO$PValue <- as.numeric(GO$PValue)
  GO$ER <- as.numeric(GO$ER)
  
  GeneOnt <- as.list(GO.db::GOTERM)
  
  # add the complete GO term name
  GO$extTerm <- stringr::str_wrap(sapply(as.character(GO[['GO.ID']]),
                                         function(x) X = GeneOnt[[x]]@Term))
  GO$GO_Term <- paste(GO$extTerm, paste0("(", GO$GO.ID, ")")) #Add IDs to terms
  
  
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  GO$PValue <- ifelse(grepl('<', GO$PValue), 1e-30, GO$PValue)
  GO<- GO %>%
    dplyr::arrange(desc(as.numeric(PValue))) %>% #This sort the dataframe but NOT the factor levels
    dplyr::slice_tail(n=terms)
  
  # Transform PValue values by -log10('PValue values')
  GO$'|log10(PValue)|' <- -(log10(GO$PValue))
  
  if (is.null(plotTitle)) { #Title
    plotTitle = "TopGO results"
  }
  
  if (Ont == "BP"){x_label <- "GO Biological Processes"}
  else  if (Ont == "MF"){x_label <- "GO Molecular Functions"}
  else if (Ont == "CC"){x_label <- "GO Cellular Components"}
  else {x_label <- "GO Terms"}
  
  ggplot(GO, aes(x = GO_Term, y = `|log10(PValue)|`)) +
    geom_hline(yintercept = -log10(pvalTh), linetype="dashed", 
               color = "darkgreen", size=.5)+
    geom_point(data=GO,aes(x=GO_Term, y=`|log10(PValue)|`,size = GeneRatio, colour = ER), alpha=.7)+
    scale_x_discrete(limits= GO$GO_Term)+
    scale_color_gradient(low="grey90",high=fillCol,limits=c(0, NA))+
    coord_flip()+
    theme_bw()+
    theme(axis.ticks.length=unit(-0.1, "cm"),
          axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
          axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
          axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5)+
    xlab(x_label)+
    ylab("-log(PValue)")+
    labs(color="ER", size="Gene ratio:\nSign/Ann", title = plotTitle) + #Replace by your variable names; \n allow a new line for text
    guides(shape = guide_legend(order=2),colour = guide_colourbar(order=1)) 
  
}

plotGenesInTerm_v2 <- function(TopGOResults, GOdata, DMRs, nterms=12, ngenes=8,
                               plotTitle=NULL, Interactive=FALSE, fillCol='forestgreen'){
  
  #Prepare the DF for plotting
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGOResults$Statistics <- ifelse(grepl('<', TopGOResults$Statistics), 1e-30,
                                    TopGOResults$Statistics)
  # then I order the results
  ResOrdered <- transform(TopGOResults,
                          GO.ID=reorder(GO.ID,
                                        -as.numeric(Statistics)))[1:nterms,]
  #Create the data frame
  finalDF <- data.frame()
  
  for (i in seq(nterms)){
    if (!is.na(ResOrdered$GO.ID[i])) {
      
      #Create the basic dataframe
      GOid <- as.character(ResOrdered$GO.ID[i])
      GOterm <- as.character(ResOrdered$Term[i])
      Genes <- topGO::genesInTerm(GOdata, GOid)
      Scores <- topGO::scoresInTerm(GOdata, GOid)
      
      GOdf <- as.data.frame(cbind(Gene=unlist(Genes), Score=unlist(Scores)),
                            row.names = FALSE) %>%
        filter(Score == 2) %>%  #select genes that are present
        mutate(
          GOid = case_when(TRUE ~ GOid),   #repeat GOid
          GOterm = case_when(TRUE ~ GOterm)) %>%
        mutate(across(c(GOterm, Gene), as.factor)) %>%  #change types
        mutate(across(c(Score), as.integer))
      
      
      #Add Score for DMG significance
      #The same gene could be assigned to more regions (that could have different areaStat), 
      #the most significant areaStat is prioritized. COnsidering that they are ordered by areaStat I just remove the duplicates.
      GOdf['areaStat'] <-  DMRs[DMRs$hgnc_symbol %in% GOdf$Gene, ] %>% filter(!duplicated(hgnc_symbol)) %>% select(areaStat)
      
      #Consider only tha absolute value of the areaStat
      GOdf['areaStat'] <- log(abs(GOdf['areaStat']))
      
      #Sort and add rank for x axis in plot
      GOdf <- GOdf %>%
        slice_head(n = ngenes)
      
      #Bind to dataframe
      finalDF <- as.data.frame(rbind(finalDF, GOdf))
      
    } else{next}
  }
  
  ## Plot
  #Graphical stuff
  if (is.null(plotTitle)) { #Title
    plotTitle = "Genes in Term"
  }
  sub <- paste("Top", nterms, "GO terms | Top", ngenes, "genes by areaStat")
  #palette <- colorRampPalette(c('#dfe5ef', '#0FBBE6', '#0635EE'))(nterms)
  
  BP <- ggplot(finalDF, aes(y=GOterm, x=.data$Score, fill=.data$areaStat, label=.data$Gene, label2=GOid)) +
    geom_col(col = 'white', alpha = 0.6, width = 1) +
    #viridis::scale_fill_viridis(begin = 0.3) +
    scale_fill_gradient(high = fillCol, low = "grey90", name = "logAreaStat") + #mid = "white", #midpoint = .02
    scale_y_discrete(limits = rev(levels(finalDF$GOterm))) + #otherwise most significant below
    labs(x='Genes in Term', y='GO Terms',
         title=plotTitle, subtitle=sub) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
          plot.subtitle = element_text(size=12, hjust=0.5), #colour='darkred',
          text = element_text(size = 15),
          axis.text = element_text(size = 14),
          axis.text.x = element_blank(), #axis.text.y = element_blank(),
          axis.ticks.x =element_blank(), axis.ticks.y = element_blank())
  
  #4. Return interactive or static plot
  if (Interactive == TRUE){  #With ggplotly
    #ggplotly doesn't show subtitles: we use this workaround
    BP <- plotly::ggplotly(BP) %>%
      plotly::layout(title= list(text=paste0(plotTitle, '<br>', '<sup>', sub),
                                 tooltip=c("Gene", "GOid", "areaStat"))) # selects the aesthetics to include in the tooltip
    return(BP)
    
  } else {  #Static plot
    BP <- BP + geom_text(size = 3.5, position = position_stack(vjust = 0.5), )
  }
  
  return(BP)
  
}

TableGenesInTerm <- function(TopGOResults, GOdata, nterms=12){
  
  #Prepare the DF for plotting
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGOResults$Statistics <- ifelse(grepl('<', TopGOResults$Statistics), 1e-30,
                                    TopGOResults$Statistics)
  # then I order the results
  ResOrdered <- transform(TopGOResults,
                          GO.ID=reorder(GO.ID,
                                        -as.numeric(Statistics)))[1:nterms,]
  GeneOnt <- as.list(GO.db::GOTERM)
  
  #Create the data frame
  finalDF <- data.frame()
  
  for (i in seq(nterms)){
    if (!is.na(ResOrdered$GO.ID[i])) {
      
      #Create the basic dataframe
      GOid <- as.character(ResOrdered$GO.ID[i])
      
      # add the complete GO term name
      GOextTerm <- stringr::str_wrap(sapply(as.character(ResOrdered$GO.ID[i]),
                                             function(x) X = GeneOnt[[x]]@Term))
      GOterm <- as.character(paste(GOextTerm, paste0("(", GOid, ")"))) #Add IDs to terms
      
      Genes <- topGO::genesInTerm(GOdata, GOid)
      Scores <- topGO::scoresInTerm(GOdata, GOid)
      
      GOdf <- as.data.frame(cbind(Gene=unlist(Genes), Score=unlist(Scores)),
                            row.names = FALSE) %>%
        filter(Score == 2) %>%  #select genes that are present
        mutate(
          GOid = case_when(TRUE ~ GOid),   #repeat GOid
          GOterm = case_when(TRUE ~ GOterm)) %>%
        mutate(across(c(GOterm, Gene), as.factor)) %>%  #change types
        mutate(across(c(Score), as.integer))
    
      #Bind to dataframe
      finalDF <- as.data.frame(rbind(finalDF, GOdf))
    } else{next}
    }
  return(finalDF)
  }

topGOGeneVectors_meth_v2 <- function(gene_symbols, genomic_type = "Promoter", gene_type = "protein_coding", Universe){
  
  Genes_df <- data.frame(Universe)
  colnames(Genes_df) <- "GeneUniv"
  
  if (genomic_type == "all" && gene_type == "all") {
    print("All genomic regions and gene types will be kept")
  } 
  else if (genomic_type == "all") {
    DMRs_df <- DMRs_df %>% filter(gene_biotype == gene_type)
  }
  else if (gene_type == "all"){
    DMRs_df <- DMRs_df[grep(genomic_type, DMRs_df$ChipSeekerAnn),]
  }
  else{
    DMRs_df <- DMRs_df[grep(genomic_type, DMRs_df$ChipSeekerAnn),] %>% filter(gene_biotype == gene_type)
  }
  
  if (!all(unique(na.omit(gene_symbols)) %in%  Genes_df$GeneUniv)){
    print("There are some genes symbols which are not in the gene universe! They will not be considered...")} else{print("All the genes symbols are in the gene universe!") }
  
  Genes_df <- Genes_df %>% mutate(
    Genes_tested = case_when(GeneUniv %in% unique(na.omit(gene_symbols)) ~ 1, TRUE ~ 0)) %>% mutate_at(c('Genes_tested'), factor)
  #DMGenesUp = case_when(GeneUniv %in% unique(na.omit(DMRs_df[DMRs_df$diff.Methy<0, "hgnc_symbol"])) ~ 1, TRUE ~ 0),
  #DMGenesDown = case_when(GeneUniv %in% unique(na.omit(DMRs_df[DMRs_df$diff.Methy>0, "hgnc_symbol"])) ~ 1, TRUE ~ 0)) 
  #%>% #set type
  #distinct(.data$GeneUniv, .keep_all = TRUE)
  
  GeneVectors <- Genes_df[,c("Genes_tested")] %>% as.list() %>% setNames(Genes_df$GeneUniv)
}