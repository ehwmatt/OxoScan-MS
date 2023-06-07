### Functions script for OxoScan-MS analysis


###########################################################################################
###################### search fragment functions ##########################################
###########################################################################################
search.fragments <- function(index, tolerance = 0.1){ 
  
  temp <- suppressMessages(
    cid.txt[which.min(abs(cid.txt$`Mass/Charge`- index)),1:2] %>%
      subset(select=which(!duplicated(names(.)))) %>% # remove duplicate columns
      rename(CID.Peak = "Mass/Charge",
             CID.Area = "Area") %>%
      mutate(V1 = index) %>%
      left_join(hcd.txt) %>%
      rename(HCD.Peak = "V1",
             HCD.Area = "V2") %>%
      mutate(PPM.Error = 10^6 *abs(CID.Peak - HCD.Peak) / HCD.Peak,
             Abs.Error = abs(CID.Peak - HCD.Peak)) %>%
      filter(Abs.Error < tolerance))
  
  return(temp)

}


bind.fragments <- function(ion){
  temp <- rbind(cid.hcd,
                search.fragments(ion))
}

match.fragments <- function(){
  
  
  cid.hcd <-  lapply(X = hcd.txt$V1,
                     # search - returns a row from CID 
                     FUN = search.fragments()
  )
  
}



######################################################################
#####################   DDA Peak Searching   #########################
######################################################################

#### 
dda.peak.extract <- function(peak_number,
                             tolerance = 0.5,
                             Score.Filter = 150,
                             LogProb.Filter = 3){
  
  temp.mass <- filter(dia.precursors, peak_num == peak_number) %>%
    pull(`[M+H]+`)
  
  temp.table <- byonic %>%
    filter(`Calc M+H` > temp.mass - tolerance,
           `Calc M+H` < temp.mass + tolerance,
           Score > Score.Filter,
           `|Log Prob|` > LogProb.Filter) %>%
    select(Peptide, Glycan, `M+H`, `Observed m/z`, `z`, `Calc M+H`, StartPosition, ## Structure info
           `Off-by-x error`, Score, Delta, `Delta Mod`, `|Log Prob|`, ## Assignment quality metrics
           `Protein Name`, `Scan #`, `Scan Time`, Pool, Glycan.Short) %>% ## Other info
    mutate(Scan = sub(`Scan #`, 
                      pattern = ".*scan=(\\d*)",
                      replacement = "\\1"),
           Genes = ifelse(str_detect(`Protein Name`, "TRFE") == TRUE, ### TF labelled as common contaminant in 
                          yes = "TF", 
                          no = sub(`Protein Name`, 
                                   pattern = ".*GN=(.*) PE.*",
                                   replacement = "\\1")),
           peak_num = peak_number) %>%
    select(-`Scan #`)
  
  return(temp.table)
  
}



##########################################################################
######### Boxplot plotting functions for protein/glycopeptides ###########
##########################################################################



### Plot protein abundances from protein.normalised dataset

protein.boxplot <- function(gene.name, 
                            plot.title,
                            bottom.plot = F,
                            strip.size = 10.5,
                            y.axis.breaks = 0.5){
  
  temp.prot.plot <- protein.normalised %>%
    filter(Genes == gene.name,
           File.Name %in% prot.df$File.Name) %>%
    select(Genes, I.WHO, F.Severity, Protein.Intensity) %>%
    distinct() %>%
    ggplot(aes(x = factor(F.Severity, levels=c("Healthy","Mild","Moderate", "Severe")), 
               y = log2(Protein.Intensity),
               fill = F.Severity)) + 
    geom_boxplot(inherit.aes = T, fill = NA, alpha=0.75, outlier.shape = NA) + 
    geom_jitter(aes(colour = F.Severity), 
                position = position_jitter(width = 0.1, height = 0.1),
                alpha = 0.75, 
                size = 3, 
                shape = 20) +
    theme_linedraw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black", 
                                    face = "bold",
                                    size = strip.size,
                                    hjust = 0),
          axis.text.x = element_text(size = 7),
          panel.grid.major = element_blank(),
          panel.border = element_rect(size = 1.5, 
                                      colour = "black"),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_y_continuous(breaks = seq(-4,10,by=y.axis.breaks)) + 
    ylab("log2(Intensity)") +
    xlab("") + 
    scale_colour_manual(values = met.brewer("Egypt",4)) + 
    scale_fill_manual(values = met.brewer("Egypt",4)) 
    #scale_fill_viridis_d() +
    #scale_colour_viridis_d()
  
  # remove facet if bottom.plot == T
  if(bottom.plot == F){
    
    temp.prot.plot <- temp.prot.plot + 
      facet_wrap(~paste(#gene.name, "\n", 
      "Protein Abundance", sep = ""), scales = "free_y")
    
  }
  
  return(temp.prot.plot)
  
}


#### Plot single glycopeptide feature

gp.boxplot <- function(peak.number,
                       end.plot = F,
                       end.facet.name = NA,
                       strip.size = 10.5,
                       y.axis.breaks = 0.5,
                       input_colour){
  
  temp1 <- filter(oxoscan.quant, peak_num == peak.number) %>%
    distinct() %>%
    left_join(validated.hits, by = "peak_num")
  
  
  temp.plot <- temp1 %>%
    ggplot(aes(x = factor(F.Severity, levels=c("Healthy","Mild","Moderate", "Severe")),
               y = log2(Intensity), ### Glycopeptide feature intensity
               fill = input_colour)) + 
    geom_boxplot(inherit.aes = T, 
                 fill = NA, 
                 alpha=0.75, 
                 outlier.shape = NA) + 
    geom_jitter(aes(colour = F.Severity), 
                position = position_jitter(width = 0.1, height = 0.1),
                alpha = 0.75, 
                size = 3, 
                shape = 20) +
    theme_linedraw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(colour = "black", 
                                    face = "bold",
                                    size = strip.size,
                                    hjust = 0),
          strip.text.y = element_text(colour = "black", 
                                      face = "bold",
                                      size = 12),
          axis.text.x = element_text(size = 7),
          panel.grid.major = element_blank(),    
          panel.border = element_rect(size = 1.5, 
                                      colour = "light grey"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
          legend.position = "none") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) + 
    ylab("") +
    xlab("") + 
    scale_colour_manual(values = met.brewer("Egypt",4)) + 
    scale_fill_manual(values = met.brewer("Egypt",4)) 
 

  
    if(end.plot == T){ ### facet grid including gene name on right hand side
      
      temp.plot <- temp.plot +
        facet_grid(~paste(end.facet.name,sep="")~paste0(Genes, 
                                                        ", ",
                                                        AA, Position,
                                                        ", ",
                                                        Glycan.Short),
                   scales = "free_y") 
      
    }else{ ### just top labelling facet wrap
      
      temp.plot <- temp.plot +
        facet_wrap(~paste0(Genes, 
                           ", ",
                           AA, Position,
                           ", ",
                           Glycan.Short),
                   scales = "free_y", 
                   nrow = 1)
    }
  
  return(temp.plot)
  
}


### Protein-normalised glycopeptide abundance

gp.norm.boxplot <- function(peak.number,
                       end.plot = F,
                       end.facet.name = NA,
                       strip.size = 10.5,
                       y.axis.breaks = 0.5){
  
  temp1 <- filter(protein.normalised, peak_num == peak.number) %>%
    distinct() %>%
    mutate(GP.Normalised = Glycopeptide.Intensity / Protein.Intensity)
  
  
  temp.plot <- temp1 %>%
    ggplot(aes(x = factor(F.Severity, levels=c("Healthy","Mild","Moderate", "Severe")),
               y = log2(GP.Normalised), ### Normalised GP feature intensity
               fill = F.Severity)) + 
    geom_boxplot(inherit.aes = T, 
                 fill = NA, 
                 alpha=0.75, 
                 outlier.shape = NA) + 
    geom_jitter(aes(colour = F.Severity), 
                position = position_jitter(width = 0.1, height = 0.1),
                alpha = 0.75, 
                size = 3, 
                shape = 20) +
    theme_linedraw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(colour = "black", 
                                      face = "bold",
                                      size = strip.size,
                                      hjust = 0),
          strip.text.y = element_text(colour = "black", 
                                      face = "bold",
                                      size = 12,
                                      #hjust = 0.0
          ),
          axis.text.x = element_text(size = 7),
          panel.grid.major = element_blank(),    
          panel.border = element_rect(size = 1.5, 
                                      colour = "light grey"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
          legend.position = "none") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) + 
    ylab("") +
    xlab("") + 
    scale_colour_manual(values = met.brewer("Egypt",4)) + 
    scale_fill_manual(values = met.brewer("Egypt",4)) 
  
  if(end.plot == T){ ### facet grid including gene name on right hand side
    
    temp.plot <- temp.plot +
      facet_grid(paste(end.facet.name,sep="")~Glycan.Short,
                 scales = "free_y") 
    
  }else{ ### just top labelling facet wrap
    
    temp.plot <- temp.plot +
      facet_wrap(~paste0(Genes, 
                        ", ",
                        AA, Position,
                        ", ",
                        Glycan.Short),
                 scales = "free_y", 
                 nrow = 1)
  }
  
  return(temp.plot)
  
}



### Plot all glycoforms of protein (normalised) - faceted

gp.norm.boxplot.facet <- function(Protein,
                                  strip.size = 10.5,
                                  y.axis.breaks = 0.5,
                                  N.Row = 1){
  
  temp1 <- filter(protein.normalised, Genes == Protein) %>%
    distinct() %>%
    mutate(GP.Normalised = Glycopeptide.Intensity / Protein.Intensity)
  
  
  temp.plot <- temp1 %>%
    ggplot(aes(x = factor(F.Severity, levels=c("Healthy","Mild", "Moderate","Severe")),
               y = log2(GP.Normalised), ### Normalised GP feature intensity
               fill = F.Severity)) + 
    geom_boxplot(inherit.aes = T, 
                 fill = NA, 
                 alpha=0.75, 
                 outlier.shape = NA) + 
    geom_jitter(aes(colour = F.Severity), 
                position = position_jitter(width = 0.1, height = 0.1),
                alpha = 0.75, 
                size = 3, 
                shape = 20) +
    theme_linedraw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(colour = "black", 
                                      face = "bold",
                                      size = strip.size,
                                      hjust = 0),
          strip.text.y = element_text(colour = "black", 
                                      face = "bold",
                                      size = 12),
          axis.text.x = element_text(size = 7),
          panel.grid.major = element_blank(),    
          panel.border = element_rect(size = 1.5, 
                                      colour = "light grey"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
          legend.position = "none") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) + 
    ylab("Glycopeptide abundance (protein-normalised)") +
    xlab("Severity") + 
    scale_colour_manual(values = met.brewer("Egypt",4)) + 
    scale_fill_manual(values = met.brewer("Egypt",4)) +
    facet_wrap(~paste0(Genes, 
                       ": ",
                       AA, Position,
                       ", ",
                       Glycan.Short,
                       ", ",
                       peak_num),
                 scales = "free_y", 
                 nrow = N.Row)
  
  return(temp.plot)
  
}



############################################################################################################
#### B2B plotting - reading in files and reformatting tables as required for plotting functions below ######
############################################################################################################


plot.b2b.fromfolder <- function(peak,
                                base.folder.dir = paste0(raw.data.dir, "Validation/"), ## Folder with each peak folder contained
                                Max.PPM.Error = 20,
                                CID.HCD.Tolerance = 0.1, ## Tolerance for CID/HCD fragments (before later filtering by PPM Error)
                                remove_oxonium = F,
                                remove_immonium =  T,
                                label_text_size = 3,
                                facet_text_size = 17,
                                intensity_threshold_for_labelling = 0){  
  
  ## Get folder name
  single.folder.name <- to.plot %>%
    filter(Peak == peak) %>%
    pull(Folder)
  
  ### Glycan info
  temp.glycan.info <- true.pos %>%
    filter(peak_num == peak)
  
  ### Import annotation
  annot.raw <- dir_ls(paste0(base.folder.dir, single.folder.name),
                      regexp = "DirectExport") %>%
    fread()
  
  ### Intact peptide + oxonium - table with all theoretical from Byonic
  intact.fragments <- annot.raw %>%
    filter(str_detect(`1`, "mz")) %>% ## Intact peptides with regexp
    select(name, `1`) %>%
    mutate(`1` = as.numeric(sub(`1`, 
                                pattern = "mz=",
                                replacement = ""))) %>%
    rename(X.Peak = "1") ## Rename m/z as X.Peak
  
  # b/y ion m/z values - table with all theoretical from Byonic
  by.fragments <- annot.raw %>%
    filter(str_detect(name, "Theo")) %>% ## b/y ions with regexp
    pivot_longer(cols = 2:ncol(annot.raw),
                 values_to = "mz",
                 names_to = "position") %>%
    mutate(`name` = sub(`name`, pattern = ".*\\[(.*)\\].*", replacement = "\\1"), ## get b/y ion name
           `name` = gsub(`name`, pattern = "1d|iso0\\*|-$|\\%", replacement = ""), ## remove '1d', y1 iso, % and `-`
           mz = sub(mz, pattern = "^(\\d*\\.\\d\\d\\d\\d).*", replacement = "\\1")) %>% ## take iso0 only
    filter(str_detect(mz, "\\d")) ## remove empty rows
  
  
  # b/y ions reformat and combine with intact/oxonium to full list
  theoretical.fragments <- annot.raw %>%
    filter(str_detect(name, "Observed")) %>%
    mutate(`name` = gsub(`1`, pattern = "\\%", replacement = ""),
           `name` = gsub(`name`, pattern = "1d", replacement = ""), ## remove '1d' extra
           `1` = gsub(`1`, pattern = ".*(iso0\\*)", replacement = "\\1"), ## take y1 over to 1 column
           `name` = gsub(`name`, pattern = "iso0\\*", replacement = ""), ## remove y1 iso from name
           `name` = gsub(`name`, pattern = "-$", replacement = "")) %>% ## remove `-` at end
    pivot_longer(cols = 2:ncol(annot.raw),
                 values_to = "observed",
                 names_to = "position") %>%
    filter(is.na(observed) == F, ## remove empty rows
           str_detect(observed, "iso")) %>% ## 
    left_join(by.fragments, 
              by = c("name", "position")) %>%
    mutate(charge = ifelse(str_detect(name, "\\+\\+"), ## reformat to readable
                           yes = "++",
                           no = "+"),
           neutral.loss = sub(name, pattern = ".*([+-]\\d+)", replacement = "\\1"), ## reform netrual loss text
           neutral.loss = ifelse(str_detect(neutral.loss,"\\d"),
                                 yes=neutral.loss,
                                 no = ""),
           name = sub(name, pattern = "[+-]\\d*", replacement = ""),
           ion.label = paste(name,position,neutral.loss,charge,sep="")) %>% ## readable ion descriptions
    select(ion.label, mz) %>%
    rename(name = "ion.label",
           X.Peak = "mz") %>%
    rbind(intact.fragments) %>% ## join b/y fragments with intact peptide / peptide + glycan fragments
    mutate(X.Peak = as.numeric(X.Peak),
           iontype = ifelse(str_detect(name, "Pep|M"), ## use ifelse statements to categorise each type of fragment
                            yes = "Y-type",
                            no = ifelse(str_detect(name, "b|y|a\\d*"),
                                        yes = "b/y",
                                        no = ifelse(str_detect(name, "Hex|Neu"),
                                                    yes = "oxonium",
                                                    no = ifelse(str_detect(name, "C.*N."),
                                                                yes = "oxonium",
                                                                no = "immonium")))),
           name = gsub(name,
                       pattern = "HexNAc",
                       replacement = "N"), ## Shorten glycan names to make plotting/annotation easier
           name = gsub(name,
                       pattern = "Hex",
                       replacement = "H"),
           name = gsub(name,
                       pattern = "NeuAc",
                       replacement = "S"),
           name = gsub(name,
                       pattern = "M(_\\d\\+)(.*)",
                       replacement = "M\\2\\1"),
           name = gsub(name,
                       pattern = " ",
                       replacement = "")) %>%
    group_by(X.Peak) %>%
    arrange(X.Peak) %>% ## arrange by mz, keep only the first
    filter(row_number()==1) %>% ## remove duplicated labels / ions
    ungroup() 
  
  ## import HCD MS/MS spectrum txt files
  hcd.txt <- paste0(base.folder.dir, single.folder.name) %>%
    dir_ls(regexp = "HCD_peaks_1.csv") %>% ## Byonic outputs multiple files, only first contains peak info
    fread()
  
  
  ## import CID MS/MS spectrum txt files
  cid.txt <- paste0(base.folder.dir, single.folder.name) %>%
    dir_ls(regexp = "Q1.*txt") %>% ## Directly exported from PeakView (Sciex)
    fread() 
  
  ### get search rows only
  cid.hcd <- data.frame()
  
  ### loop throug each fragment in HCD spectrum txt file, compare against CID fragments with given Da tolerance
  for(j in seq_along(hcd.txt$V1)){
    
    temp.search.mz <- hcd.txt$V1[j]
    
    cid.hcd <- suppressMessages(rbind(cid.hcd,
                                      cid.txt[which.min(abs(cid.txt$`Mass/Charge`- temp.search.mz)),1:2] %>%
                                        subset(select=which(!duplicated(names(.)))) %>% # remove duplicate columns
                                        rename(CID.Peak = "Mass/Charge",
                                               CID.Area = "Area") %>%
                                        mutate(V1 = temp.search.mz) %>% ### matched m/z added here in loop
                                        left_join(hcd.txt) %>%
                                        rename(HCD.Peak = "V1",
                                               HCD.Area = "V2") %>%
                                        mutate(PPM.Error = 10^6 *abs(CID.Peak - HCD.Peak) / HCD.Peak,
                                               Abs.Error = abs(CID.Peak - HCD.Peak)) %>%
                                        filter(Abs.Error < CID.HCD.Tolerance)))
    
  }
  
  ## Fragments which match between HCD and CID spectra only - to be left joined
  cid.hcd.match <- cid.hcd %>%
    group_by(CID.Peak) %>%
    filter(Abs.Error == min(Abs.Error)) %>% ### keep only closest match
    ungroup() %>%
    filter(PPM.Error < Max.PPM.Error) %>% ### Tolerance for matching
    rename(CID.Abs.Error = "Abs.Error")
  
  
  ### Fragments matched between HCD spectrum and theoretical fragments
  hcd.frag.matched <- hcd.txt %>%
    rename(HCD.Peak = "V1",
           HCD.Area = "V2") %>%
    difference_left_join(theoretical.fragments %>%
                           rename(HCD.Peak = "X.Peak"),
                         by = "HCD.Peak",
                         max_dist = 0.1) %>% ## join with maximum tolerance acceptable (wide)
    as.data.frame() %>%
    rename(Annot.mz = "HCD.Peak.y",
           HCD.Peak = "HCD.Peak.x") %>%
    mutate(Annot.Abs.Error = abs(HCD.Peak - Annot.mz),
           Annot.PPM.Error = 10^6 * abs(HCD.Peak - Annot.mz) / HCD.Peak) %>% ## mass error for peak annotation
    group_by(Annot.mz) %>%
    filter(Annot.PPM.Error == min(Annot.PPM.Error)) %>% ### Keep only best match per annotation, remove others
    ungroup() %>%
    filter(Annot.PPM.Error < Max.PPM.Error)
  
  
  ## join matched cid peaks and matched theoretical fragment tables for plotting
  hcd.all.matches <- hcd.txt %>%
    rename(HCD.Peak = "V1",
           HCD.Area = "V2") %>%
    left_join(hcd.frag.matched) %>%
    left_join(cid.hcd.match) %>%
    distinct() %>%
    mutate(cut = cut(HCD.Peak, seq(100,2000, by = 40))) %>% ## take local maximum for annotation y value
    group_by(cut) %>%
    mutate(Annot.Y = max(HCD.Area)) %>%
    ungroup() %>%
    select(-cut) %>%
    mutate(CID.Peak = as.numeric(CID.Peak),
           Annot.Alpha = 1,
           Colour = "grey",
           Colour = ifelse(CID.Peak > 0 & !is.na(CID.Peak),
                           yes = "dark blue",
                           no = Colour),
           Colour = ifelse(str_detect(iontype, "oxonium") & !is.na(iontype),
                           yes = "dark red",
                           no = Colour)) %>% ## colour code by ion type
    mutate(name = ifelse(HCD.Area > intensity_threshold_for_labelling * max(HCD.Area), ## intensity threshold for annotation 
                         yes=name,
                         no=NA))
  
  ### Fragments matched between CID spectrum and theoretical fragments
  cid.frag.matched <- cid.txt %>%
    subset(select=which(!duplicated(names(.)))) %>% # remove duplicate columns
    rename(CID.Peak = "Mass/Charge",
           CID.Area = "Area") %>%
    select(CID.Peak, CID.Area) %>%
    difference_left_join(theoretical.fragments %>%
                           rename(CID.Peak = "X.Peak"),
                         by = "CID.Peak",
                         max_dist = 0.1) %>% ## join with maximum tolerance acceptable
    rename(Annot.mz = "CID.Peak.y",
           CID.Peak = "CID.Peak.x") %>%
    mutate(Annot.Abs.Error = abs(CID.Peak - Annot.mz), ## mass error
           Annot.PPM.Error = 10^6 * abs(CID.Peak - Annot.mz) / CID.Peak) %>% ## mass error for peak annotation
    group_by(Annot.mz) %>%
    filter(Annot.PPM.Error == min(Annot.PPM.Error)) %>% ### Keep only best match per annotation, remove others
    ungroup() %>%
    filter(Annot.PPM.Error < Max.PPM.Error)    
  
  
  ### Join matched hcd peaks and matched theoretical fragment tables for plotting
  cid.all.matches <- cid.txt %>%
    subset(select=which(!duplicated(names(.)))) %>% # remove duplicate columns
    rename(CID.Peak = "Mass/Charge",
           CID.Area = "Area") %>%
    select(CID.Peak, CID.Area) %>%
    left_join(cid.frag.matched) %>%
    left_join(cid.hcd.match) %>%
    distinct() %>%
    mutate(cut = cut(CID.Peak, seq(100,2000, by = 40))) %>% ## take local maximum for annotation y value
    group_by(cut) %>%
    mutate(Annot.Y = max(CID.Area)) %>%
    ungroup() %>%
    select(-cut) %>%
    mutate(Annot.Alpha = 1,
           Colour = "grey",
           Colour = ifelse(HCD.Peak > 0 & !is.na(HCD.Peak),
                           yes = "dark blue",
                           no = Colour),
           Colour = ifelse(str_detect(iontype, "oxonium") & !is.na(iontype),
                           yes = "dark red",
                           no = Colour)) %>% ## colour code by ion type
    mutate(name = ifelse(CID.Area > intensity_threshold_for_labelling * max(CID.Area), ## intensity threshold for annotation 
                         yes=name,
                         no=NA))
  
  
  
  ### Checks for PPM threshold of matches
  if(cid.all.matches$Annot.PPM.Error %>% max(na.rm = T) > Max.PPM.Error){
    return("Error with CID annotation mass accuracy thresholding - max(PPM.Error) too high")
  }
  
  if(cid.all.matches$PPM.Error %>% max(na.rm = T) > Max.PPM.Error){
    return("Error with CID (to HCD) mass accuracy thresholding - max(PPM.Error) too high")
  }
  
  if(hcd.all.matches$Annot.PPM.Error %>% max(na.rm = T) > Max.PPM.Error){
    return("Error with HCD annotation mass accuracy thresholding - max(PPM.Error) too high")
  }
  
  if(hcd.all.matches$PPM.Error %>% max(na.rm = T) > Max.PPM.Error){
    return("Error with HCD (to CID) mass accuracy thresholding - max(PPM.Error) too high")
  }
  
  
  if(remove_immonium == TRUE){
    
    hcd.all.matches <- hcd.all.matches %>%
      filter(is.na(iontype) | iontype != "immonium")
    
    cid.all.matches <- cid.all.matches %>%
      filter(is.na(iontype) | iontype != "immonium")
    
  }
  
  hcd.plot <- hcd.all.matches %>%
    ggplot(aes(x = HCD.Peak, 
               y = HCD.Area,
               colour = factor(Colour))) + 
    geom_col(size = 0.3) +
    geom_col(data = hcd.all.matches %>% filter(Colour == "dark blue"),
             colour = "dark blue",
             size = 0.3, 
             width = 0.3) + 
    theme_linedraw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),        
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = 0, size = facet_text_size),
          panel.grid = element_blank(),       
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(100,2000),
                       breaks = seq(0,2000,200),
                       expand = expansion(mult = c(0.01, 0.01)))  +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) + 
    scale_colour_identity() +
    facet_wrap(~paste0(#temp.glycan.info$Genes[1],
                       #", ",
                       temp.glycan.info$AA[1],
                       temp.glycan.info$Position[1], 
                       ", ",
                       temp.glycan.info$Glycan.Short[1],
                       "\n",
                       temp.glycan.info$Peptide[1])) +
    labs(x = "m/z",
         y = "HCD\nRelative Intensity") +
    geom_text_repel(aes(x = HCD.Peak,
                        y = HCD.Area,
                        label = name,
                        alpha = Annot.Alpha,
                        segment.alpha = Annot.Alpha),
                    size = label_text_size,
                    max.overlaps = 10,
                    segment.linetype = "dotted",
                    segment.size = 0.25,
                    force_pull = 0.1,
                    force = 4,
                    direction = "y") + 
    scale_alpha_identity()
  
  
  
  cid.plot <- cid.all.matches %>%
    ggplot(aes(x = CID.Peak, 
               y = CID.Area,
               colour = factor(Colour))) + 
    geom_col(size = 0.3) +
    geom_col(data = cid.all.matches %>% 
               filter(Colour == "dark blue"),
             inherit.aes = T,
             width = 0.3,
             size = 0.3
    ) + 
    theme_linedraw() + 
    theme(axis.line.x = element_line(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = 0, size = 14),
          strip.background = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(100,2000),
                       breaks = seq(0,2000,200),
                       position = "bottom",
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_reverse(expand = expansion(mult = c(0.05, 0.0))) + 
    scale_colour_identity() +
    #facet_wrap(~"", switch = "x") + ## To equalise plot borders top and bottom
    labs(x = "m/z",
         y = "CID\nRelative Intensity") +
    geom_text_repel(data = cid.all.matches,
                    aes(x = CID.Peak,
                        y = CID.Area,
                        label = name,
                        alpha = Annot.Alpha,
                        segment.alpha = Annot.Alpha),
                    size = label_text_size,
                    segment.linetype = "dotted",
                    segment.size = 0.25,
                    max.overlaps = 10,
                    force_pull = 0.1,
                    force = 4,
                    direction = "y") + 
    scale_alpha_identity()
  
  plot_grid(hcd.plot, 
            cid.plot,
            ncol = 1,
            label_size = 10,
            label_x = 0.45,
            label_y = 0.85)
  
  # Write out table for Source Data
  write_csv(cid.all.matches, file=paste0("/Users/whitem/Documents/GitHub/OxoScan-MS/paper_analysis/Figure3-5/SourceData/FigureS4/peak_num-",
                                         single.folder.name, 
                                         "_CID.csv"))
  
  write_csv(hcd.all.matches, file=paste0("/Users/whitem/Documents/GitHub/OxoScan-MS/paper_analysis/Figure3-5/SourceData/FigureS4/peak_num-",
                                         single.folder.name, 
                                         "_HCD.csv"))
  
}




#### B2B plotting matched only ######

plot.b2b.matched.only <- function(){
  
  plot <- plot_grid(
    
    cid.hcd.match %>%
      ggplot(aes(x = HCD.Peak, y = HCD.Area)) + 
      geom_col(width = 0.8) +
      theme_linedraw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),        
            panel.grid = element_blank(),   
            strip.background = element_blank(),
            strip.text = element_text(colour = "black", face = "bold", hjust = 0),            #panel.border = element_blank(),
            #panel.grid.major.x = element_line(),
            axis.text.y = element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.length = unit(0, "pt"),
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_x_continuous(limits = c(100,2000),
                         breaks = seq(0,2000,200),
                         expand = expansion(mult = c(0.01, 0.01)))  +
      scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) + 
      ylab("HCD\nRelative Intensity") +
      facet_wrap(~paste("Gene: ", 
                        pull(temp.glycan.info[1,"Genes"]), ## Gene
                        ", ",
                        pull(temp.glycan.info[1,"AA"]), ## Asn or Ser
                        pull(temp.glycan.info[1,"Position"]), ## Position
                        ", ", 
                        pull(temp.glycan.info[1,"Glycan.Short"]),
                        sep = ""))
    
    ,
    
    cid.hcd.match %>%
      ggplot(aes(x = CID.Peak, y = CID.Area))  + 
      geom_col(width = 0.8) +
      theme_linedraw() + 
      theme(#axis.title.x=element_blank(),
        #axis.text.x=element_blank(),   
        axis.line.x = element_line(),
        panel.border = element_rect(size=0.5),
        panel.grid = element_blank(),
        #panel.grid.major.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_x_continuous(limits = c(100,2000),
                         breaks = seq(0,2000,200),
                         position = "bottom",
                         expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_reverse(expand = expansion(mult = c(0.05, 0.0))) + 
      ylab("CID\nRelative Intensity")
    
    , 
    
    ncol = 1
    
    ,
    
    labels = 
      ,
    
    label_size = 10,
    label_x = 0.45,
    label_y = 0.85
    
  )
  
  return(plot)
  
}


#### Plot all fragments, coloured by matched / oxonium to display quality of match

plot.b2b.all <- function(remove_oxonium,
                         hcd_annot = T,
                         cid_annot = T,
                         remove_immonium =  T,
                         plot_matched_only = T) {
  
  if(plot_matched_only == T){
    
    hcd.matched <- hcd.matched %>%
      mutate(Annot.Alpha = ifelse(is.na(CID.Peak) == F & is.na(name) == F,
                                  yes = 1,
                                  no = 0))
    
    cid.matched <- cid.matched %>%
      mutate(Annot.Alpha = ifelse(is.na(HCD.Peak) == F & is.na(name) == F,
                                  yes = 1,
                                  no = 0))
    
    
  }
  
  
  if(remove_oxonium == TRUE){
    
    hcd.matched <- hcd.matched %>%
      filter(is.na(iontype) | iontype != "oxonium")
    
    cid.matched <- cid.matched %>%
      filter(is.na(iontype) | iontype != "oxonium")
    
  }
  
  if(remove_immonium == TRUE){
    
    hcd.matched <- hcd.matched %>%
      filter(is.na(iontype) | iontype != "immonium")
    
    cid.matched <- cid.matched %>%
      filter(is.na(iontype) | iontype != "immonium")
    
  }
  
  
  hcd.plot <- hcd.matched %>%
    ggplot(aes(x = HCD.Peak, 
               y = HCD.Area,
               colour = factor(Colour))) + 
    geom_col(size = 0.3) +
    geom_col(data = hcd.matched %>% filter(Colour == "dark green"),
             colour = "dark green",
             size = 0.3, 
             width = 0.3) + 
    theme_linedraw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),        
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = 0, size = 14),
          panel.grid = element_blank(),       
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(100,2000),
                       breaks = seq(0,2000,200),
                       expand = expansion(mult = c(0.01, 0.01)))  +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) + 
    scale_colour_identity() +
    facet_wrap(~paste(pull(temp.glycan.info[1,"Genes"]), ## Gene
                      ", ",
                      pull(temp.glycan.info[1,"AA"]), ## Asn or Ser
                      pull(temp.glycan.info[1,"Position"]), ## Position
                      ", ", 
                      pull(temp.glycan.info[1,"Glycan.Short"]),
                      sep = "")) +
    labs(x = "m/z",
         y = "HCD\nRelative Intensity")
  
  cid.plot <- cid.matched %>%
    ggplot(aes(x = CID.Peak, 
               y = CID.Area,
               colour = factor(Colour))) + 
    geom_col(size = 0.3) +
    geom_col(data = cid.matched %>% 
               filter(Colour == "dark green"),
             inherit.aes = T,
             width = 0.3,
             size = 0.3
    ) + 
    theme_linedraw() + 
    theme(axis.line.x = element_line(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(100,2000),
                       breaks = seq(0,2000,200),
                       position = "bottom",
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_reverse(expand = expansion(mult = c(0.05, 0.0))) + 
    scale_colour_identity() +
    labs(x = "m/z",
         y = "CID\nRelative Intensity")
  
  ### if phrases to add annotations
  if(hcd_annot == T){
    
    hcd.plot <- hcd.plot +
      geom_text_repel(aes(x = HCD.Peak,
                          y = HCD.Area,
                          label = name,
                          alpha = Annot.Alpha,
                          segment.alpha = Annot.Alpha),
                      size = 2,
                      max.overlaps = 25,
                      segment.linetype = "dotted",
                      segment.size = 0.25,
                      #segment.alpha = ,
                      #min.segment.length = 100,
                      force_pull = 5,
                      force = 1,
                      direction = "y") + 
      scale_alpha_identity()
    
    
  }
  
  if(cid_annot == T){
    
    cid.plot <- cid.plot +
      geom_text_repel(aes(x = CID.Peak,
                          y = CID.Area,
                          label = name,
                          alpha = Annot.Alpha,
                          segment.alpha = Annot.Alpha),
                      size = 2,
                      segment.linetype = "dotted",
                      segment.size = 0.25,
                      #min.segment.length = 100,
                      max.overlaps = 25,
                      force_pull = 5,
                      force = 1,
                      direction = "y") + 
      scale_alpha_identity()
  }
  
  return(plot_grid(hcd.plot, 
                   cid.plot,
                   ncol = 1,
                   label_size = 10,
                   label_x = 0.45,
                   label_y = 0.85))
  
  
  
  
}





#################################################################
#################### XIC plotting functions ######################
#################################################################


########## 2D XIC plotting function
plot.xic.2d <- function(peak,
                        remove_ions = c(),
                        plot_exclusion = T,
                        PlasmaPooledRep = 2){

  regex.filter <- "PlasmaPooled"
  
  regex.sub <- ".*PlasmaPooled_(\\d*).*"
  
  file.name <- true.pos %>%
    filter(peak_num == peak) %>%
    pull(Genes)
  
  ## RT for that peak in OxoScan output
  # merged_rt = reference RT
  ref.rt <- oxoscan.raw %>%
    filter(peak_num == peak,
           str_detect(spectrum_name, "PlasmaPooled_5")) %>%
    pull(merged_rt) %>%
    unique() 
  
  ## m/z for peak in OxoScan output
  ref.mz <- oxoscan.raw %>%
    filter(peak_num == peak,
           str_detect(spectrum_name, "PlasmaPooled_5")) %>%
    pull(merged_mz) %>%
    unique() 
  
  ## get accurate precursor m/z for horizontal line from dia.precursors
  acc.mz <- dia.precursors %>%
    filter(peak_num == peak) %>%
    pull(`m/z`) %>%
    round(3)
  
  
  ### RT alignment for ellipse
  align <- oxoscan.raw %>%
    mutate(Sample = sub(spectrum_name,
                        pattern = regex.sub,
                        replacement = "\\1")) %>%
    filter(Sample %in% PlasmaPooledRep, ## Take only single QC sample
           peak_num == peak) %>% ## filter for peak of choice
    select(Sample, merged_rt, aligned_rt, merged_mz) %>%
    distinct() %>%
    mutate(Sample = as.character(Sample))
  
  ### make reference for alpha to be +/- 0.1min around peak centre 
  ### (prevent nearby peaks drowning out signal)
  xic.alpha.ref <- xic.raw %>%
    filter(mz > ref.mz - 10, ## take only local area around peak of interest
           mz < ref.mz + 10,
           RT > ref.rt - 0.1,  
           RT < ref.rt + 0.1) %>% ## scale intensity to within plotting area
    pull(Intensity) %>%
    max()

  
  xic <- xic.raw %>%
    filter(mz > ref.mz - 60, ## take only local area around peak of interest
           mz < ref.mz + 60,
           RT > ref.rt - 0.5,  
           RT < ref.rt + 0.5) %>%
    mutate(Sample = as.character(PlasmaPooledRep)) %>%
    left_join(align, by = "Sample") %>%
    mutate(RT.Normalised = RT + (merged_rt - aligned_rt), ## sample alignment
           Scaled = Intensity / xic.alpha.ref, ## scale intensity to within quant area
           Scaled = ifelse(Scaled > 1, ## optional: improves visualisation of low abundant
                           1,  ## peaks with nearby high abundant peaks
                           Scaled)) 
  
  max.rt <- xic %>%
    pull(merged_rt) %>%
    unique() 

  
  rt.adjustment <- xic$merged_rt[1] - xic$aligned_rt[1]
  
  
  
  #### Plot
  
  xic.plot <- xic %>%
    ggplot(aes(x = RT.Normalised,
               y = mz,
               fill = Scaled,
               alpha = Scaled)) + 
    geom_tile(height = 2,
              width = 0.025) +
    geom_hline(yintercept = acc.mz, 
               linetype = "dashed",
               colour = "#859b6c",
               size = 1.5,
               alpha = 0.5) + 
    theme_linedraw() + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = -0.01, size = 14),
          plot.margin = unit(c(0, 0.2, 0, 0.2), "cm")) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.0)),
                       limits = c(ref.mz-25, ref.mz+25)) +
    facet_wrap(~"") + 
    xlab("Retention Time") + 
    ylab("Precursor m/z") + 
    scale_alpha_identity() + 
    scale_x_continuous(expand = expansion(mult = c(0.0, 0.0)),
                       limits = c(ref.rt - 0.5, ref.rt + 0.5)) +
    scale_fill_gradientn(values = c(0,0.5,1), 
                         colours = c("white", "grey", "black")) +
    geom_ellipse(data = data.frame(),
                 inherit.aes = F,
                 aes(x0 = max.rt,
                     y0 = ref.mz,
                     a = 3.5*0.025, b = 11,
                     angle = 0),
                 colour = "#61599d",
                 fill = NA,
                 alpha = 0.4, 
                 size = 1,
                 linetype = "dashed") + 
    geom_ellipse(data = data.frame(),
                 inherit.aes = F,
                 aes(x0 = max.rt,
                     y0 = ref.mz,
                     a = 0.2625, 
                     b = 22,
                     angle = 0),
                 colour = "#c36377", ##
                 fill = NA,
                 alpha = 0.4, 
                 size = 1,
                 linetype = "dashed")
  
  return(xic.plot)
  
}


########################################################################
############ MS1 / Q1 overlay ##########################################
########################################################################



plot.q1.ms1 <- function(peak,
                        Inset.x = 0.45){
  
  ## Pull folder for dir
  folder <- to.plot %>%
    filter(Peak == peak) %>%
    pull(Folder)
  
  ## Get glycopeptide info for peak of interest
  glycan.info <- true.pos %>%
    filter(peak_num == peak)
  
  ## Set folder for peak of interest
  dir <- paste0(raw.data.dir, "Validation/", folder) ## replace with [1,i] in loop
  
  ## List the individual ion traces exported from PeakView
  fragments <- dir_ls(dir, regexp = "Fragment Ion")
  
  ## Number of fragments to define colour scale
  n.frag <- length(fragments)
  
  ## Read in Q1 traces from each fragment and combine into one dataframe
  for(i in seq_along(fragments)){
    
    if(i==1){ ## For first loop, make empty data.frame to add q1.temp to
      q1 <- data.frame()
    }
    
    ### Fragment name/mz
    frag <- fragments[i] %>%
      as.character() %>%
      sub(pattern = ".*\\((.*)\\).*",
          replacement = "\\1") %>%
      as.numeric() %>%
      round(2)
    
    ### Read in and add label
    q1.temp <- fragments[i] %>%
      fread() %>%
      mutate(frag = frag)
    
    
    q1 <- rbind(q1,
                q1.temp)
    
    rm(q1.temp) ## Leave no trace
    
    if(i == length(fragments)){ ## Remove temp data.frames
      
      q1.scale <- q1 %>%
        mutate(Scaled.Intensity = 100*intensity / max(intensity)) ### Intensity == %age intensity
      
      rm(q1) ## Remove previous data.frame
      
    }
    
  }
  
  
  ## Read in MS1 spectrum
  ms1.raw <- dir_ls(dir, regexp = "MS1_500msMethod") %>%
    fread() %>%
    subset(select=which(!duplicated(names(.)))) %>% ## remove duplicate columns
    mutate(Scaled.Height = 100*Height/max(Height)) ## scale intensities
  
  
  ## Get accurate mass for peak of interest
  acc.mz <- dia.precursors %>%
    filter(peak_num == peak) %>%
    pull(`m/z`) %>%
    round(2)
  
  
  ### For scaling MS1 and Q1 plots: 
  ## Take scaled MS1 height (as percentage) for peak of interest
  scaled.ms1.height <- ms1.raw %>%
    filter(`Mass/Charge` > acc.mz - 1, 
           `Mass/Charge` < acc.mz + 1) %>%
    pull(Scaled.Height) %>%
    max()
  
  
  
  
  ## Find maximum oxonium intensity around acc.mz
  max.ox.height <- q1.scale %>%
    filter(Precursor < acc.mz + 5, 
           Precursor > acc.mz - 5) %>%
    pull(Scaled.Intensity) %>%
    max()
  
  ## Use 100 / max.ox.height to account for other taller oxonium maxima in plotting
  
  
  ## Scale Q1 intensities to be either 20% (if scaled.ms1.height < 15) or 1.25 * scaled.ms1.height
  # Keeps plotting clear
  
  q1.scaling.factor <- ifelse(scaled.ms1.height > 15,
                              yes = scaled.ms1.height * 1.5,
                              no = 20) / 100
  
  q1.inset.factor <- 1.2 * scaled.ms1.height / 100
  
  box.height <- 1.35 * scaled.ms1.height / 100
  
  
  ### Plot whole spectrum plot
  main.plot <- q1.scale %>%
    ggplot(aes(x = Precursor,
               y = 1.25 * (100 / max.ox.height) * q1.scaling.factor * Scaled.Intensity,
               colour = factor(frag),
               fill = factor(frag))) + 
    geom_line(size = 0.25,
              alpha = 0.75) + 
    geom_col(data = ms1.raw,
             aes(x = `Mass/Charge`,
                 y = `Scaled.Height`),
             width = 0.005,
             colour = "black") +
    theme_linedraw() + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", size = 18),
          strip.text.x = element_text(hjust = -0.01)) + 
    scale_colour_manual(values = pnw_palette("Sunset2", n.frag)) +
    scale_fill_manual(values = pnw_palette("Sunset2", n.frag)) +
    scale_x_continuous(expand = expansion(mult = c(0.0, 0.0)),
                       #breaks = seq(0,70,10),
                       limits=c(800,1400)) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) + 
    xlab("m/z") + 
    ylab("Relative Intensity") +
    labs(colour = "Ion",
         fill = "Ion") + 
    facet_wrap(~paste0(glycan.info[1, "Genes"],
                       ", ",
                       glycan.info[1, "AA"],
                       glycan.info[1, "Position"],
                       ", Precursor m/z = ",
                       acc.mz)) + 
    geom_rect(data = data.frame(),
              inherit.aes = F,
              aes(xmin = acc.mz - 10,
                  xmax = acc.mz + 10,
                  ymin = 0, 
                  ymax = box.height),
              alpha = 0.25, 
              fill = "dark grey",
              colour = "black",
              size = 0.2)
  
  
  
  
  ### inset #1 
  inset <- q1.scale %>%
    ggplot(aes(x = Precursor,
               y = q1.inset.factor * Scaled.Intensity,
               colour = factor(frag),
               fill = factor(frag))) + 
    geom_line(size = 0.75,
              alpha = 0.75) + 
    geom_col(inherit.aes = F,
             data = ms1.raw,
             aes(x = `Mass/Charge`,
                 y = `Scaled.Height`),
             width = 0.005,
             colour = "black") +
    theme_linedraw() + 
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "#F4F4F4"),
          legend.position = "none",
          axis.title = element_blank()) + 
    scale_colour_manual(values = pnw_palette("Sunset2", n.frag)) +
    scale_fill_manual(values = pnw_palette("Sunset2", n.frag)) + 
    scale_x_continuous(expand = expansion(mult = c(0.0, 0.0)),
                       breaks = acc.mz,
                       limits=c(acc.mz - 8, acc.mz + 8)) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.05)),
                       breaks = c(0, round(scaled.ms1.height,0)),
                       limits = c(0, 1.25 * scaled.ms1.height)
    ) + 
    xlab("") + 
    ylab("") 
  
  
  ggdraw() + 
    draw_plot(main.plot) + 
    draw_plot(inset, x = Inset.x, y = 0.35, height = 0.5, width = 0.35)
  
}


########################################################################
################## Skyline/OxoScan MSMS B2B Plot #######################
########################################################################


prm.oxoscan.b2b.msms <- function(JoinedDFName,
                                 TopPlotLabel,
                                 SetWidth,
                                 LabelXPos,
                                 LabelYPos,
                                 LabelSize,
                                 PeakAnnotationSize = 4,
                                 MinLabelIntensity){
  
  ## Works with a pre-made table
  temp.prm.oxoscan.b2b.df <- get(JoinedDFName) %>%
    mutate(name = ifelse(Rel.Intensity > MinLabelIntensity,
                         yes=name,
                         no=NA))
  
  prm.plot <- temp.prm.oxoscan.b2b.df %>%
    filter(Method == "PRM") %>%
    ggplot(aes(x = `m/z`, 
               y = Rel.Intensity)) + 
    geom_col(width=SetWidth) +
    geom_text_repel(aes(x = `m/z`,
                        y = Rel.Intensity,
                        label = name),
                    colour = "grey",
                    size = PeakAnnotationSize,
                    #min.segment.length = 0.1,
                    #max.overlaps = 25,
                    #point.padding = 2,
                    nudge_y = 40,
                    segment.linetype = "dotted",
                    #segment.size = 0.5,
                    #force_pull = 1,
                    force = 3,
                    direction = "y") +
    theme_linedraw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),        
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = 0, size = LabelSize*3),
          panel.grid = element_blank(),       
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(0,2000),
                       breaks = seq(0,2000,200),
                       expand = expansion(mult = c(0.01, 0.01)))  +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.15))) + 
    geom_text(inherits.aes = F,
              aes(x = LabelXPos, 
                  y = LabelYPos,
                  label = "HR-MRM, ZenoTOF 7600"),
              colour = "black",
              size=LabelSize,
              hjust=1) + 
    facet_wrap(~paste(TopPlotLabel)) + 
    labs(x = "m/z",
         y = "HR-MRM\nIntensity / %") +
    scale_colour_identity()
  
  
  oxoscan.plot <- temp.prm.oxoscan.b2b.df %>%
    filter(Method == "OxoScan") %>%
    ggplot(aes(x = `m/z`, 
               y = Rel.Intensity)) + 
    geom_col(width=SetWidth) +
    theme_linedraw() +
    theme(axis.title.x=element_text(size=15),
          #axis.text.x=element_blank(),        
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", face = "bold", hjust = 0, size = LabelSize*3),
          panel.grid = element_blank(),       
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 15),
          axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = unit(c(0, 1, 0, 1), "cm")) +
    scale_x_continuous(limits = c(0,2000),
                       breaks = seq(0,2000,250),
                       expand = expansion(mult = c(0.01, 0.01)))  +
    scale_y_reverse(expand = expansion(mult = c(0.15, 0.0))) + 
    labs(x = "m/z",
         y = "OxoScan\nIntensity / %") + 
    geom_text(aes(x = LabelXPos, 
                  y = LabelYPos,
                  label = "OxoScan-MS, TripleTOF 6600"),
              size = LabelSize,
              hjust=1) + 
    scale_colour_identity()
  
  
  plot_grid(prm.plot, 
            oxoscan.plot,
            ncol = 1,
            rel_heights = c(1,0.95))
  
  

}
