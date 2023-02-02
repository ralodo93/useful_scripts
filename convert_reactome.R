library(ReactomePA)
library(reactome.db)

get_Reactome_Env <- function() {
  if (!exists(".ReactomePA_Env", envir = .GlobalEnv)) {
    assign(".ReactomePA_Env", new.env(), .GlobalEnv)
  }
  get(".ReactomePA_Env", envir= .GlobalEnv)
}

get_Reactome_DATA <- function(organism = "human", list_ann_id) {
  
  ALLEG <- getALLEG(organism)
  
  EXTID2PATHID <- as.list(reactomeEXTID2PATHID)
  EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% ALLEG]
  
  PATHID2EXTID <- as.list(reactomePATHID2EXTID) ## also contains reactions
  PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% list_ann_id]
  
  PATHID2NAME <- as.list(reactomePATHID2NAME)
  PI <- names(PATHID2NAME)
  ## > PATHID2NAME[['68877']]
  ## [1] "Homo sapiens: Mitotic Prometaphase" "Homo sapiens: Mitotic Prometaphase"
  PATHID2NAME <- lapply(PATHID2NAME, function(x) x[1])
  names(PATHID2NAME) <- PI
  
  PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME)]
  PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))]
  PATHID2EXTID <- lapply(PATHID2EXTID, function(x) intersect(x, ALLEG))
  
  PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]
  PATHID2NAME <- lapply(PATHID2NAME, function(x){gsub("^\\w+\\s\\w+:\\s+", "", x)})
  export_results <- data.frame(Term = names(PATHID2NAME), Description = unlist(PATHID2NAME))
}

getALLEG <- function(organism) {
  annoDb <- getDb(organism)
  require(annoDb, character.only = TRUE)
  annoDb <- eval(parse(text=annoDb))
  eg=keys(annoDb, keytype="ENTREZID")
  return(eg)
}

getDb <- function(organism) {
  if (organism == "worm") {
    organism = "celegans"
    warning("'worm' is deprecated, please use 'celegans' instead...")
  }
  
  annoDb <- switch(organism,
                   anopheles   = "org.Ag.eg.db",
                   arabidopsis = "org.At.tair.db",
                   bovine      = "org.Bt.eg.db",
                   canine      = "org.Cf.eg.db",
                   celegans    = "org.Ce.eg.db",
                   chicken     = "org.Gg.eg.db",
                   chimp       = "org.Pt.eg.db",
                   coelicolor  = "org.Sco.eg.db", 
                   ecolik12    = "org.EcK12.eg.db",
                   ecsakai     = "org.EcSakai.eg.db",
                   fly         = "org.Dm.eg.db",
                   gondii      = "org.Tgondii.eg.db",
                   human       = "org.Hs.eg.db",
                   malaria     = "org.Pf.plasmo.db",
                   mouse       = "org.Mm.eg.db",
                   pig         = "org.Ss.eg.db",
                   rat         = "org.Rn.eg.db",
                   rhesus      = "org.Mmu.eg.db",
                   xenopus     = "org.Xl.eg.db",
                   yeast       = "org.Sc.sgd.db",
                   zebrafish   = "org.Dr.eg.db",
  )
  return(annoDb)
}

Reactome_DATA <- get_Reactome_DATA(organism = "human", list_ann_id = c("R-HSA-1989781"))

                                   