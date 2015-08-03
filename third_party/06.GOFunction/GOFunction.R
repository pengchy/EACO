#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

library(GOFunction)
library(org.Lmigratoria.eg.db)

interestGenes <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/G5-VS-S5.down",what="character",sep="\n")
refGenes <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/refids",what="character",sep="\n")

sigTerm <- GOFunction(interestGenes, refGenes, organism="org.Lmigratoria.eg.db",
		 ontology="BP", fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05,
		 poth=0.05, peth=0.05, bmpSize=2000, filename="sigTerm")

data(exampledata)
sigTerm <- GOFunction(interestGenes, refGenes, organism="org.Hs.eg.db",
		 ontology="BP", fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05,
		 poth=0.05, peth=0.05, bmpSize=2000, filename="sigTerm")



> GOFunction
function (interestGenes, refGenes, organism = "org.Hs.eg.db", 
		    ontology = "BP", fdrmethod = "BY", fdrth = 0.05, ppth = 0.05,
				pcth = 0.05, poth = 0.05, peth = 0.05, bmpSize = 2000, filename = "sigTerm") 
{
	require(organism, character.only = TRUE) || stop(paste("package",
				organism, "is required", sep = " "))
	.sql <- paste("select distinct t1.gene_id,t2.go_id", " from genes as t1 inner join",
			paste("go", tolower(ontology), "all", sep = "_"), " as t2 on t1._id=t2._id", seq = "")
	.sql <- "show tables "
	organism <- strsplit(organism, ".db")
	organism <- organism[[1]]
	conn <- get(paste(organism, "_dbconn", sep = ""))()
	generalAnn <- dbGetQuery(conn, .sql)
	annRef <- generalAnn[generalAnn[, 1] %in% refGenes, ]
	annInterest <- generalAnn[generalAnn[, 1] %in% interestGenes,]
	cat("Finding statistically significant terms...\n")
	termInfo <- enrichmentFunction(annRef, annInterest, fdrmethod,
			fdrth)
	sigTerm <- termInfo$sigTerm
	if (nrow(sigTerm) == 0) {
		warning("There is no significant term! \n")
			return(NULL)
	}

	allTerm <- termInfo$allTerm
	    require("GO.db") || stop("package GO.db is required")
	    conn <- get("GO_dbconn")()
	    .sql <- paste("select distinct go_id goid,term name from go_term where ontology='", 
	        toupper(ontology), "'", sep = "")
	    allTermName <- dbGetQuery(conn, .sql)
	    sigTermName <- allTermName[allTermName[, 1] %in% sigTerm[, 
	        1], ]
	    sigTermName <- sigTermName[order(sigTermName[, 1]), ]
	    sigTerm <- sigTerm[order(sigTerm[, 1]), ]
	    sigTerm$name <- sigTermName[, 2]
	    sigTerm <- sigTerm[, c(1, 6, 2, 3, 4, 5)]
	    .sql <- paste("select distinct t1.go_id parentid,t2.go_id childid from ", 
	        paste("go", tolower(ontology), "offspring", sep = "_"), 
	        " as t3 inner join  go_term as t1 on t1._id=t3._id inner join go_term as t2", 
	        " on t2._id=t3._offspring_id", sep = "")
	    allTermRelation <- dbGetQuery(conn, .sql)
	    sigTermRelation <- allTermRelation[(allTermRelation[, 1] %in% 
	        sigTerm[, 1]) & (allTermRelation[, 2] %in% sigTerm[, 
	        1]), ]
	    rm(allTermRelation)
	    cat("Treating for local redundant terms...\n")
	    sigTerm_LocalRedun <- localRedundancy(sigTerm, generalAnn, 
	        sigTermRelation, annRef, annInterest, ppth, pcth)
	    cat("Treating for global redundant terms...\n")
	    sigTerm_GlobalRedun <- globalRedundancy(generalAnn, sigTermRelation, 
	        annRef, annInterest, sigTerm_LocalRedun, poth, peth)
	    cat("Visualizing the GO DAG...\n")
	    require("graph") || stop("package graph is required")
	    sigDAG <- createGODAG(as.character(sigTerm[, 1]), ontology)
	    allDAGTerm <- allTerm[allTerm[, 1] %in% nodes(sigDAG), ]
	    dagTermName <- allTermName[allTermName[, 1] %in% allDAGTerm[, 
	        1], ]
	    dagTermName <- dagTermName[order(dagTermName[, 1]), ]
	    allDAGTerm <- allDAGTerm[order(allDAGTerm[, 1]), ]
	    allDAGTerm$name <- dagTermName[, 2]
	    allDAGTerm <- allDAGTerm[, c(1, 6, 2, 3, 4, 5)]
	    sigTermID <- as.character(sigTerm[, 1])
	    sigTerm_LocalRedunID <- as.character(sigTerm_LocalRedun[, 
	        1])
	    sigTerm_GlobalRedunID <- as.character(sigTerm_GlobalRedun[, 
	        1])
	    showSigNodes(sigDAG, sigTermID, sigTerm_LocalRedunID, sigTerm_GlobalRedunID, 
	        allDAGTerm, bmpSize, filename)
	    label <- array("", dim = c(nrow(sigTerm), 1))
	    rmsigLocalTerm <- setdiff(sigTermID, sigTerm_LocalRedunID)
	    label[sigTermID %in% rmsigLocalTerm, 1] <- "Local"
	    rmsigGlobalTerm <- setdiff(sigTerm_LocalRedunID, sigTerm_GlobalRedunID)
	    label[sigTermID %in% rmsigGlobalTerm, 1] <- "Global"
	    label[sigTermID %in% sigTerm_GlobalRedunID, 1] <- "Final"
	    sigTerm$FinalResult <- label
	    tablename <- paste(filename, ".csv", sep = "")
	    write.csv(sigTerm, tablename, row.names = F)
	    cat("\n********************  Results  ********************\n")
	    cat("The number of annotated interesting genes:", length(unique(annInterest[, 
	        1])), "\n")
	    cat("The number of annotated reference genes:", length(unique(annRef[, 
	        1])), "\n")
	    cat("The number of statistically significant terms:", nrow(sigTerm), 
	        "\n")
	    cat("The number of terms after treating local redundancy:", 
	        nrow(sigTerm_LocalRedun), "\n")
	    cat("The number of terms after treating global redundancy:", 
	        nrow(sigTerm_GlobalRedun), "\n")
	    note = paste("Please see details about the significant terms in the files ", 
	        filename, ".csv and ", filename, ".bmp!\n", sep = "")
	    cat(note)
	    return(sigTerm_GlobalRedun)
}
