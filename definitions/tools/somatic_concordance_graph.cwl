#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Generate somatic concordance report graph"
baseCommand: ["/usr/local/bin/Rscript", "somatic_concordance_graph.R"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/mgi-cle/somatic-concordance-report:v1"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'somatic_concordance_graph.R'
        entry: |
            args <- commandArgs(TRUE)
            outDir <- args[1]
            sompy <- args[2]
            vaf <- args[3]
            
            require(ggplot2, quietly = T)
            require(dplyr, quietly = T)
            require(reshape2, quietly = T)
            require(gridExtra, quietly = T)
            
            cmp <- read.table(sompy, header = T, stringsAsFactors = F)
            vafs <- read.table(vaf, header = T, stringsAsFactors = F)
            vafs$Concordance <- do.call(paste0, vafs[,c("Base", "Query")])
            vafs$type <- ifelse(nchar(vafs$Ref) == nchar(vafs$Alt), "snv", "indel")
            vafs <- vafs %>% mutate(Concordance = case_when(Concordance == "11" ~ "Concordant",
                                                            Concordance == "01" ~ "New",
                                                            Concordance == "10" ~ "Missed"))
            
            vafs$Concordance <- factor(vafs$Concordance, levels = c("Concordant", "Missed", "New"), ordered = T)
            
            p1 <- ggplot() + geom_point(data=vafs, aes(x=Base_Tumor_Vaf, y=Query_Tumor_Vaf, color=Concordance, shape=type)) +
              geom_hline(yintercept = 0.05, linetype = "dashed") +
              geom_vline(xintercept = 0.05, linetype = "dashed") +
              scale_color_manual(values = c("black", "blue", "red"))
            
            p2 <- ggplot() + geom_jitter(data=vafs, aes(x=Base_Normal_Vaf, y=Query_Normal_Vaf, color=Concordance, shape=type)) +
              geom_hline(yintercept = 0.01,linetype = "dashed") +
              geom_vline(xintercept = 0.01,linetype = "dashed") +
              scale_color_manual(values = c("black", "blue", "red"))
            
            dat <- melt(vafs[,grep("Vaf|Concordance", colnames(vafs))])
            p3 <- ggplot() + geom_boxplot(data=dat, aes(x=variable, y=value, fill=variable)) + facet_grid(. ~ Concordance)
            
            dat <- melt(vafs[,grep("Dp|Concordance",colnames(vafs))])
            p4 <- ggplot() + geom_boxplot(data=dat, aes(x=variable, y=value, fill=variable)) + facet_grid(. ~ Concordance)
           
            outpdf <- file.path(outDir, "somatic_concordance_graph.pdf") 
            pdf(file=outpdf, width=16, height=12)
            grid.arrange(p1,p2,p3,p4, layout_matrix=matrix(1:4, nrow = 2, byrow = T),widths = c(1,1), heights = c(.5,.7))
            
            dev.off()

arguments: [$(runtime.outdir)]

inputs:
    sompy_file:
        type: File
        inputBinding:
            position: 1
    vaf_file:
        type: File
        inputBinding:
            position: 2
outputs:
    out_pdf:
        type: File
        outputBinding:
            glob: "somatic_concordance_graph.pdf"

