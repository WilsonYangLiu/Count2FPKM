rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/2+.Count2FPKM/')

gencode <- read.csv(file = './GRCm38.p5/GCF_000001635.25_GRCm38.p5_genomic.tsv', sep = '\t', comment.char = '#')
exon <- gencode[gencode$feature == 'exon', ]
counts <- read.csv(file = './count.txt', sep = '\t', row.names = 1)

exon <- exon[exon$gene %in% rownames(counts), ]

# Find gene in gene's synonym name
not.in.gencode <- rownames(counts)[!(rownames(counts) %in% gencode$gene)]
#not.in.exon <- rownames(counts)[!(rownames(counts) %in% exon$gene)]
tmp <- sapply(not.in.gencode, function(x) {grep(x, gencode$gene_synonym)})
selt.gencode <- data.frame()
for (i in 1:length(tmp)) {
  if (length(tmp[[i]]) > 0) {
    tmp.frame <- gencode[tmp[[i]], ]
    tmp.frame$match <- names(tmp[i])
    selt.gencode <- rbind(selt.gencode, tmp.frame)
  }
}
flag <- sapply(1:dim(selt.gencode)[1], function(x) {selt.gencode$match[x] %in% stringr::str_split(selt.gencode$gene_synonym[x], ',')[[1]]})
selt.gencode <- selt.gencode[flag, ]
selt.gencode$feature <- paste('match', selt.gencode$gene, sep = '.')
selt.gencode$gene <- selt.gencode$match; selt.gencode$match <- NULL
exon <- rbind(exon, selt.gencode)
exon$lenght <- exon$end - exon$start + 1
exon <- exon[order(exon$seqname), ]

not.in.gencode <- not.in.gencode[!(not.in.gencode %in% unique(sort(selt.gencode$gene)))]
write.csv(not.in.gencode, file = 'count.not.in.gencode.csv')

# Calculate FPKM
require(DESeq2)
sequence.name <- read.csv(file = './GRCm38.p5/sequence-region', header = FALSE, sep = ' ', row.names = 1)
colnames(sequence.name) <- c('start', 'end')
sequence.name$lenght <- (sequence.name$end - sequence.name$start + 1)
seq.info <- Seqinfo(seqnames = rownames(sequence.name), seqlengths = sequence.name$lenght, genome = 'GRCm38.p5')

gr <- GRangesList()
gene.set <- as.character(unique(exon$gene))
for (gene in gene.set) {
  gene.info <- exon[exon$gene == gene, ]
  x <- as.data.frame(table(gene.info$seqname))
  gene.gr <- GRanges(Rle(values = x$Var1, lengths = x$Freq), 
                     IRanges(start = gene.info$start, end = gene.info$end, names = gene.info$gene),
                     seqinfo = seq.info)
  gr[[gene]] <- gene.gr
}
# x <- as.data.frame(table(exon$seqname))
# gr <- GRanges(Rle(values = x$Var1, lengths = x$Freq), 
#               IRanges(start = exon$start, end = exon$end, names = exon$gene),
#               seqinfo = seq.info)

design <- data.frame(condition = colnames(counts), row.names = colnames(counts))
design$condition <- sub('X(.*)W[CN].*', 'w\\1C', design$condition)
design$condition <- sub('X(.*)WT.*', 'w\\1t', design$condition)
design$condition <- as.factor(design$condition)
dds <- DESeqDataSetFromMatrix(countData = counts[gene.set, ], 
                              colData = design, 
                              design = ~condition)
rowRanges(dds) <- gr
save(dds, file = 'dds(count.in.gencode).RData')

fpkmData <- fpkm(dds)
write.csv(fpkmData, file = 'fpkm.csv')
detach('package:DESeq2')


