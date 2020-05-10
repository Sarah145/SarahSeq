pal1 <- c('#DAF7A6', '#FF5733', '#900C3F', '#43077B')
pal2 <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
pal3 <- c("#0D0887FF", "#9C179EFF", "#ED7953FF", "#F0F921FF")
pal4 <- c('#FFFF00', '#FF0000', '#00FF00', '#0000FF', '#000000')
pals <- c(list(pal1), list(pal2), list(pal3), list(pal4))
hil_pal1 <- colorRampPalette(c("#0c0026", "#56106EFF", "#BB3754FF", "#BB3754FF", "#F98C0AFF", "#F98C0AFF", "#F98C0AFF", "#FCFFA4FF", "#FCFFA4FF", "#FCFFA4FF", "#FCFFA4FF"))
hil_pal2 <- colorRampPalette(c("#440154FF", "#3B528BFF",  "#21908CFF", "#21908CFF", "#5DC863FF", "#5DC863FF", "#5DC863FF", "#FDE725FF", "#FDE725FF", "#FDE725FF", "#FDE725FF", "#FDE725FF"))
hil_pals <- c(list(hil_pal1), list(hil_pal2))

goi_plot <- function(plot_df, gene_pos, goi){
  ggplotly(
  ggplot(plot_df, aes(x = position, y= num_alleles, 
                      text  = paste('Position:', position, 
                                    '\nVariant:', variant,
                                    '\nNumber of Alleles:', num_alleles,
                                    '\nrsID:', id))) + 
    geom_point(aes(col = factor(ifelse(str_detect(as.character(variant), '^[ACGT]$'), as.character(variant), 'Other'), levels = c('A', 'C', 'G', 'T', 'Other')))) +
    geom_segment(aes(x=position, xend=position, y=0, yend=num_alleles, text = NULL), size = 0.1, col = 'gray50', show.legend = F) +
    geom_hline(yintercept = 0, size = 10, col = 'gray5', lty = 1) +
    scale_x_continuous(limits = c(gene_pos$start_position, gene_pos$end_position)) +
    scale_y_continuous(limits = c(0,2)) +
    scale_color_manual(values = pal4) +
    labs(x= paste('Chromosome', gene_pos$chromosome_name), y = NULL, col = NULL, title = paste0(goi, '\n', nrow(plot_df), ' variants')) +
    theme_classic() +
    theme(legend.position = 'bottom',
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()), tooltip = 'text')}

genomes_1k_pca <- function(vcf){
  data('genomes_1k_df')
  data('ref_allele')
  ancestry_snps <- rownames(ref_allele)
  my_snps <- subset(vcf, paste0(str_extract(vcf@fix[,'CHROM'], '[0-9]+'), '_', vcf@fix[,'POS']) %in% as.character(row.names(ref_allele)))
  my_snps@fix[,'ID'] <- NA
  my_snps_gt <- extract.gt(my_snps, return.alleles = T)
  rownames(my_snps_gt) <- str_extract(rownames(my_snps_gt), '[^a-z]+')
  my_gt <- sapply(as.character(ref_allele$varid), function(x) ifelse(x %in% rownames(my_snps_gt), str_count(my_snps_gt[x,1], ref_allele[x, 'ref']) , 2))
  my_gt <- my_gt[colnames(genomes_1k_df[1:53])]
  snp_df <- rbind(genomes_1k_df, c(my_gt, 'Me', 'Me'))
  pca_df <- mutate_all(snp_df[,1:53], function(x) as.numeric(as.character(x)))
  pca <- prcomp(pca_df)
  plot_df <- data.frame(pca$x[,1:3], pop = snp_df$pop, superpop = snp_df$superpop)
  plot_df$pop <- factor(plot_df$pop, levels = c(sort(unique(genomes_1k_df$pop)), 'Me'))
  plot_df$superpop <- factor(plot_df$superpop, levels = c(sort(unique(genomes_1k_df$superpop)), 'Me'))
  
  ggplotly(
    ggplot(plot_df, aes(x = PC1, y= PC2, col = superpop, 
                               alpha = superpop, 
                               text = paste('Super-population:', superpop, 
                                            '\nPopulation:', pop))) + 
             geom_point(size = 0.9) +
             scale_color_manual(values = c("#ffff4d", "#FF5733", "#900C3F", "#43077B", "#66ccff", "black"),
                                name = "Population") +
             scale_alpha_manual(values = c(rep(0.5,5), 1), guide = F) +
             theme_classic(), tooltip = 'text')
}

get_ukbb_overlap <- function(vcf){
  data('ukb')
  my_vars <- paste0(str_extract(vcf@fix[,'CHROM'], '[0-9]+'), ':', vcf@fix[,'POS'], '_', vcf@fix[,'REF'], '_', vcf@fix[,'ALT'])
  overlap <- subset(ukb, ukb$varid %in% my_vars)
  return(overlap)
}

print_rare_vars <- function(ukbb_overlap_df, p_threshold, pal = pal1){
  rare_vars <- subset(ukbb_overlap_df, ukbb_overlap_df$percent_people < p_threshold)
  rare_vars$alt <- str_extract(rare_vars$varid, '[^_]+$')
  g_size <- floor(sqrt(nrow(rare_vars)))
  plot_df <- data.frame(x = rep(seq(1,g_size), g_size), y = rep(seq(g_size, 1, -1), each = g_size), var = as.factor(rare_vars$alt)[1:g_size**2])
  print(ggplot(plot_df, aes(x, y, fill = var)) +     
    geom_raster(interpolate = F, show.legend = F) + 
    scale_fill_manual(values = pal) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme_nothing())
}

hilbert_plot <- function(vcf, hil_pal, legend = NULL){
  hg19_seqlengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)
  names(hg19_seqlengths) <- paste0('chr', seq(1:22))
  sub_vcf <- subset(vcf, str_extract(vcf@fix[,'CHROM'], '[0-9]+') %in% seq(1:22))
  gr <- GRanges(seqnames = paste0('chr', str_extract(sub_vcf@fix[,'CHROM'], '[0-9]+')), ranges = IRanges(start = as.numeric(sub_vcf@fix[,'POS']), width = 1))
  seqlengths(gr) <- hg19_seqlengths[seqinfo(gr)@seqnames]
  mcols(gr)$density <- 100
  bins <- IRangesList(lapply(seqlengths(gr),
                             function(seqlen)
                               IRanges(breakInChunks(seqlen, 500))))
  avg_per_bin <- as(bins, "GRanges")
  seqinfo(avg_per_bin) <- seqinfo(gr)
  averageMCol <- function(colname){
    cvg <- coverage(gr, weight=colname)
    views_list <- RleViewsList(
      lapply(names(cvg),
             function(seqname)
               Views(cvg[[seqname]], bins[[seqname]])))
    unlist(viewMeans(views_list), use.names=FALSE)
  }
  mcols(avg_per_bin) <- DataFrame(lapply(mcols(gr)['density'], averageMCol))
  norm <- (avg_per_bin$density-min(avg_per_bin$density))/(max(avg_per_bin$density)-min(avg_per_bin$density))
  scale_size <- norm + 1.2
  norm <- norm * 99
  norm <- norm + 1
  norm <- round(norm)
  pal <- hil_pal(100)
  col <- pal[norm]
  hc = GenomicHilbertCurve(chr = paste0("chr", 1:22), level = 7, arrow = FALSE, background_col = 'black', padding = unit(1, 'mm'), legend = legend)
  hc_segments(hc, avg_per_bin, gp = gpar(lwd = scale_size, col = col))
}


SarahSeq <- function(){
  appDir <- system.file("app", package = "SarahSeq")
  shiny::runApp(appDir, display.mode = "normal")
}

readVCF <- function(filepath){
  vcf <- vcfR::read.vcfR(filepath)
  return(vcf)
}
