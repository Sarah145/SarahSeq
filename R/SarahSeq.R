pal1 <- c('#DAF7A6', '#FF5733', '#900C3F', '#43077B')
pal2 <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
pal3 <- c("#0D0887FF", "#9C179EFF", "#ED7953FF", "#F0F921FF")
pal4 <- c('#FFFF00', '#FF0000', '#00FF00', '#0000FF', '#000000')
pals <- c(list(pal1), list(pal2), list(pal3), list(pal4))

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
  #genomes_1k_df <- read.csv('~/Documents/Sarah-Seq/new_Sarah_Seq/www/ancestry_snps.csv', sep  = '\t', colClasses = c("pop" = "character", "superpop" = "character"))
  #ref_allele <- read.csv('~/Documents/Sarah-Seq/new_Sarah_Seq/www/ancestry_snps_ra.csv', sep  = '\t', colClasses = c("ref" = "character"))
  data('genomes_1k_df')
  data('ref_allele')
  ancestry_snps <- colnames(genomes_1k_df)[1:53]
  my_snps_gt <- extract.gt(subset(vcf, vcf@fix[,'ID'] %in% ancestry_snps), return.alleles = T)
  my_gt <- sapply(ancestry_snps, function(x) ifelse(x %in% rownames(my_snps_gt), str_count(my_snps_gt[x,1], ref_allele[x, 'ref']) , 2))
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
  ukb <- fread('ukbb_variants_small.tsv.gz')
  ukb$varid <- as.character(ukb$varid)
  my_vars <- paste0(str_sub(vcf@fix[,'CHROM'], 4, -1), ':', vcf@fix[,'POS'], '_', vcf@fix[,'REF'], '_', vcf@fix[,'ALT'])
  overlap <- subset(ukb, ukb$varid %in% my_vars)
  overlap$percent_people <- (overlap$n_non_ref/overlap$n_called) * 100
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
