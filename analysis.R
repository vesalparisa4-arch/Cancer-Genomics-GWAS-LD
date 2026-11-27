# Tutorial 2 - Cancer Genomics: GWAS and LD
# Computational Biomedicine WS 2026/27
# Run this script from any working directory. It searches for snps_10.dat
# in the data/ or tutorial_2_data/ folder next to the script and in the cwd.
# REQUIRES: R, BiocManager, snpStats and hexbin.
# The analysis follows the tutorial's specified snpStats::single.snp.tests()
# and snpStats::ld() procedures. Outputs are placed in "R_outputs".

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("snpStats", quietly=TRUE)) BiocManager::install("snpStats", ask=FALSE, update=FALSE)
if (!requireNamespace("hexbin", quietly=TRUE)) install.packages("hexbin")
suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(hexbin))

outdir <- "R_outputs"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

candidates <- c("data/snps_10.dat", "tutorial_2_data/snps_10.dat",
                "snps_10.dat", "../tutorial_2_data/snps_10.dat")
found <- candidates[file.exists(candidates)]
if (length(found) == 0) {
  stop("Cannot find snps_10.dat. Place it in 'data/' next to this script and set the working directory to the project folder.")
}
input <- found[1]
cat("Loading:", normalizePath(input), "\n")
loaded <- load(input)
cat("Loaded objects:", paste(loaded, collapse=", "), "\n")
print(ls())
stopifnot(all(c("snps.10", "subject.support", "snp.support") %in% loaded))

# (1) DATA INSPECTION
cat("\n=== DATA INSPECTION ===\n")
str(snps.10, max.level=1)
str(subject.support)
str(snp.support)
cat("Genotype dimensions (subjects x SNPs):", paste(dim(snps.10), collapse=" x "), "\n")
cat("Phenotype table columns:", paste(names(subject.support), collapse=", "), "\n")
cat("Annotation columns:", paste(names(snp.support), collapse=", "), "\n")
print(table(subject.support$cc, useNA="ifany"))
print(table(subject.support$stratum, useNA="ifany"))

# The SNP matrix rows, phenotype rows, and SNP annotation must be aligned.
stopifnot(identical(rownames(snps.10), rownames(subject.support)))
stopifnot(identical(colnames(snps.10), rownames(snp.support)))

# (2) GWAS
cat("\n=== HETEROZYGOSITY FILTER ===\n")
qc <- row.summary(snps.10)
heterozygosity <- qc$Heterozygosity  # fraction among called SNPs
keep <- is.finite(heterozygosity) & heterozygosity > 0
cat("Excluded subjects with zero/undefined heterozygosity:", sum(!keep), "\n")
cat("IDs excluded:", paste(rownames(snps.10)[!keep], collapse=", "), "\n")
cat("Remaining individuals:", sum(keep), "\n")
write.csv(data.frame(subject_id=rownames(snps.10),
                     heterozygosity=heterozygosity,
                     call_rate=qc$Call.rate,
                     retained=keep),
          file.path(outdir, "sample_qc.csv"), row.names=FALSE)

snps.filtered <- snps.10[keep, , drop=FALSE]
subject.filtered <- subject.support[keep, , drop=FALSE]
subject.filtered$stratum <- factor(subject.filtered$stratum)  # avoid numeric-stratum ambiguity
cat("Filtered case/control counts:\n")
print(table(subject.filtered$cc))

cat("\n=== SINGLE-SNP ASSOCIATION TESTS ===\n")
# Primary (unadjusted) 1-df Cochran-Armitage test for binary cc.
tests <- single.snp.tests(cc, data=subject.filtered,
                          snp.data=snps.filtered)
gwas_df <- as(tests, "data.frame")
stopifnot("P.1df" %in% colnames(gwas_df))

# Optional sensitivity analysis accounting for population group (stratum).
# This is NOT a replacement for the primary unadjusted exercise analysis.
tests.stratified <- single.snp.tests(cc, stratum=stratum,
                                     data=subject.filtered,
                                     snp.data=snps.filtered)
gwas_strat <- as(tests.stratified, "data.frame")

ann <- snp.support[rownames(gwas_df), , drop=FALSE]
stopifnot(identical(rownames(ann), rownames(gwas_df)))
result <- data.frame(
  snp_id=rownames(gwas_df),
  chromosome=ann$chromosome,
  position_bp=ann$position,
  N=gwas_df$N,
  chi_squared_1df=gwas_df[["Chi.squared.1.df"]],
  p_value_1df=gwas_df[["P.1df"]],
  p_value_stratified=gwas_strat[["P.1df"]],
  stringsAsFactors=FALSE)
valid <- is.finite(result$p_value_1df) & result$p_value_1df > 0
m <- sum(valid)  # Number of *actually valid* tests
alpha <- .05
bonferroni_p <- alpha/m
bonferroni_log <- -log10(bonferroni_p)
result$minus_log10_p <- -log10(result$p_value_1df)
result$bonferroni_significant <- valid & result$p_value_1df < bonferroni_p
write.csv(result, file.path(outdir,"gwas_results.csv"),row.names=FALSE)
cat("Tested SNPs:", ncol(snps.filtered), "; valid tests:",m,"\n")
cat("Bonferroni alpha/m:",format(bonferroni_p, digits=8),"\n")
cat("Bonferroni -log10(p) cutoff:",round(bonferroni_log,3),"\n")
cat("Significant SNPs:",sum(result$bonferroni_significant),"\n")
print(head(result[order(result$p_value_1df), ], 10))
write.csv(result[order(result$p_value_1df), ][1:20, ],
          file.path(outdir,"top20_gwas_snps.csv"),row.names=FALSE)

# Standard Manhattan scatter: appropriate because these data contain only chr10.
plotdf <- result[valid & is.finite(result$position_bp), ]
pdf(file.path(outdir,"manhattan_chr10.pdf"),width=11,height=5.5)
par(mar=c(4.5,5,3.2,1))
plot(plotdf$position_bp/1e6, plotdf$minus_log10_p,
     pch=16,cex=.44,col=adjustcolor("#227d89",alpha.f=.55),
     xlab="Chromosome 10 position (Mb)",ylab=expression(-log[10](p)),
     main="GWAS for case/control phenotype (cc) - chr10")
abline(h=bonferroni_log,col="#bc5f3c",lwd=2,lty=2)
legend("topright",legend=paste0("Bonferroni: p = ",signif(bonferroni_p,3)),
       col="#bc5f3c",lty=2,bty="n",cex=.82)
dev.off()

# Optional Manhattan HEXBIN figure: cells represent numbers of SNPs per bin.
# Uses hexbin::hexbin() plus base plotting so the significance line is exact.
hb <- hexbin(plotdf$position_bp/1e6,plotdf$minus_log10_p,xbins=105)
centers <- hcell2xy(hb)
bin_density <- log1p(hb@count)
bin_colors <- grDevices::colorRampPalette(c("#d6e9e9", "#24798a", "#0f2d4a"))(100)
ids <- 1L + floor(99*(bin_density-min(bin_density))/max(.Machine$double.eps, diff(range(bin_density))))
pdf(file.path(outdir,"manhattan_chr10_hexbin.pdf"),width=11,height=5.5)
par(mar=c(4.5,5,3.2,1))
plot(plotdf$position_bp/1e6,plotdf$minus_log10_p,type="n",
     xlab="Chromosome 10 position (Mb)",ylab=expression(-log[10](p)),
     main="GWAS - hexagonal SNP density (chr10)")
rx <- diff(range(plotdf$position_bp/1e6))/125
ry <- diff(range(plotdf$minus_log10_p))/70
angles <- seq(0,2*pi,length.out=7)+pi/6
for (i in seq_along(hb@count)) {
  polygon(centers$x[i]+rx*cos(angles),centers$y[i]+ry*sin(angles),
          col=bin_colors[ids[i]],border=NA)
}
abline(h=bonferroni_log,col="#bc5f3c",lwd=2,lty=2)
legend("topright",legend="Bonferroni threshold",col="#bc5f3c",lty=2,bty="n")
dev.off()

# (3) LINKAGE DISEQUILIBRIUM
cat("\n=== LINKAGE DISEQUILIBRIUM ===\n")
# 200 consecutive SNPs around the most significant GWAS marker.
lead_id <- result$snp_id[which.min(replace(result$p_value_1df,
                   !is.finite(result$p_value_1df),Inf))]
lead_index <- match(lead_id,colnames(snps.filtered))
region_start <- max(1L,min(lead_index-99L,ncol(snps.filtered)-199L))
region <- seq.int(region_start,length.out=200L)
cat("Lead SNP:",lead_id,"at",result$position_bp[match(lead_id,result$snp_id)],"bp\n")
cat("Selected region: SNP indices",min(region),"to",max(region),"\n")
cat("Region positions:",min(snp.support$position[region]),"to",
    max(snp.support$position[region]),"bp\n")

# The data use stratum 1 for CEU and stratum 2 for JPT/CHB.
# Verify this by the sample labels (ceu.*, jpt.*).
print(table(prefix=sub("\\..*$","",rownames(subject.filtered)),
            stratum=subject.filtered$stratum))

make_heatmap <- function(M,file_name,pop_name,metric,depth=100L) {
  z <- as.matrix(M)
  z[lower.tri(z,diag=FALSE)] <- NA_real_
  z[abs(row(z)-col(z))>depth] <- NA_real_
  diag(z) <- 1
  z[is.finite(z)] <- pmax(0,pmin(1,abs(z[is.finite(z)])))
  palette <- grDevices::colorRampPalette(c("#f3f7f8","#6cabb8","#124c68"))(150)
  pdf(file_name,width=8,height=7.4)
  par(mar=c(4.5,4.5,3.4,1))
  image(seq_len(nrow(z)), seq_len(ncol(z)), z,
        zlim=c(0,1),col=palette,useRaster=TRUE,asp=1,
        xlab="SNP index within selected 200",ylab="SNP index within selected 200",
        main=paste(pop_name,"-",metric,"(depth=100)"))
  box()
  mtext("White: uncalculated pairs / not available", side=1,line=3.1,cex=.7)
  dev.off()
  invisible(z)
}

ld_summary <- data.frame()
for (code in c(1,2)) {
  pop_name <- if (code==1) "CEU" else "JPT_CHB"
  snps.pop <- snps.filtered[as.character(subject.filtered$stratum)==as.character(code),region,drop=FALSE]
  stopifnot(ncol(snps.pop)==200)
  cat("Calculating LD:",pop_name," (n=",nrow(snps.pop),")\n")
  measures <- ld(snps.pop,stats=c("D.prime","R.squared"),depth=100)
  zD <- make_heatmap(measures$D.prime,
                     file.path(outdir,paste0("LD_",pop_name,"_Dprime.pdf")),
                     pop_name,"|D'|")
  zR <- make_heatmap(measures$R.squared,
                     file.path(outdir,paste0("LD_",pop_name,"_r2.pdf")),
                     pop_name,"r-squared")
  upper <- upper.tri(zD,diag=FALSE)
  for (metric in c("Dprime","r2")) {
    z <- if (metric=="Dprime") zD else zR
    keep_pairs <- upper & is.finite(z)
    near <- keep_pairs & (abs(row(z)-col(z))<=10)
    far <- keep_pairs & (abs(row(z)-col(z))>=50)
    ld_summary <- rbind(ld_summary,data.frame(
      population=pop_name,metric=metric,
      n_pairs=sum(keep_pairs),median_all=median(z[keep_pairs]),
      median_near_1_to_10=median(z[near]),
      median_far_50_to_100=median(z[far]),
      fraction_at_least_0_8=mean(z[keep_pairs]>=.8)))
  }
}
print(ld_summary)
write.csv(ld_summary,file.path(outdir,"ld_summary.csv"),row.names=FALSE)
cat("\nALL DONE. Files saved in:",normalizePath(outdir),"\n")
