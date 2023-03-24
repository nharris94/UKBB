library(tidyverse)



## sample file
## @: the file contains at least two columns,
##         IID (individual ID) and MeanCoverage (sequencing depth per sample)
smpl = read.csv(file = paste0("tryptase.sample.tsv"),
                sep = "\t", header = T, check.names = F, stringsAsFactors = F)



## tryptase genotype
## @: concatenate all pipeline output.txt files into one file
##         no header in the file
##         columns DEPTH, RNAME are required
##         add sample ID to the first column
geno = read.csv(file = paste0("tryptase.genotype.tsv"),
                sep = "\t", header = F, check.names = F, stringsAsFactors = F)
colnames(geno) = c("IID", "QNAME", "COUNT", "DEPTH", "RNAME", "POS", "CIGAR", "DISTANCE", "MISMATCH")
geno = merge(x = geno, y = smpl %>% select(IID, MeanCoverage), by = "IID") %>% arrange(IID)



## tryptase cnv
geno = geno %>% mutate(DEPTH.NORM = DEPTH/MeanCoverage, AC.DEPTH.NORM = NA)
for (allele in unique(geno$RNAME)) {
  idx = (geno$RNAME == allele)
  th = mean(geno$DEPTH.NORM[idx])
  geno$AC.DEPTH.NORM[idx] = round(geno$DEPTH.NORM[idx] / th)
}
geno$AC.DEPTH.NORM[geno$AC.DEPTH.NORM == 0] = 1

cnv = geno %>%
  select(IID, RNAME, AC.DEPTH.NORM) %>%
  group_by(IID, RNAME) %>%
  summarise(AC.DEPTH.NORM = sum(AC.DEPTH.NORM, na.rm = T)) %>%
  as.data.frame() %>%
  pivot_wider(id_cols = IID, names_from = RNAME, values_from = AC.DEPTH.NORM,
              names_glue = "{RNAME}.{.value}")
cnv[is.na(cnv)] = 0

cnv = cnv %>% mutate(CN.ALPHA = ALT2.AC.DEPTH.NORM,
                     CN.BETA = AB1.AC.DEPTH.NORM + B2.AC.DEPTH.NORM,
                     CN.ALL = AB1.AC.DEPTH.NORM + ALT2.AC.DEPTH.NORM + B2.AC.DEPTH.NORM)
cnv = cnv %>% mutate(CNV = 1 * ( (CN.ALPHA>2 & CN.BETA>1) | (CN.ALPHA>1 & CN.BETA>2) ))

write.table(x = cnv, file = paste0("tryptase.cnv.tsv"),
            sep = "\t", quote = F, na = "-", row.names = F, col.names = T)

