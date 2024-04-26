### Karyotype Plot:
1. Prepare the ROH data used to plot the ROH density
    
  a. merge the ROH files for each sample to one file "all.roh.bed"


```python
## add all ROH files for all samples to one file
for sample in sampleList
do 
    cat ${sample}.roh.bed >> all.roh.bed
done
```

  b. remove the centromeres from all.roh.bed


```python
## remove the centromeres from all.roh.bed using bedtools subtract
bedtools subtract -a all.roh.bed -b centromeres.bed > roh.bed
```

2. Extract the frequent and infrequent regions to plot 

  a. Create a bedgraph using the ROH regions for all samples


```python
## to get a file with the lengths of the chromosomes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom,size from hg38.chromInfo" | grep -v "_" | 
grep -v size | grep -v -i chrM | sort -Vf > hg38.genome

## be sure that your bed file is sorted based on the chrom and pos
sort -k 1,1 roh.bed > roh.sorted.bed

## produce a bedgraph output by dividing the regions based on the breakpoints in the samples and then counting how many times the region is found in the samples
bedtools genomecov -bga -i roh.sorted.bed -g hg38.genome > roh.merged.bed
```

  b. Get rare and frequent regions


```python
## written in R

filepath <- "roh.merged.bed"
merged <- read.table(filepath, header = TRUE, col.names = c("chrom", "start", "end", "count"))
merged <- merged[!(merged$chrom %in% c("chrX", "chrY")), ]

tmp <- quantile(merged$count, probs = c(0.05,0.95))
rare <- as.numeric(tmp[1]) # anything less than this is rare hom
freq <- as.numeric(tmp[2]) # anything more than this is freq hom

rare_regions <- merged[merged$count <= rare, ]
freq_regions <- merged[merged$count >= freq, ]

write.table(rare_regions, file = "rare_regions.bed", sep = "\t", row.names=FALSE, quote = FALSE)
write.table(freq_regions, file = "freq_regions.bed", sep = "\t", row.names=FALSE, quote = FALSE)
```

3. Extract deletion and duplicatin CNV variants
   - In order to get the regions we create two files ( one for DEL variants and one for DUP )
   - grep DEL first to get all DEL regions and then grep DEL & 1/1 to write the regions again for HOM variants
   - the same goes for DUP
   - the code for this is as follows:



```python
for SM in `cat sample-list`
do 

    grep -i DUP cnv/gff3/${SM}.cnv.gff3 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDUP.bed
    grep -i DUP cnv/gff3/${SM}.cnv.gff3 | grep -i 1/1 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDUP.bed
    
    grep -i DEL cnv/gff3/${SM}.cnv.gff3 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDEL.bed
    grep -i DEL cnv/gff3/${SM}.cnv.gff3 | grep -i 1/1 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDEL.bed
done  
```


```python
cat cnv-bedfiles/annotatedDEL.bed | cut -f 1,2,3 > cnv-bedfiles/DEL.bed
cat cnv-bedfiles/annotatedDUP.bed | cut -f 1,2,3 > cnv-bedfiles/DUP.bed
```

4. Extract deletion and insertion SV variants
   - In order to get the regions we create two files ( one for DEL variants and one for INS )
   - grep DEL first to get all DEL regions and then grep DEL & 1/1 to write the regions again for HOM variants
   - the same goes for INS
   - the code for this is as follows:


```python
## INS
for SM in `cat sample-list`
do 
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 1,2 > sv-bedfiles/INS1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/INS2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/INS3.bed
    paste -d '\t' INS1.bed INS2.bed INS3.bed >> annotatedINS.bed
    
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 1,2 > sv-bedfiles/INS1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/INS2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/INS3.bed
    paste -d '\t' INS1.bed INS2.bed INS3.bed >> annotatedINS.bed
done

## DEL
for SM in `cat sample-list`
do 
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 1,2 > sv-bedfiles/DEL1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/DEL2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/DEL3.bed
    paste -d '\t' DEL1.bed DEL2.bed DEL3.bed >> annotatedDEL.bed
    
    echo $SM >> annotatedDEL.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 1,2 > sv-bedfiles/DEL1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/DEL2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/DEL3.bed
    paste -d '\t' DEL1.bed DEL2.bed DEL3.bed >> annotatedDEL.bed
done
```


```python
cat sv-bedfiles/annotatedINS.bed | cut -f 1,2,3 > sv-bedfiles/INS.bed
cat sv-bedfiles/annotatedDEL.bed | cut -f 1,2,3 > sv-bedfiles/DEL.bed
```

5. Karyotype script in R


```python
suppressWarnings(library(karyoploteR, quietly = TRUE))

## ROH
filepath <- "roh.bed"
roh <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
roh <- toGRanges(roh)

## rare and freq regions
filepath <- "rare_regions.bed"
rare_regions <- read.table(filepath, header = TRUE)
roh_rare <- toGRanges(rare_regions)
roh_rare <- roh_rare[roh_rare$count > 0,]
roh_rare <- reduce(roh_rare)
roh_rare <- roh_rare[width(roh_rare)> 100000,]

filepath <- "freq_regions.bed" 
freq_regions <- read.table(filepath, header = TRUE)
roh_freq <- toGRanges(freq_regions)
roh_freq <- reduce(roh_freq)
roh_freq <- roh_freq[width(roh_freq)> 100000,]

### PLOT A - Autosomal ROH with rare and freq regions
pdf("../karyotype_plot_a.pdf", 10, 10) #width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")
kpPlotDensity(kp, roh, col = 'gray')
kpPlotRegions(kp, data=roh_rare, data.panel = 2, col = colors[[2]])
kpPlotRegions(kp, data=roh_freq, data.panel = 1, col = colors[[1]] )
legend(x = 0.7, y = 0.7, fill = c ('gray',colors[[2]], colors[[1]]), 
                        legend = c("ROH",'infrequent','frequent'), cex = 1)
dev.off()



## CNV Files
filepath <- "cnv-bedfiles/DEL.bed"
cnv_del <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
cnv_del <- toGRanges(cnv_del)

filepath <- "cnv-bedfiles/DUP.bed"
cnv_dup <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
cnv_dup <- toGRanges(cnv_dup)

### PLOT B - Autosomal CNV
pdf("../karyotype_plot_b.pdf", 10, 10)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")#, plot.params = pp)
kpPlotDensity(kp, cnv_dup, col = colors[[4]])
kpPlotDensity(kp, cnv_del,  data.panel = 2, col = colors[[5]])
legend(x = 0.7, y = 0.7, fill = c (colors[[4]],colors[[5]]), 
                        legend = c("CNV DUP","CNV DEL"), cex = 1)
dev.off()

### PLOT d - Sex CNV
pdf("../karyotype_plot_d.pdf", 7, 7)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = c("chrX", "chrY"))#, plot.params = pp)
kpPlotDensity(kp, cnv_dup, col = colors[[4]])
kpPlotDensity(kp, cnv_del,  data.panel = 2, col = colors[[5]])
legend(x = 0.7, y = 0.3, fill = c (colors[[4]],colors[[5]]), 
                        legend = c("CNV DUP","CNV DEL"), cex = 1)
dev.off()

## SV Files
filepath <- "sv-bedfiles/DEL.bed"
sv_del <- read.table(filepath, header = FALSE, fill = TRUE, col.names = c("chrom", "start", "end"))
sv_del <- toGRanges(sv_del)

## SV Files
filepath <- "INS.bed"
sv_ins <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
sv_ins <- toGRanges(sv_ins)

### PLOT C - Autosomal SV
pdf("../karyotype_plot_c.pdf", 10, 10)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")
kpPlotDensity(kp, sv_ins,  col = colors[[1]])
kpPlotDensity(kp, sv_del,  col = colors[[3]], data.panel=2)
legend(x = 0.7, y = 0.7, fill = c (colors[[1]],colors[[3]]), 
                        legend = c("SV INS","SV DEL"), cex = 1)
dev.off()

### PLOT e - Sex SV
pdf("../karyotype_plot_e.pdf", 7, 7)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = c("chrX", "chrY"))
kpPlotDensity(kp, sv_ins,  col = colors[[1]])
kpPlotDensity(kp, sv_del,  col = colors[[3]], data.panel=2)
legend(x = 0.7, y = 0.3, fill = c (colors[[1]],colors[[3]]), 
                        legend = c("SV INS","SV DEL"), cex = 1)
dev.off()
```
