R Notebook
================
Yman Mokrani

``` bash
wget https://github.com/ANF-MetaBioDiv/course-material/archive/refs/heads/main.zip
unzip main.zip
```

#### Cette condition permet de ne pas remplacer un dossier déjà existant.

``` r
refdb_folder <- here::here("data", "refdb")
refdb_folder
```

    ## [1] "/home/rstudio/ADM2023/data/refdb"

``` r
if (!dir.exists(refdb_folder)) dir.create(refdb_folder, recursive = TRUE)
```

``` bash
cp -R course-material-main/data/raw ./data
```

``` r
# so we change timeout to be 20 minutes
options(timeout = 1200)

# we save in variable the path to the refdb
# in the working space
silva_train_set <- file.path (refdb_folder, "silva_nr99_v138.1_train_set.fa.gz")

silva_species_assignment <- file.path(refdb_folder,
                                      "silva_species_assignment_v138.1.fa.gz")

# then we download the files if they don't already exist

if (!file.exists(silva_train_set)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
    silva_train_set,
    quiet = TRUE
  )
}

if (!file.exists(silva_species_assignment)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
    silva_species_assignment,
    quiet = TRUE
  )
}
```

#### Chercher toutes les fonctions dans r

``` r
devtools::load_all(path="/home/rstudio/ADM2023_tutoriel/course-material-main/R")
```

    ## ℹ Loading ANF_metaB

``` r
path_to_fastqs <- here::here("data", "raw")
```

#### le liste des fichiers qui contient tous

``` r
fnFs <- sort(list.files(path_to_fastqs,
                        pattern = "_R1.fastq.gz",
                        full.names = TRUE))
```

``` r
fnRs <- sort(list.files(path_to_fastqs,
                        pattern = "_R2.fastq.gz",
                        full.names = TRUE))
```

#### base name==\> nom du fichier

#### string split ==\> découpage de chaine de caractère avec comme séparateur le undersocre “*” Il y a dessus elements découpés par * dans la chaine de caractère

#### Si 3 valeurs alors 2 \_ pour séparer

#### sapply ==\> récupère le 1er element

``` r
sample_names <- basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1)
```

``` r
basename(fnFs) |>
  head()
```

    ## [1] "S11B_R1.fastq.gz" "S1B_R1.fastq.gz"  "S2B_R1.fastq.gz"  "S2S_R1.fastq.gz" 
    ## [5] "S3B_R1.fastq.gz"  "S3S_R1.fastq.gz"

``` r
basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1) |>
  head()
```

    ## [1] "S11B" "S1B"  "S2B"  "S2S"  "S3B"  "S3S"

### 3. Contrôle de qualité des séquences

``` r
# create a directory for the outputs
quality_folder <- here::here("outputs",
                             "dada2",
                             "quality_plots")

if (!dir.exists(quality_folder)) {
  dir.create(quality_folder, recursive = TRUE)
}

qualityprofile(fnFs,
               fnRs,
               file.path(quality_folder, "quality_plots.pdf"))
```

    ## png 
    ##   2

``` r
path_to_trimmed_reads <- here::here(
  "outputs",
  "dada2",
  "trimmed"
)

if (!dir.exists(path_to_trimmed_reads)) dir.create(path_to_trimmed_reads, recursive = TRUE)
```

``` r
primer_fwd  <- "CCTACGGGNBGCASCAG"
primer_rev  <- "GACTACNVGGGTATCTAAT"
```

``` r
Biostrings::readDNAStringSet(
  fnFs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   293 CCTACGGGGGGCAGCAGTAGGGA...ACATCGGCTTAACCGATGAAGT M01522:260:000000...
    ##  [2]   293 CCTACGGGTGGCACCAGTAGGGA...CGGGGCTTAACCTCGGAACTGC M01522:260:000000...
    ##  [3]   292 CCTACGGGGCGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGAGGT M01522:260:000000...
    ##  [4]   293 CCTACGGGGTGCAGCAGTAGGGA...TCAAAACTCCCAGTCTAGAGTT M01522:260:000000...
    ##  [5]   291 CCTACGGGTGGCAGCAGTGGGGA...GCAGTGGAAACTGTTGGGCTTG M01522:260:000000...
    ##  [6]   293 CCTACGGGATGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGGGGG M01522:260:000000...
    ##  [7]   292 CCTACGGGATGCAGCAGTGGGGA...TTTAATCCTGATGAGCTAGAAA M01522:260:000000...
    ##  [8]   293 CCTACGGGGCGCAGCAGTAGGGA...TTAAAACTTTTGTTCTGGAATT M01522:260:000000...
    ##  [9]   292 CCTACGGGTTGCAGCAGTGGGGA...ATTAAAACTTTTCAGCTAGAGT M01522:260:000000...
    ## [10]   293 CCTACGGGAGGCAGCAGTGGGGA...CCCGGGCTCAACCTGGGAACGG M01522:260:000000...

``` r
Biostrings::readDNAStringSet(
  fnRs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   301 GACTACCAGGGTATCTAATCCTG...GGCTGCTGGCACGAAGTTCGCC M01522:260:000000...
    ##  [2]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [3]   300 AATCCGGTTCGTGCCCCTAGGCT...TCTTTCCCAGCCCTTATTCCAA M01522:260:000000...
    ##  [4]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [5]   301 GACTACCGGGGTATCTAATCCCT...GGCTGCTGGCCCGGAATTAGCC M01522:260:000000...
    ##  [6]   301 GGTATCTAATCCGGTTCGTGCCC...CACCGTCCTTACCCCCCCCTTT M01522:260:000000...
    ##  [7]   301 GGTATCTAATCTTGTTTGCTCCC...CCCGACGTTAGCCGGGGCTTCT M01522:260:000000...
    ##  [8]   301 GACTACGAGGGTATCTAATCCCG...GGCTGCTGGCACGGAATTAGCC M01522:260:000000...
    ##  [9]   301 GGTATCTAATCCTCTTCGCTACC...CACGAAGTTAGCCGGACCTTCT M01522:260:000000...
    ## [10]   301 GACTACGGGGGTATCTAATCCTG...GGCTGCCGGCACGGGGTTAGCC M01522:260:000000...

``` bash
pwd
cp -R /home/rstudio/ADM2023_tutoriel/course-material-main/bash .
```

    ## /home/rstudio/ADM2023

``` r
(primer_log <- primer_trim(
  forward_files = fnFs,
  reverse_files = fnRs,
  primer_fwd = primer_fwd,
  primer_rev = primer_rev,
  output_dir = path_to_trimmed_reads,
  min_size = 200
))
```

    ##    sample status in_reads   in_bp too_short too_long too_many_n out_reads
    ## 1    S11B     OK     2000 1186767         0        0          0      1863
    ## 2     S1B     OK     2000 1186613         1        0          0      1855
    ## 3     S2B     OK     2000 1186942         0        0          0      1839
    ## 4     S2S     OK     2000 1186868         0        0          0      1833
    ## 5     S3B     OK     2000 1186650         0        0          0      1860
    ## 6     S3S     OK     2000 1186475         1        0          0      1880
    ## 7     S4B     OK     2000 1186331         2        0          0      1867
    ## 8     S4S     OK     2000 1186681         0        0          0      1872
    ## 9     S5B     OK     2000 1186386         1        0          0      1841
    ## 10    S5S     OK     2000 1186501         1        0          0      1861
    ## 11    S6B     OK     2000 1186261         2        0          0      1839
    ## 12    S6S     OK     2000 1187078         1        0          0      1835
    ## 13    S7B     OK     2000 1186888         0        0          0      1825
    ## 14    S7S     OK     2000 1186299         3        0          0      1845
    ## 15    S8B     OK     2000 1186354         3        0          0      1840
    ## 16    S8S     OK     2000 1186610         1        0          0      1848
    ## 17    S9B     OK     2000 1187038         0        0          0      1834
    ## 18    S9S     OK     2000 1186867         0        0          0      1835
    ##    w/adapters qualtrim_bp out_bp w/adapters2 qualtrim2_bp out2_bp
    ## 1        1986           0 513149        1876            0  528595
    ## 2        1975           0 511096        1877            0  525893
    ## 3        1987           0 506659        1850            0  521371
    ## 4        1989           0 504998        1843            0  519979
    ## 5        1989           0 512326        1870            0  527518
    ## 6        1989           0 517598        1891            0  532758
    ## 7        1980           0 514342        1884            0  529379
    ## 8        1987           0 515511        1884            0  530555
    ## 9        1984           0 506972        1856            0  522013
    ## 10       1991           0 512539        1869            0  527592
    ## 11       1981           0 506577        1857            0  521787
    ## 12       1982           0 505929        1851            0  520562
    ## 13       1987           0 503033        1836            0  517931
    ## 14       1987           0 508524        1857            0  523039
    ## 15       1993           0 507178        1847            0  522137
    ## 16       1982           0 509177        1865            0  524085
    ## 17       1983           0 505424        1851            0  520706
    ## 18       1979           0 505519        1853            0  520103

``` r
nopFw <- sort(list.files(path_to_trimmed_reads, pattern = "R1", full.names = TRUE))
nopRv <- sort(list.files(path_to_trimmed_reads, pattern = "R2", full.names = TRUE))
```

#### 5.

``` r
path_to_filtered_reads <- here::here("outputs", "dada2", "filtered")
if (!dir.exists(path_to_filtered_reads)) dir.create(path_to_filtered_reads, recursive = TRUE)
```

``` r
filtFs <- file.path(path_to_filtered_reads, basename(fnFs))
filtRs <- file.path(path_to_filtered_reads, basename(fnRs))
```

``` r
names(filtFs) <- sample_names
names(filtRs) <- sample_names
```

#### 5

``` r
(out <- dada2::filterAndTrim(
  fwd = nopFw,
  filt = filtFs,
  rev = nopRv,
  filt.rev = filtRs,
  minLen = 150,
  matchIDs = TRUE,
  maxN = 0,
  maxEE = c(3, 3),
  truncQ = 2
))
```

    ##                  reads.in reads.out
    ## S11B_R1.fastq.gz     1863      1200
    ## S1B_R1.fastq.gz      1855      1251
    ## S2B_R1.fastq.gz      1839      1255
    ## S2S_R1.fastq.gz      1833      1244
    ## S3B_R1.fastq.gz      1860      1244
    ## S3S_R1.fastq.gz      1880      1312
    ## S4B_R1.fastq.gz      1867      1262
    ## S4S_R1.fastq.gz      1872      1328
    ## S5B_R1.fastq.gz      1841      1255
    ## S5S_R1.fastq.gz      1861      1244
    ## S6B_R1.fastq.gz      1839      1251
    ## S6S_R1.fastq.gz      1835      1239
    ## S7B_R1.fastq.gz      1825      1203
    ## S7S_R1.fastq.gz      1845      1182
    ## S8B_R1.fastq.gz      1840      1169
    ## S8S_R1.fastq.gz      1848      1267
    ## S9B_R1.fastq.gz      1834      1195
    ## S9S_R1.fastq.gz      1835      1249

#### 6.1

``` r
errF <- dada2::learnErrors(filtFs,
                           randomize = TRUE,
                           multithread = TRUE)
```

    ## 6157072 total bases in 22350 reads from 18 samples will be used for learning the error rates.

``` r
errR <- dada2::learnErrors(filtRs,
                           randomize = TRUE,
                           multithread = TRUE)
```

    ## 6337638 total bases in 22350 reads from 18 samples will be used for learning the error rates.

``` r
dada2::plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](TD1_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### 6.2

#### on enleve les séquences identiques et on garde que les séquences représentatives.

``` r
derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S11B_R1.fastq.gz

    ## Encountered 754 unique sequences from 1200 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S1B_R1.fastq.gz

    ## Encountered 779 unique sequences from 1251 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S2B_R1.fastq.gz

    ## Encountered 789 unique sequences from 1255 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S2S_R1.fastq.gz

    ## Encountered 762 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S3B_R1.fastq.gz

    ## Encountered 772 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S3S_R1.fastq.gz

    ## Encountered 763 unique sequences from 1312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S4B_R1.fastq.gz

    ## Encountered 738 unique sequences from 1262 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S4S_R1.fastq.gz

    ## Encountered 638 unique sequences from 1328 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S5B_R1.fastq.gz

    ## Encountered 782 unique sequences from 1255 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S5S_R1.fastq.gz

    ## Encountered 663 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S6B_R1.fastq.gz

    ## Encountered 696 unique sequences from 1251 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S6S_R1.fastq.gz

    ## Encountered 657 unique sequences from 1239 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S7B_R1.fastq.gz

    ## Encountered 691 unique sequences from 1203 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S7S_R1.fastq.gz

    ## Encountered 675 unique sequences from 1182 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S8B_R1.fastq.gz

    ## Encountered 697 unique sequences from 1169 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S8S_R1.fastq.gz

    ## Encountered 714 unique sequences from 1267 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S9B_R1.fastq.gz

    ## Encountered 685 unique sequences from 1195 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S9S_R1.fastq.gz

    ## Encountered 677 unique sequences from 1249 total sequences read.

``` r
derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S11B_R2.fastq.gz

    ## Encountered 928 unique sequences from 1200 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S1B_R2.fastq.gz

    ## Encountered 948 unique sequences from 1251 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S2B_R2.fastq.gz

    ## Encountered 968 unique sequences from 1255 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S2S_R2.fastq.gz

    ## Encountered 925 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S3B_R2.fastq.gz

    ## Encountered 948 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S3S_R2.fastq.gz

    ## Encountered 967 unique sequences from 1312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S4B_R2.fastq.gz

    ## Encountered 953 unique sequences from 1262 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S4S_R2.fastq.gz

    ## Encountered 904 unique sequences from 1328 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S5B_R2.fastq.gz

    ## Encountered 975 unique sequences from 1255 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S5S_R2.fastq.gz

    ## Encountered 887 unique sequences from 1244 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S6B_R2.fastq.gz

    ## Encountered 914 unique sequences from 1251 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S6S_R2.fastq.gz

    ## Encountered 846 unique sequences from 1239 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S7B_R2.fastq.gz

    ## Encountered 881 unique sequences from 1203 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S7S_R2.fastq.gz

    ## Encountered 874 unique sequences from 1182 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S8B_R2.fastq.gz

    ## Encountered 879 unique sequences from 1169 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S8S_R2.fastq.gz

    ## Encountered 967 unique sequences from 1267 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S9B_R2.fastq.gz

    ## Encountered 892 unique sequences from 1195 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/ADM2023/outputs/dada2/filtered/S9S_R2.fastq.gz

    ## Encountered 911 unique sequences from 1249 total sequences read.

### 6.3

``` r
dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE)
```

    ## Sample 1 - 1200 reads in 754 unique sequences.
    ## Sample 2 - 1251 reads in 779 unique sequences.
    ## Sample 3 - 1255 reads in 789 unique sequences.
    ## Sample 4 - 1244 reads in 762 unique sequences.
    ## Sample 5 - 1244 reads in 772 unique sequences.
    ## Sample 6 - 1312 reads in 763 unique sequences.
    ## Sample 7 - 1262 reads in 738 unique sequences.
    ## Sample 8 - 1328 reads in 638 unique sequences.
    ## Sample 9 - 1255 reads in 782 unique sequences.
    ## Sample 10 - 1244 reads in 663 unique sequences.
    ## Sample 11 - 1251 reads in 696 unique sequences.
    ## Sample 12 - 1239 reads in 657 unique sequences.
    ## Sample 13 - 1203 reads in 691 unique sequences.
    ## Sample 14 - 1182 reads in 675 unique sequences.
    ## Sample 15 - 1169 reads in 697 unique sequences.
    ## Sample 16 - 1267 reads in 714 unique sequences.
    ## Sample 17 - 1195 reads in 685 unique sequences.
    ## Sample 18 - 1249 reads in 677 unique sequences.

``` r
dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE)
```

    ## Sample 1 - 1200 reads in 928 unique sequences.
    ## Sample 2 - 1251 reads in 948 unique sequences.
    ## Sample 3 - 1255 reads in 968 unique sequences.
    ## Sample 4 - 1244 reads in 925 unique sequences.
    ## Sample 5 - 1244 reads in 948 unique sequences.
    ## Sample 6 - 1312 reads in 967 unique sequences.
    ## Sample 7 - 1262 reads in 953 unique sequences.
    ## Sample 8 - 1328 reads in 904 unique sequences.
    ## Sample 9 - 1255 reads in 975 unique sequences.
    ## Sample 10 - 1244 reads in 887 unique sequences.
    ## Sample 11 - 1251 reads in 914 unique sequences.
    ## Sample 12 - 1239 reads in 846 unique sequences.
    ## Sample 13 - 1203 reads in 881 unique sequences.
    ## Sample 14 - 1182 reads in 874 unique sequences.
    ## Sample 15 - 1169 reads in 879 unique sequences.
    ## Sample 16 - 1267 reads in 967 unique sequences.
    ## Sample 17 - 1195 reads in 892 unique sequences.
    ## Sample 18 - 1249 reads in 911 unique sequences.

### 7.

``` r
mergers <- dada2::mergePairs(
  dadaF = dadaFs,
  derepF = derepFs,
  dadaR = dadaRs,
  derepR = derepRs,
  maxMismatch = 0,
  verbose = TRUE
)
```

    ## 879 paired-reads (in 28 unique pairings) successfully merged out of 970 (in 51 pairings) input.

    ## 835 paired-reads (in 33 unique pairings) successfully merged out of 943 (in 63 pairings) input.

    ## 783 paired-reads (in 30 unique pairings) successfully merged out of 944 (in 59 pairings) input.

    ## 929 paired-reads (in 32 unique pairings) successfully merged out of 1040 (in 59 pairings) input.

    ## 786 paired-reads (in 26 unique pairings) successfully merged out of 927 (in 60 pairings) input.

    ## 920 paired-reads (in 36 unique pairings) successfully merged out of 1040 (in 60 pairings) input.

    ## 808 paired-reads (in 29 unique pairings) successfully merged out of 971 (in 62 pairings) input.

    ## 1050 paired-reads (in 32 unique pairings) successfully merged out of 1130 (in 56 pairings) input.

    ## 905 paired-reads (in 24 unique pairings) successfully merged out of 1036 (in 40 pairings) input.

    ## 898 paired-reads (in 27 unique pairings) successfully merged out of 1039 (in 56 pairings) input.

    ## 970 paired-reads (in 31 unique pairings) successfully merged out of 1061 (in 51 pairings) input.

    ## 900 paired-reads (in 23 unique pairings) successfully merged out of 1062 (in 62 pairings) input.

    ## 823 paired-reads (in 31 unique pairings) successfully merged out of 988 (in 67 pairings) input.

    ## 852 paired-reads (in 30 unique pairings) successfully merged out of 968 (in 48 pairings) input.

    ## 842 paired-reads (in 26 unique pairings) successfully merged out of 944 (in 58 pairings) input.

    ## 849 paired-reads (in 31 unique pairings) successfully merged out of 1031 (in 62 pairings) input.

    ## 787 paired-reads (in 25 unique pairings) successfully merged out of 976 (in 55 pairings) input.

    ## 873 paired-reads (in 29 unique pairings) successfully merged out of 1044 (in 57 pairings) input.

#### assemblage ==\> contig ==\> plus de rid 1 et 2

#### si mismatch il enleve les séquences

### 8

``` r
seqtab <- dada2::makeSequenceTable(mergers)
```

### 9

``` r
seqtab_nochim <- dada2::removeBimeraDenovo(seqtab,
                                           method = "consensus",
                                           multithread = TRUE,
                                           verbose = TRUE)
```

    ## Identified 2 bimeras out of 162 input sequences.

### 10

``` r
taxonomy <- dada2::assignTaxonomy(
  seqs = seqtab_nochim,
  refFasta = silva_train_set,
  taxLevels = c("Kingdom", "Phylum", "Class",
                "Order", "Family", "Genus",
                "Species"),
  multithread = TRUE,
  minBoot = 60
)
```

``` r
taxonomy <- dada2::addSpecies(
  taxonomy,
  silva_species_assignment,
  allowMultiple = FALSE
)
```

### 11.1 EXPORT

``` r
export_folder <- here::here("outputs", "dada2", "asv_table")

if (!dir.exists(export_folder)) dir.create(export_folder, recursive = TRUE)

saveRDS(object = seqtab_nochim,
        file = file.path(export_folder, "seqtab_nochim.rds"))

saveRDS(object = taxonomy,
        file = file.path(export_folder, "taxonomy.rds"))
```

### 11.2

``` r
asv_seq <- colnames(seqtab_nochim)
```

``` r
ndigits <- nchar(length(asv_seq))
asv_id <- sprintf(paste0("ASV_%0", ndigits, "d"), seq_along(asv_seq))
```

``` r
row.names(taxonomy) <- colnames(seqtab_nochim) <- names(asv_seq) <- asv_id
```

``` r
taxonomy_export <- df_export(taxonomy, new_rn = "asv")

seqtab_nochim_export <- t(seqtab_nochim)
seqtab_nochim_export <- df_export(seqtab_nochim_export, new_rn = "asv")
```

``` r
write.table(taxonomy_export,
            file = file.path(export_folder, "taxonomy.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
```

``` r
write.table(seqtab_nochim_export,
            file = file.path(export_folder, "asv_table.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
```

#### Résultat asv tsv

``` r
Resultatasvtsv <- read.delim("outputs/dada2/asv_table/asv_table.tsv")
Resultatasvtsv
```

    ##         asv S11B S1B S2B S2S S3B S3S S4B S4S S5B S5S S6B S6S S7B S7S S8B S8S
    ## 1   ASV_001  111  60  45 104  60  78  13  66 180 178 155  70  48  79  52  67
    ## 2   ASV_002  110  22  52  78  48   0  46  41  62  36  60  32  57  49  60  91
    ## 3   ASV_003   24   0   0  84   9   0   0   7  16   0  90 295 108  43  60  10
    ## 4   ASV_004   88  55   0  96  49  39  47  61   0  51  81  33   0   0  55  61
    ## 5   ASV_005   75   0   0  12   0  26  31   0  70  60  28  23  80  64  77  55
    ## 6   ASV_006   58   0   0  35   0   0   0   0  56  36  73  47  26  75  59  55
    ## 7   ASV_007   32   0   0  74   0   0  10  10  14   0  56  73  78  36  66  24
    ## 8   ASV_008   39  31  54  42  67   0  45   0   0   0  33  19  30  45  45  34
    ## 9   ASV_009    0  57  39   0  50  45  42  47  94  38   0   0   0  48   0   0
    ## 10  ASV_010    0  58  60   0  67  50  43 150   0  62   0   0   0   0   0   0
    ## 11  ASV_011    0  26   0  46   0  35   0  30  31  28  44  25  37  27  30  40
    ## 12  ASV_012    0  50  32   0  59  96  59 120   0   0   0   0   0   0   0   0
    ## 13  ASV_013   24   0   0   0   0   0   0   0   0 100  69   0   0 104   0   0
    ## 14  ASV_014    0  29   0   0   0  92  52  94  34  31   0   0   0   0   0  22
    ## 15  ASV_015    0  56   0   0  49  54   0  70  54  29   0   0   0   0   0   0
    ## 16  ASV_016   20   0   0  61   0   0   0   0   0   0  22  18  67  19  47   0
    ## 17  ASV_017    0  43  57   0  59  42   0  63   0   0   0   0   0   0   0   0
    ## 18  ASV_018    0  60  64   0   0  20   0  76   0   0   0   0   0   0   0   0
    ## 19  ASV_019   42   0   0  10   0   0   0   0   0   0   8  21  11  14  22  24
    ## 20  ASV_020   58   0   0  37   0   0   0   0   0   0   0   0  23  25   0  20
    ## 21  ASV_021    6   0   0   2   0   6  22   0  19   9  13   9  16  13  21  13
    ## 22  ASV_022   25   0  16  20   0   0   0   0   0   0  32  13   0  21  13  18
    ## 23  ASV_023    0   0  32   0  25  27  24  31   0   0   0   0   0   0   0   0
    ## 24  ASV_024    0   0   0   0   0   0   0   0   0   0   0   0   0   0  77  96
    ## 25  ASV_025   20   0  27  24  33   0   0   0   0   0   0   3  16   0   0  23
    ## 26  ASV_026    0   0   0   0   0   0   0   0   0   9   0  82  33   0   0   0
    ## 27  ASV_027    0   0   0   0   0   0  38   0   0   0   0   0  25   0  34  23
    ## 28  ASV_028    0   0  27   0  26  34  31  19   0   0   0   0   0   0   0   0
    ## 29  ASV_029    0   0   0  15   0   0   0   0  34  33   0   0   9  10  17  24
    ## 30  ASV_030    0   0   0   0   0   0   0   0  84  40   0   0   0   0   0   0
    ## 31  ASV_031   13   5   0   4   5  53   5  11   0  12   0   0   3   0   0   0
    ## 32  ASV_032    0  32  14   0  20   7  13  23   0   0   0   0   0   0   0   0
    ## 33  ASV_033    0  28   9   0  10  31   0  19   6   0   0   0   5   0   0   0
    ## 34  ASV_034   16   0   0  29   0   0   0   0  20   0  20   0   0   0   0  17
    ## 35  ASV_035    0   0  67   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 36  ASV_036   23   0   0  14   0  15   0   6   0   0  18   0   0   5   0  11
    ## 37  ASV_037    0  46   0   0   0   0   0   0   0   0  16   0   0  18  14   0
    ## 38  ASV_038    0   0   0   0   0   0  89   0   0   0   0   0   0   0   0   0
    ## 39  ASV_039    9   0  20  16   0   0   0   0   0   0   0   0  14  16   0   0
    ## 40  ASV_040   34   0   0   0   0   0   0   0   0   0   0  22   0  30   0   0
    ## 41  ASV_041    0   0   0   0   0   0   0   0   0  19  18   0   0  22   0   0
    ## 42  ASV_042    0   0   0  39   0   0   0   0   0   0   0  40   0   0   0   0
    ## 43  ASV_043    0   0   0   0   0   0   0   0   0   0  11   0  12  11  11  18
    ## 44  ASV_044    0   0   0   0   0   0   0   0  28  32   0   0   0   0   0   0
    ## 45  ASV_045    0   0   0  11   0   0   0   0   0   0   0   0  16  14   0  21
    ## 46  ASV_046    0   0   0   0   0   0  55   5   0   0   0   0   0   0   0   0
    ## 47  ASV_047    0   0   0   0   0   0  30   0   0   0   0   0   0   0   0   0
    ## 48  ASV_048    0   0  19   0  14  13   9   0   0   0   0   0   0   0   0   0
    ## 49  ASV_049    0   0   0   6   0   0   0   0   7   7   0   0   5   9   8   7
    ## 50  ASV_050    0   0   0   0  49   0   0   0   0   0   0   0   0   0   0   0
    ## 51  ASV_051    0  11  28   0   8   0   0   0   0   0   0   0   0   0   0   0
    ## 52  ASV_052    6   5   0   0   0   0   0   0   0   0   6   0  14   0   8   7
    ## 53  ASV_053    0   0   0   0   0   0   0   0   0   0  24   0   0   0   0   0
    ## 54  ASV_054   14   0   0  16   0   0   0   0   0   0   5   0   0   0   0   5
    ## 55  ASV_055    0  11  11   0   0  18   0   0   0   0   0   0   0   0   0   0
    ## 56  ASV_056    0   0   0   0   0   0   0   0  16  15   0   0   0   0   0   0
    ## 57  ASV_057    0   9  14   0   0   0   0  11   4   0   0   0   0   0   0   0
    ## 58  ASV_058    0   0   0   0   0   0   0   0   0   0   9   0   0   0   9  13
    ## 59  ASV_059    0   0  14   3  15   0   0   4   0   0   0   0   0   0   0   0
    ## 60  ASV_060    0   0   0   0   0   0   0   0   0   0   0   0   0   0  35   0
    ## 61  ASV_061    0   0   0   0  13  11   0  10   0   0   0   0   0   0   0   0
    ## 62  ASV_062    0   0   0  10   0   0   0   0   0   0   0   0   0  11   0  12
    ## 63  ASV_063    0   0   0   0   0   0  19   0   0   0   0   0   0   0   0   0
    ## 64  ASV_064    0   0   0   0   0   0   0   0   0   0   0   0  31   0   0   0
    ## 65  ASV_065    0   0   0   0   0   0   0   0   0   0  30   0   0   0   0   0
    ## 66  ASV_066    0  27   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 67  ASV_067    0   0   0  10   0   0   0   0   0   0   0   0   5   6   0   0
    ## 68  ASV_068    0   0   0   0   0   0   0   0   0   0   0  27   0   0   0   0
    ## 69  ASV_069    0  11   0   0  15   0   0   0   0   0   0   0   0   0   0   0
    ## 70  ASV_070    0   0   0   0   0   0  26   0   0   0   0   0   0   0   0   0
    ## 71  ASV_071    0   0   0   0   0   0   0   0   0   0   0   0   0   0  13  13
    ## 72  ASV_072    0   0  12   0   0  13   0   0   0   0   0   0   0   0   0   0
    ## 73  ASV_073    0   0   0   0   0   0   0   0  15   0   0   0   0   0   0   0
    ## 74  ASV_074    0   0   0   0   0   0   0   0   0   0  15   0   0  10   0   0
    ## 75  ASV_075    0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 76  ASV_076    0   0   0   0   0   0   0   0   0  24   0   0   0   0   0   0
    ## 77  ASV_077    0   0   0   0   0   0   0   0   0   0   0   0   0  16   0   0
    ## 78  ASV_078    0   0   7   0   0   5   0  10   0   0   0   0   0   0   0   0
    ## 79  ASV_079    0   0   0   0   8  14   0   0   0   0   0   0   0   0   0   0
    ## 80  ASV_080    0   0   0   0   0   0   0  20   0   0   0   0   0   0   0   0
    ## 81  ASV_081    5   5   0   5   0   0   0   0   0   0   0   0   0   0   0   0
    ## 82  ASV_082    0  19   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 83  ASV_083    0   0   0   7   0   0   0   0   0   0   0   0   0   0   0  12
    ## 84  ASV_084    0   0   0   0   0   0   0   0  12   7   0   0   0   0   0   0
    ## 85  ASV_085   13   0   0   0   0   0   0   0   0   0   0   0   0   0   5   0
    ## 86  ASV_086    0   0   0   0  12   0   0   0   0   0   0   0   0   0   0   0
    ## 87  ASV_087    0   0   0   0   0   0   0   0  18   0   0   0   0   0   0   0
    ## 88  ASV_088    0   7   0   0   0   0   0   0   0   0   0   0  10   0   0   0
    ## 89  ASV_089    0   0   0   0   0   0  17   0   0   0   0   0   0   0   0   0
    ## 90  ASV_090    5   0   0   0   0   0   0   0   0   0  11   0   0   0   0   0
    ## 91  ASV_091    0  16   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 92  ASV_092    0   0  11   0   0   0   0   0   0   0   0   0   5   0   0   0
    ## 93  ASV_093    0   0   8   0   0   0   0   8   0   0   0   0   0   0   0   0
    ## 94  ASV_094    0   0   0   0   0   0   0   0  16   0   0   0   0   0   0   0
    ## 95  ASV_095    0   0   0   0   0   0   0   0   0  16   0   0   0   0   0   0
    ## 96  ASV_096    0   0   0   0   0   0   0   0   0   0   0   0  16   0   0   0
    ## 97  ASV_097    0   0  15   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 98  ASV_098    0   0   0   0   0  15   0   0   0   0   0   0   0   0   0   0
    ## 99  ASV_099    0   0   0   0   0   6   0   9   0   0   0   0   0   0   0   0
    ## 100 ASV_100    0   0   0   0   0   0   0   0  15   0   0   0   0   0   0   0
    ## 101 ASV_101    0   0   0   0   0   0   0   0   0  15   0   0   0   0   0   0
    ## 102 ASV_102    0   0   0   0   0  14   0   0   0   0   0   0   0   0   0   0
    ## 103 ASV_103    0   0   0   0   0  13   0   0   0   0   0   0   0   0   0   0
    ## 104 ASV_104    0   0   0   0   0   0  13   0   0   0   0   0   0   0   0   0
    ## 105 ASV_105    0   0   0   0   0   0   0   0   0   0   0  13   0   0   0   0
    ## 106 ASV_106    0   0   7   0   0   0   0   5   0   0   0   0   0   0   0   0
    ## 107 ASV_107    0   6   0   0   0   5   0   0   0   0   0   0   0   0   0   0
    ## 108 ASV_108    0   0   0   0  11   0   0   0   0   0   0   0   0   0   0   0
    ## 109 ASV_109    0   0   0   0   0   0   0   0   0   0  11   0   0   0   0   0
    ## 110 ASV_110    0   0   0   0   0   0   0   0   0   0   0  11   0   0   0   0
    ## 111 ASV_111    0   0   0   0   0   0   0   0   0   0   0   0   0   7   0   0
    ## 112 ASV_112    0   0   0   0   0   0   0   0   0   0   0   0   0   5   0   0
    ## 113 ASV_113    0  10   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 114 ASV_114    0   0   0   0   0   0  10   0   0   0   0   0   0   0   0   0
    ## 115 ASV_115    0   9   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 116 ASV_116    0   0   0   0   0   0   5   4   0   0   0   0   0   0   0   0
    ## 117 ASV_117    0   0   0   0   0   0   0   9   0   0   0   0   0   0   0   0
    ## 118 ASV_118    0   0   0   0   0   0   0   0   0   0   0   9   0   0   0   0
    ## 119 ASV_119    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 120 ASV_120    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 121 ASV_121    0   8   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 122 ASV_122    0   0   0   8   0   0   0   0   0   0   0   0   0   0   0   0
    ## 123 ASV_123    0   0   0   0   0   8   0   0   0   0   0   0   0   0   0   0
    ## 124 ASV_124    0   0   0   0   0   0   8   0   0   0   0   0   0   0   0   0
    ## 125 ASV_125    0   0   0   0   0   0   0   0   0   0   0   8   0   0   0   0
    ## 126 ASV_126    0   0   0   0   0   0   0   0   0   0   0   0   8   0   0   0
    ## 127 ASV_127    0   0   0   0   0   0   0   0   0   0   0   0   8   0   0   0
    ## 128 ASV_128    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   8
    ## 129 ASV_129    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 130 ASV_130    0   0   7   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 131 ASV_131    0   0   0   0   0   0   0   7   0   0   0   0   0   0   0   0
    ## 132 ASV_132    0   0   0   0   0   0   0   0   0   0   0   7   0   0   0   0
    ## 133 ASV_133    0   0   0   0   0   0   0   0   0   0   0   0   7   0   0   0
    ## 134 ASV_134    0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 135 ASV_135    0   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0
    ## 136 ASV_136    0   0   0   0   0   6   0   0   0   0   0   0   0   0   0   0
    ## 137 ASV_137    0   5   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 138 ASV_138    0   0   0   5   0   0   0   0   0   0   0   0   0   0   0   0
    ## 139 ASV_139    0   0   0   0   5   0   0   0   0   0   0   0   0   0   0   0
    ## 140 ASV_140    0   0   0   0   0   0   0   0   0   5   0   0   0   0   0   0
    ## 141 ASV_141    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   5
    ## 142 ASV_142    4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 143 ASV_143    0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0
    ## 144 ASV_144    0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0
    ## 145 ASV_145    0   0   0   0   0   0   4   0   0   0   0   0   0   0   0   0
    ## 146 ASV_146    0   0   0   0   0   0   0   4   0   0   0   0   0   0   0   0
    ## 147 ASV_147    0   0   0   0   0   0   0   0   0   4   0   0   0   0   0   0
    ## 148 ASV_148    0   0   0   0   0   0   0   0   0   0   4   0   0   0   0   0
    ## 149 ASV_149    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 150 ASV_150    3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 151 ASV_151    0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 152 ASV_152    0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0
    ## 153 ASV_153    0   0   0   0   0   0   0   0   0   0   3   0   0   0   0   0
    ## 154 ASV_154    0   0   0   0   0   0   0   0   0   0   3   0   0   0   0   0
    ## 155 ASV_155    2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## 156 ASV_156    0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0
    ## 157 ASV_157    0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0
    ## 158 ASV_158    0   0   0   0   0   0   0   0   0   0   2   0   0   0   0   0
    ## 159 ASV_159    0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0
    ## 160 ASV_160    0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0
    ##     S9B S9S
    ## 1    74  88
    ## 2    48  66
    ## 3    77  73
    ## 4    87   0
    ## 5    60  20
    ## 6    56  64
    ## 7    52  70
    ## 8    36  54
    ## 9     0  33
    ## 10    0   0
    ## 11   33  54
    ## 12    0   0
    ## 13    0  83
    ## 14    0   0
    ## 15    0   0
    ## 16   25  14
    ## 17    0   0
    ## 18    0   0
    ## 19   35  23
    ## 20   25   0
    ## 21   19   9
    ## 22   17   0
    ## 23    0  36
    ## 24    0   0
    ## 25    0  11
    ## 26   27   0
    ## 27    0  29
    ## 28    0   9
    ## 29    0   0
    ## 30    0   0
    ## 31    0   0
    ## 32    0   0
    ## 33    0   0
    ## 34    0   0
    ## 35   33   0
    ## 36    6   0
    ## 37    0   0
    ## 38    0   0
    ## 39   13   0
    ## 40    0   0
    ## 41    0  23
    ## 42    0   0
    ## 43   14   0
    ## 44    0  16
    ## 45    0   0
    ## 46    0   0
    ## 47    0  27
    ## 48    0   0
    ## 49    0   0
    ## 50    0   0
    ## 51    0   0
    ## 52    0   0
    ## 53   19   0
    ## 54    0   0
    ## 55    0   0
    ## 56    0   8
    ## 57    0   0
    ## 58    0   7
    ## 59    0   0
    ## 60    0   0
    ## 61    0   0
    ## 62    0   0
    ## 63    0  13
    ## 64    0   0
    ## 65    0   0
    ## 66    0   0
    ## 67    6   0
    ## 68    0   0
    ## 69    0   0
    ## 70    0   0
    ## 71    0   0
    ## 72    0   0
    ## 73    0  10
    ## 74    0   0
    ## 75    0   0
    ## 76    0   0
    ## 77    0   8
    ## 78    0   0
    ## 79    0   0
    ## 80    0   0
    ## 81    0   4
    ## 82    0   0
    ## 83    0   0
    ## 84    0   0
    ## 85    0   0
    ## 86    6   0
    ## 87    0   0
    ## 88    0   0
    ## 89    0   0
    ## 90    0   0
    ## 91    0   0
    ## 92    0   0
    ## 93    0   0
    ## 94    0   0
    ## 95    0   0
    ## 96    0   0
    ## 97    0   0
    ## 98    0   0
    ## 99    0   0
    ## 100   0   0
    ## 101   0   0
    ## 102   0   0
    ## 103   0   0
    ## 104   0   0
    ## 105   0   0
    ## 106   0   0
    ## 107   0   0
    ## 108   0   0
    ## 109   0   0
    ## 110   0   0
    ## 111   0   4
    ## 112   6   0
    ## 113   0   0
    ## 114   0   0
    ## 115   0   0
    ## 116   0   0
    ## 117   0   0
    ## 118   0   0
    ## 119   9   0
    ## 120   0   9
    ## 121   0   0
    ## 122   0   0
    ## 123   0   0
    ## 124   0   0
    ## 125   0   0
    ## 126   0   0
    ## 127   0   0
    ## 128   0   0
    ## 129   0   8
    ## 130   0   0
    ## 131   0   0
    ## 132   0   0
    ## 133   0   0
    ## 134   0   0
    ## 135   0   0
    ## 136   0   0
    ## 137   0   0
    ## 138   0   0
    ## 139   0   0
    ## 140   0   0
    ## 141   0   0
    ## 142   0   0
    ## 143   0   0
    ## 144   0   0
    ## 145   0   0
    ## 146   0   0
    ## 147   0   0
    ## 148   0   0
    ## 149   4   0
    ## 150   0   0
    ## 151   0   0
    ## 152   0   0
    ## 153   0   0
    ## 154   0   0
    ## 155   0   0
    ## 156   0   0
    ## 157   0   0
    ## 158   0   0
    ## 159   0   0
    ## 160   0   0

``` r
cat(paste0(">", names(asv_seq), "\n", asv_seq),
    sep = "\n",
    file = file.path(export_folder, "asv.fasta"))
```

#### Résultat asv.fasta

``` r
path_asv.fasta <- "outputs/dada2/asv_table/asv.fasta"
Resultatasvfasta <- readDNAStringSet(path_asv.fasta)
Resultatasvfasta
```

    ## DNAStringSet object of length 160:
    ##       width seq                                             names               
    ##   [1]   406 TGGGGAATTTTCCGCAATGGGC...GAAAGCCAGGGGAGCGAAAGGG ASV_001
    ##   [2]   404 TGGGGAATCTTGCACAATGGAG...GAAAGCATGGGTAGCGAAGAGG ASV_002
    ##   [3]   429 TGGGGAATATTGCACAATGGGC...GAAAGCGTGGGGAGCAAACAGG ASV_003
    ##   [4]   404 TGGGGAATCTTGCACAATGGAG...GAAAGCATGGGTAGCGAAGAGG ASV_004
    ##   [5]   387 GCGCGAAAACTTGACAATGCGA...GAAGCCTAGGGGCACGAACCGG ASV_005
    ##   ...   ... ...
    ## [156]   429 TGGGGAATATTGCGCAATGGGC...GAAAGCGTGGGGAGCAAACAGG ASV_156
    ## [157]   422 TTTCGAATCATTCACAATGGGC...GAAAGCATGGGGAGCGAAAGGG ASV_157
    ## [158]   404 TGGGGAATATTGGACAATGGGC...GAAAGCGTGGGGAGCAAACAGG ASV_158
    ## [159]   404 TAGGGAATATTGGACAATGGGG...GAAAGCATGGGTAGCGAAGAGG ASV_159
    ## [160]   404 TGGGGAATATTGGACAATGGGC...GAAAGCGTGGGTAGCAAACAGG ASV_160

### 11.3

``` r
getN <- function(x) sum(dada2::getUniques(x))

log_table <- data.frame(
  input = primer_log$in_reads,
  with_fwd_primer = primer_log$`w/adapters`,
  with_rev_primer = primer_log$`w/adapters2` ,
  with_both_primers = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim),
  perc_retained = rowSums(seqtab_nochim) / out[, 1] * 100
)

rownames(log_table) <- sample_names
```

``` r
df_export(log_table, new_rn = "sample") |>
  write.table(file = file.path(export_folder, "log_table.tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
```

#### Résultat log_table.tsv

``` r
Resultatlogtabletsv <- read.delim("outputs/dada2/asv_table/log_table.tsv")
Resultatlogtabletsv
```

    ##    sample input with_fwd_primer with_rev_primer with_both_primers filtered
    ## 1    S11B  2000            1986            1876              1863     1200
    ## 2     S1B  2000            1975            1877              1855     1251
    ## 3     S2B  2000            1987            1850              1839     1255
    ## 4     S2S  2000            1989            1843              1833     1244
    ## 5     S3B  2000            1989            1870              1860     1244
    ## 6     S3S  2000            1989            1891              1880     1312
    ## 7     S4B  2000            1980            1884              1867     1262
    ## 8     S4S  2000            1987            1884              1872     1328
    ## 9     S5B  2000            1984            1856              1841     1255
    ## 10    S5S  2000            1991            1869              1861     1244
    ## 11    S6B  2000            1981            1857              1839     1251
    ## 12    S6S  2000            1982            1851              1835     1239
    ## 13    S7B  2000            1987            1836              1825     1203
    ## 14    S7S  2000            1987            1857              1845     1182
    ## 15    S8B  2000            1993            1847              1840     1169
    ## 16    S8S  2000            1982            1865              1848     1267
    ## 17    S9B  2000            1983            1851              1834     1195
    ## 18    S9S  2000            1979            1853              1835     1249
    ##    denoisedF denoisedR merged nonchim perc_retained
    ## 1       1066      1015    879     879      47.18196
    ## 2       1004      1053    835     835      45.01348
    ## 3       1051      1023    783     783      42.57749
    ## 4       1093      1110    929     929      50.68194
    ## 5       1020      1001    786     786      42.25806
    ## 6       1108      1132    920     904      48.08511
    ## 7       1090      1032    808     808      43.27799
    ## 8       1181      1170   1050    1050      56.08974
    ## 9       1106      1099    905     905      49.15807
    ## 10      1107      1088    898     896      48.14616
    ## 11      1134      1108    970     970      52.74606
    ## 12      1094      1118    900     900      49.04632
    ## 13      1085      1032    823     823      45.09589
    ## 14      1038      1011    852     852      46.17886
    ## 15      1006      1000    842     842      45.76087
    ## 16      1113      1082    849     849      45.94156
    ## 17      1056      1022    787     787      42.91167
    ## 18      1101      1108    873     873      47.57493
