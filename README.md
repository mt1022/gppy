# Genomic Positioning with Python

`GPP` is a python package for genomic interval conversions to facilitate related transcriptome or translatome analysis. `GPP` can convert transcript/CDS intervals to genomic intervals in `bed12` format and vice versa, while taking well care of the presence of introns. Besides, `GPP` can extract mRNA/CDS/UTR from gtf and export in `bed12` format and generate summary table of basic transcript information (inlcuding ids and transcript/CDS/UTR lengths). More related features will be included in the future.

### Dependency and Installation
Scripts in this package rely only on the standard python (tested with version >= 3.7). No third party dependency is required. All the scripts can be run from the command line without installation after downloading.

```bash
wget https://raw.githubusercontent.com/mt1022/GPP/main/gpp/gtf.py
```

### Examples
Extract transcript length stats and metadata:
```bash
python gpp/gtf.py txinfo -g test/human.chrY.gtf >test/human.chrY.txinfo.tsv
cut -f1-9,12,15,19-22 test/human.chrY.txinfo.tsv | head
# tx_name	gene_id	chrom	strand	nexon	tx_len	cds_len	utr5_len	utr3_len	gene_name	transcript_biotype	ccds	ensembl_canonical	mane_select	basic
# ENST00000431340	ENSG00000215601	Y	+	4	443	0	0	0	TSPY24P	unprocessed_pseudogene	False	True	False	True
# ENST00000415010	ENSG00000215603	Y	-	1	1191	0	0	0	ZNF92P1Y	processed_pseudogene	False	True	False	True
# ENST00000449381	ENSG00000231436	Y	-	8	1145	0	0	0	RBMY3AP	unprocessed_pseudogene	False	True	False	True
# ENST00000436888	ENSG00000225878	Y	-	1	1164	0	0	0	SERBP1P2	processed_pseudogene	False	True	False	True
# ENST00000421279	ENSG00000236435	Y	-	5	868	0	0	0	TSPY12P	unprocessed_pseudogene	False	True	False	True
# ENST00000430032	ENSG00000278478	Y	+	1	279	0	0	0		processed_pseudogene	False	True	False	True
# ENST00000557448	ENSG00000258991	Y	+	1	1267	0	0	0	DUX4L19	unprocessed_pseudogene	False	True	False	True
# ENST00000651670	ENSG00000237048	Y	+	4	1123	0	0	0	TTTY12	lncRNA	False	False	False	True
# ENST00000413466	ENSG00000237048	Y	+	3	1046	0	0	0	TTTY12	lncRNA	False	True	False	False
```

### Usage
List utilities
```
$ usage: gtf.py [-h] {txinfo,convert2bed,t2g,g2t,tiv2giv,giv2tiv,extract_thick} ...

GTF file manipulation

options:
  -h, --help            show this help message and exit

GTF operations:
  {txinfo,convert2bed,t2g,g2t,tiv2giv,giv2tiv,extract_thick}
                        supported operations
    txinfo              summary information of each transcript
    convert2bed         convert GTF to bed12 format
    t2g                 convert tpos to gpos
    g2t                 convert gpos to tpos
    tiv2giv             convert tiv to giv
    giv2tiv             convert giv to tiv
    extract_thick       Extract nested thick regions from bed12
```

Extract basic transcript information
```
$ python gpp/gtf.py txinfo -h
usage: gtf.py txinfo [-h] [-g GTF]

options:
  -h, --help         show this help message and exit
  -g GTF, --gtf GTF  input gtf file (default: -)
```

Extract transcript/CDS/UTR features in GTF as bed12 format
```
$ python gpp/gtf.py convert2bed -h
usage: gtf.py convert2bed [-h] [-g GTF] [-t {exon,cds,utr5,utr3}] [-e EXTEND]

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     input gtf file (default: -)
  -t {exon,cds,utr5,utr3}, --type {exon,cds,utr5,utr3}
                        types of intervals to be converted to bed for each transcript (default: exon)
  -e EXTEND, --extend EXTEND
                        number of bases to extend at both sides (default: 0)
```

Convert transcript positions to genomic positions
```
$ python gpp/gtf.py t2g -h
usage: gtf.py t2g [-h] [-g GTF] [-i INFILE]

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     input gtf file (default: -)
  -i INFILE, --infile INFILE
                        tab-delimited file with the first two columns composed of tx_id and transcript coordinates (default: None)
```

Convert transcript intervals to genomic intervals (allow spliced regions)
```
$ python gpp/gtf.py tiv2giv -h
usage: gtf.py tiv2giv [-h] [-g GTF] [-i INFILE] [-a]

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     input gtf file (default: -)
  -i INFILE, --infile INFILE
                        tab-delimited file with the first three columns composed of tx_id, start and end coordinates (default: None)
  -a, --append          whether to append input at the end of the ouput (default: False)
```

Convert genomic positions to transcript positions
```
$ python gpp/gtf.py g2t -h
usage: gtf.py g2t [-h] [-g GTF] [-i INFILE]

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     input gtf file (default: -)
  -i INFILE, --infile INFILE
                        tab-delimited file with the first two columns composed of tx_id and genomic coordinates (default: None)
```

Convert genomic intervals to transcript intervals
```
$ python gpp/gtf.py giv2tiv -h
usage: gtf.py giv2tiv [-h] [-g GTF] [-i INFILE]

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     input gtf file (default: -)
  -i INFILE, --infile INFILE
                        tab-delimited file with the first three columns composed of tx_id, start and end coordinates (default: None)
```

### TODO
- [ ] File format conversion between gtf/bed/SAF;
- [x] improved transcrip info extraction in `txinfo` subcommand (for example, parse all possible tags);
- [x] check txinfo results compatilibality with previous R code relying on `GenomicFeatures`;

Please use the issues section to report if you have spotted any bug or want a feature to be implemented :)

