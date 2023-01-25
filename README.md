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

Extract CDS regions of each protein-coding transcript and export in bed12 format
```bash
python gpp/gtf.py convert2bed -g test/human.chrY.gtf -t cds >test/human.chrY.cds.bed12
head test/human.chrY.cds.bed12
# Y	22501564	22514067	ENST00000303728	ENSG00000169789	+	0	0	0	3	69,116,256,	0,2644,12247,
# Y	22501564	22512665	ENST00000477123	ENSG00000169789	+	0	0	0	3	69,116,28,	0,2644,11073,
# Y	12709447	12859413	ENST00000651177	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,
# Y	12709447	12859413	ENST00000338981	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,
# Y	12847044	12859413	ENST00000453031	ENSG00000114374	+	0	0	0	5	109,157,213,96,135,	0,9295,9588,10521,12234,
# Y	22072325	22084839	ENST00000303804	ENSG00000169807	-	0	0	0	3	256,116,69,	0,9754,12445,
# Y	22073730	22084839	ENST00000472391	ENSG00000169807	-	0	0	0	3	28,116,69,	0,8349,11040,
# Y	20575871	20592343	ENST00000361365	ENSG00000198692	+	0	0	0	7	16,84,104,51,82,92,3,	0,3736,6718,8602,12152,13612,16469,
# Y	20575871	20592343	ENST00000382772	ENSG00000198692	+	0	0	0	6	16,84,104,82,92,3,	0,3736,6718,12152,13612,16469,
# Y	22992343	22992376	ENST00000602732	ENSG00000183753	+	0	0	0	1	33,	0,
```

Convert CDS regions in genome coordinates to transcriptome coordinates
```bash
awk -v OFS="\t" '{print $4, $2 + 1, $3, $6}' test/human.chrY.cds.bed12 >test/human.chrY.cds.giv.tsv
python gpp/gtf.py giv2tiv -g test/human.chrY.gtf -i test/human.chrY.cds.giv.tsv >test/human.chrY.cds.tiv.tsv

head test/human.chrY.cds.giv.tsv
# ENST00000303728	22501565	22514067	+
# ENST00000477123	22501565	22512665	+
# ENST00000651177	12709448	12859413	+
# ENST00000338981	12709448	12859413	+
# ENST00000453031	12847045	12859413	+
# ENST00000303804	22072326	22084839	-
# ENST00000472391	22073731	22084839	-
# ENST00000361365	20575872	20592343	+
# ENST00000382772	20575872	20592343	+
# ENST00000602732	22992344	22992376	+

head test/human.chrY.cds.tiv.tsv
# ENST00000303728	22501565	22514067	+	228	668	exon	exon
# ENST00000477123	22501565	22512665	+	228	440	exon	exon
# ENST00000651177	12709448	12859413	+	587	8251	exon	exon
# ENST00000338981	12709448	12859413	+	946	8610	exon	exon
# ENST00000453031	12847045	12859413	+	1	710	exon	exon
# ENST00000303804	22072326	22084839	-	228	668	exon	exon
# ENST00000472391	22073731	22084839	-	228	440	exon	exon
# ENST00000361365	20575872	20592343	+	97	528	exon	exon
# ENST00000382772	20575872	20592343	+	79	459	exon	exon
# ENST00000602732	22992344	22992376	+	527	559	exon	exon
```
Convert CDS regions in transcriptome coordinates to genome coordinates
```bash
cut -f1,5,6 test/human.chrY.cds.tiv.tsv >test/human.chrY.cds.tiv2.tsv
python gpp/gtf.py tiv2giv -g test/human.chrY.gtf -i test/human.chrY.cds.tiv2.tsv -a >test/human.chrY.cds.giv2.bed12

head test/human.chrY.cds.tiv2.tsv
# ENST00000303728	228	668
# ENST00000477123	228	440
# ENST00000651177	587	8251
# ENST00000338981	946	8610
# ENST00000453031	1	710
# ENST00000303804	228	668
# ENST00000472391	228	440
# ENST00000361365	97	528
# ENST00000382772	79	459
# ENST00000602732	527	559

head test/human.chrY.cds.giv2.bed12
# Y	22501564	22514067	ENST00000303728	ENSG00000169789	+	0	0	0	3	69,116,256,	0,2644,12247,	ENST00000303728	228	668
# Y	22501564	22512665	ENST00000477123	ENSG00000169789	+	0	0	0	3	69,116,28,	0,2644,11073,	ENST00000477123	228	440
# Y	12709447	12859413	ENST00000651177	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,	ENST00000651177	587	8251
# Y	12709447	12859413	ENST00000338981	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,	ENST00000338981	946	8610
# Y	12847044	12859413	ENST00000453031	ENSG00000114374	+	0	0	0	5	109,157,213,96,135,	0,9295,9588,10521,12234,	ENST00000453031	1	710
# Y	22072325	22084839	ENST00000303804	ENSG00000169807	-	0	0	0	3	256,116,69,	0,9754,12445,	ENST00000303804	228	668
# Y	22073730	22084839	ENST00000472391	ENSG00000169807	-	0	0	0	3	28,116,69,	0,8349,11040,	ENST00000472391	228	440
# Y	20575871	20592343	ENST00000361365	ENSG00000198692	+	0	0	0	7	16,84,104,51,82,92,3,	0,3736,6718,8602,12152,13612,16469,	ENST00000361365	97	528
# Y	20575871	20592343	ENST00000382772	ENSG00000198692	+	0	0	0	6	16,84,104,82,92,3,	0,3736,6718,12152,13612,16469,	ENST00000382772	79459
# Y	22992343	22992376	ENST00000602732	ENSG00000183753	+	0	0	0	1	33,	0,	ENST00000602732	527	559

# the above should be identical to the CDS regions we extracted from GTF with `convert2bed`
head test/human.chrY.cds.bed12
# Y	22501564	22514067	ENST00000303728	ENSG00000169789	+	0	0	0	3	69,116,256,	0,2644,12247,
# Y	22501564	22512665	ENST00000477123	ENSG00000169789	+	0	0	0	3	69,116,28,	0,2644,11073,
# Y	12709447	12859413	ENST00000651177	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,
# Y	12709447	12859413	ENST00000338981	ENSG00000114374	+	0	0	0	44	96,149,80,113,219,116,252,139,153,105,207,137,134,88,343,96,212,241,150,121,131,279,126,126,170,109,147,147,223,221,191,174,142,751,124,226,130,186,221,89,157,213,96,135,	0,11141,12660,15665,17127,26164,26550,26963,28709,30077,47744,49058,51036,61620,64135,66023,67201,68571,69158,70078,76755,77074,80959,82051,83584,100731,101224,102187,103382,106676,108972,124240,128463,130416,131562,132792,133616,136885,137485,137791,146892,147185,148118,149831,
# Y	12847044	12859413	ENST00000453031	ENSG00000114374	+	0	0	0	5	109,157,213,96,135,	0,9295,9588,10521,12234,
# Y	22072325	22084839	ENST00000303804	ENSG00000169807	-	0	0	0	3	256,116,69,	0,9754,12445,
# Y	22073730	22084839	ENST00000472391	ENSG00000169807	-	0	0	0	3	28,116,69,	0,8349,11040,
# Y	20575871	20592343	ENST00000361365	ENSG00000198692	+	0	0	0	7	16,84,104,51,82,92,3,	0,3736,6718,8602,12152,13612,16469,
# Y	20575871	20592343	ENST00000382772	ENSG00000198692	+	0	0	0	6	16,84,104,82,92,3,	0,3736,6718,12152,13612,16469,
# Y	22992343	22992376	ENST00000602732	ENSG00000183753	+	0	0	0	1	33,	0,
```

Converison between genomic and transcriptomic positions for individual sites
```bash
cut -f1,2 test/human.chrY.cds.tiv2.tsv >test/human.chrY.cds.start.tpos.tsv
python gpp/gtf.py t2g -g test/human.chrY.gtf -i test/human.chrY.cds.start.tpos.tsv >test/human.chrY.cds.start.gpos.tsv

head test/human.chrY.cds.start.tpos.tsv
# ENST00000303728	228
# ENST00000477123	228
# ENST00000651177	587
# ENST00000338981	946
# ENST00000453031	1
# ENST00000303804	228
# ENST00000472391	228
# ENST00000361365	97
# ENST00000382772	79
# ENST00000602732	527

head test/human.chrY.cds.start.gpos.tsv
# ENST00000303728	228	Y	+	22501565
# ENST00000477123	228	Y	+	22501565
# ENST00000651177	587	Y	+	12709448
# ENST00000338981	946	Y	+	12709448
# ENST00000453031	1	Y	+	12847045
# ENST00000303804	228	Y	-	22084839
# ENST00000472391	228	Y	-	22084839
# ENST00000361365	97	Y	+	20575872
# ENST00000382772	79	Y	+	20575872
# ENST00000602732	527	Y	+	22992344

cut -f1,5 test/human.chrY.cds.start.gpos.tsv >test/human.chrY.cds.start.gpos2.tsv
python gpp/gtf.py g2t -g test/human.chrY.gtf -i test/human.chrY.cds.start.gpos2.tsv >test/human.chrY.cds.start.tpos2.tsv

head test/human.chrY.cds.start.gpos2.tsv
# ENST00000303728	22501565
# ENST00000477123	22501565
# ENST00000651177	12709448
# ENST00000338981	12709448
# ENST00000453031	12847045
# ENST00000303804	22084839
# ENST00000472391	22084839
# ENST00000361365	20575872
# ENST00000382772	20575872
# ENST00000602732	22992344

head test/human.chrY.cds.start.tpos2.tsv
# ENST00000303728	22501565	228	exon
# ENST00000477123	22501565	228	exon
# ENST00000651177	12709448	587	exon
# ENST00000338981	12709448	946	exon
# ENST00000453031	12847045	1	exon
# ENST00000303804	22084839	228	exon
# ENST00000472391	22084839	228	exon
# ENST00000361365	20575872	97	exon
# ENST00000382772	20575872	79	exon
# ENST00000602732	22992344	527	exon
```


### Usage
List utilities
```
$ python gpp/gtf.py -h
usage: gtf.py [-h] {txinfo,convert2bed,t2g,g2t,tiv2giv,giv2tiv,extract_thick} ...

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
- [x] improved transcrip info extraction in `txinfo` subcommand (for example, parse all possible tags);
- [x] check txinfo results compatilibality with previous R code relying on `GenomicFeatures`;
- [ ] File format conversion between gtf/bed/SAF;

Please use the issues section to report if you have spotted any bug or want a feature to be implemented :)

