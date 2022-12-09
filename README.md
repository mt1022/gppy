# Genomic Positioning with Python

`GPP` is a python package for genomic interval conversions to facilitate related transcriptome or translatome analysis. `GPP` can convert transcript/CDS intervals to genomic intervals in bed12 format and vice versa, while taking well care of the presence of introns. Besides, `GPP` can extract mRNA/CDS/UTR from gtf and export in bed12 format and generate summary table of basic transcript information (inlcuding ids and lengths).

### Dependency and Installation
Scripts in this package rely only on the standard python (tested with version >= 3.7). No third party dependency is required. All the scripts can be run from the command line without installation after downloading.

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
- File format conversion between gtf/bed/SAF;
- improved transcrip info extraction in `txinfo` subcommand (for example, parse all possible tags);

