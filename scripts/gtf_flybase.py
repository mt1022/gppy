# -*- coding: utf-8 -*-
"""
Operations on genomic intervals stored in GTF file

note:
- all the exons of a gene should be on the same strand. Genes with exons trans-
  spliced from the other strand, like mod(mdg4) in D. melanogaster, should be
  excluded (before or after).
- stop codon is not part of the CDS, at least for Ensembl GTF
-----------------------------------------------------------------
@author: Hong Zhang (mt1022)
"""
import sys
import argparse
import re
import gzip
import csv
from dataclasses import dataclass
from collections import defaultdict
from typing import Any

###############################################################################
# class definiton #############################################################
###############################################################################

@dataclass
class Region:
    start: int
    end: int
    
    def __post_init__(self):
        if self.start > self.end:
            raise ValueError('Invalid region boundary!')

    def __len__(self):
        return self.end - self.start + 1


@dataclass
class Gene:
    gene_id  : str
    gene_name: str
    chrom    : str = ''
    strand   : str = '+'


class Transcript:
    def __init__(self, tx_id: str, gene: Gene):
        self.tx_id: str = tx_id
        self.gene: Gene = gene
        self.exons: list[Region] = []
        self.cdss: list[Region] = []
        self.stop_codon: list[Region] = []
    
    def add_region(self, region, region_type):
        if region_type == 'exon':
            self.exons.append(region)
        elif region_type == 'CDS':
            self.cdss.append(region)
        elif region_type == 'stop_codon':
            self.stop_codon.append(region)
        return
    
    def update(self):
        """
        Update the order of regions so that operations related to intervals can
        work correctly.
        """
        self.exons = sorted(self.exons, key=lambda r: r.start)
        self.cdss = sorted(self.cdss, key=lambda r: r.start)
        self.stop_codon = sorted(self.stop_codon, key=lambda r: r.start)
        return
    
    def __len__(self):
        return sum(len(i) for i in self.exons)
    
    @property
    def n_exons(self):
        return len(self.exons)

    @property
    def cds_len(self):
        return sum(len(i) for i in self.cdss)

    @property
    def tx_start(self):
        if self.gene.strand == '+':
            return self.exons[0].start
        else:
            return self.exons[-1].end

    @property
    def tx_end(self):
        if self.gene.strand == '+':
            return self.exons[-1].end
        else:
            return self.exons[0].start
    
    @property
    def cds_start(self):
        if len(self.cdss) == 0:
            return None
        elif self.gene.strand == '+':
            return self.cdss[0].start
        else:
            return self.cdss[-1].end
    
    @property
    def cds_end(self):
        if len(self.cdss) == 0:
            return None
        elif self.gene.strand == '+':
            return self.cdss[-1].end
        else:
            return self.cdss[0].start
    
    @property
    def stop_codon_start(self):
        if len(self.stop_codon) == 0:
            return None
        elif self.gene.strand == '+':
            return self.stop_codon[-1].end
        else:
            return self.stop_codon[0].start
    
    @property
    def stop_codon_end(self):
        if len(self.stop_codon) == 0:
            return None
        elif self.gene.strand == '+':
            return self.stop_codon[-1].end
        else:
            return self.stop_codon[0].start
    
    @property
    def introns(self):
        if len(self.exons) == 1:
            return []
        else:
            return [Region(self.exons[i].end + 1, self.exons[i+1].start - 1)
                    for i in range(self.n_exons - 1)]

    def tpos_to_gpos(self, pos: int):
        """
        transform transcript coordinate to genomic coordinate

        param pos: int, position on transcript, 1-based.
        """
        if pos < 1:
            return 0
        elif pos > len(self):
            return -1
        else:
            if self.gene.strand == '-':
                pos = len(self) - pos + 1
            for i in range(self.n_exons):
                if len(self.exons[i]) < pos:
                    pos -= len(self.exons[i])
                else:
                    return self.exons[i].start + pos - 1

    def gpos_to_tpos(self, pos):
        """
        transform genomic coordinate to transcript coordinate
        
        param pos: int, position on genome, 1-based.
        """
        if pos < self.exons[0].start:
            tpos = self.exons[0].start - pos
            ptype = 'upstream' if self.gene.strand == '+' else 'downstream'
            return tpos, ptype
        elif pos > self.exons[-1].end:
            tpos = pos - self.exons[-1].end
            ptype = 'downstream' if self.gene.strand == '+' else 'upstream'
            return tpos, ptype
        else:
            tpos = 0
            for i in range(self.n_exons):
                if self.exons[i].start <= pos:
                    if self.exons[i].end <= pos:
                        tpos += len(self.exons[i])
                    else:
                        tpos += pos - self.exons[i].start + 1
                else:
                    if self.exons[i-1].end < pos:
                        if self.gene.strand == '+':
                            ptype = 'intron_' + str(i)
                            tpos = pos - self.exons[i - 1].end
                        else:
                            ptype = 'intron_' + str(len(self.exons) - i)
                            tpos = self.exons[i].start - pos
                        return tpos, ptype
                    break
            ptype = 'exon'
            tpos = tpos if self.gene.strand == '+' else len(self) - tpos + 1
            return tpos, ptype
    
    def cpos_to_gpos(self, pos):
        """
        transform CDS coordinate to genomic coordinate

        param pos: int position on CDS, 1-based.
        """
        tpos = self.gpos_to_tpos(self.cds_start)[0] + pos - 1
        gpos = self.tpos_to_gpos(tpos)
        return gpos

    def gpos_to_cpos(self, pos):
        """
        transform genomic coordinate to CDS coordinate

        param: int, position on genome, 1-based.
        """
        tpos = self.gpos_to_tpos(pos)[0]
        cpos = tpos - self.gpos_to_tpos(self.cds_start)[0] + 1
        return cpos


    def tiv_to_giv(self, pos1, pos2):
        """
        given transcript region boundary:
        return one or more(for features spanning more than one exon)
        exonic region interval(s) in list of string interval

        param pos1: int, left transcript coordinate, 1-based.
        param pos2: int, right transcript coordinate, 1-based.
        """
        cod1 = self.tpos_to_gpos(pos1)
        cod2 = self.tpos_to_gpos(pos2)

        start = min(cod1, cod2)
        end = max(cod1, cod2)
        givs = []

        for i in range(self.n_exons):
            if self.exons[i].end < start:
                continue
            if self.exons[i].start > end:
                break
            if self.exons[i].start <= start:
                if self.exons[i].end <= end:
                    givs.append(Region(start, self.exons[i].end))
                else:
                    givs.append(Region(start, end))
            else:
                if self.exons[i].end <= end:
                    givs.append(Region(self.exons[i].start, self.exons[i].end))
                else:
                    givs.append(Region(self.exons[i].start, end))
        return givs
    
    @property
    def five_prime_utrs(self):
        if len(self.cdss) == 0 or self.cds_start == self.tx_start:
            return []
        else:
            return self.tiv_to_giv(1, self.gpos_to_tpos(self.cds_start)[0] - 1)
    
    @property
    def three_prime_utrs(self):
        if len(self.cdss) == 0 or self.stop_codon_end == self.tx_end or self.cds_end == self.tx_end:
            return []
        else:
            if len(self.stop_codon) > 0:
                return self.tiv_to_giv(self.gpos_to_tpos(self.stop_codon_end)[0] + 1, len(self))
            else:
                return self.tiv_to_giv(self.gpos_to_tpos(self.cds_end)[0] + 1, len(self))
    
    def format_region_bed12(self, rs, flank=0):
        """
        format a spliced region in a transcript into bed12 format

        param rs: a list of items of class Region
        """
        rs = sorted(rs, key=lambda r: r.start)
        if flank > 0:
            rs[0].start -= flank
            rs[-1].end += flank

        starts = [r.start - 1 for r in rs]
        ends = [r.end for r in rs]

        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(len(r)) for r in rs]

        s = [self.gene.chrom, starts[0], ends[-1], self.tx_id, self.gene.gene_id,
             self.gene.strand, '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']

        return s


###############################################################################
# functions ###################################################################
###############################################################################
def parse_gtf(gtf_file, parse_attrs=False):
    """
    read GTF file

    param: path to GTF file, gzipped format allowed.
    """
    gtf = {}
    if parse_attrs:
        tx_meta = {}
        metavars = dict(attrs=set(), tags=set())
        regex_attr = re.compile(r'(\w+) "?(.*?)"?;')

    if gtf_file.endswith('.gz'):
        f = gzip.open(gtf_file, 'rt')
    elif gtf_file == '-':
        f = sys.stdin
    else:
        f = open(gtf_file)

    regex_gid = re.compile(r'gene_id "(.*?)"')
    regex_tid = re.compile(r'transcript_id "(.*?)"')
    regex_name = re.compile(r'gene_symbol "(.*?)"')

    for line in f:
        if line[0] == '#':
            continue
        ary = line.strip().split('\t')
        m_tid = regex_tid.search(ary[8])
        m_gid = regex_gid.search(ary[8])
        if m_tid:
            tx_name = m_tid.group(1)
            if tx_name in gtf:
                gtf[tx_name].add_region(region = Region(int(ary[3]), int(ary[4])), region_type=ary[2])
            else:
                m_name = regex_name.search(ary[8])
                if m_name:
                    gene_name = m_name.group(1)
                else:
                    gene_name = m_gid.group(1)
                gene = Gene(gene_id=m_gid.group(1), gene_name=gene_name, chrom=ary[0], strand=ary[6])
                tx = Transcript(tx_id=tx_name, gene=gene)
                tx.add_region(region = Region(int(ary[3]), int(ary[4])), region_type=ary[2])
                gtf[tx_name] = tx
            if parse_attrs and ary[2] in ['mRNA', 'tRNA', 'miRNA', 'ncRNA', 'pre_miRNA', 'pseudogene', 'rRNA', 'snRNA', 'snoRNA']:
                attrs = dict()
                tags = dict()
                for m in regex_attr.finditer(ary[8]):
                    if m.group(1) == 'tag':
                        tags[m.group(2)] = True
                        metavars['tags'].add(m.group(2))
                    else:
                        attrs[m.group(1)] = m.group(2)
                        metavars['attrs'].add(m.group(1))
                tx_meta[tx_name] = dict(attrs=attrs, tags=tags)
            if parse_attrs and ary[2] == 'CDS':
                try:
                     _ = tx_meta[tx_name]['attrs']['protein_id']
                except KeyError:
                    for m in regex_attr.finditer(ary[8]):
                        if m.group(1) in ['protein_id', 'protein_version']:
                            metavars['attrs'].add(m.group(1))
                            tx_meta[tx_name]['attrs'][m.group(1)] = m.group(2)


    f.close()
    
    for tx in gtf:
        gtf[tx].update()
    if parse_attrs:
        return gtf, tx_meta, metavars
    else:
        return gtf


def exon_to_bed(gtf_file, extend=0):
    """
    print exons of each transcript in bed12 format

    param: path to GTF file, gzipped format allowed.
    """
    gtf = parse_gtf(gtf_file)
    for tx_id in gtf:
        tx = gtf[tx_id]
        items = tx.format_region_bed12(tx.exons, flank=extend)
        print('\t'.join(str(i) for i in items))
    return


def cds_to_bed(gtf_file, extend=0):
    """
    print CDSs of each transcript in bed12 format

    param: path to GTF file, gzipped format allowed.
    """
    gtf = parse_gtf(gtf_file)
    for tx_id in gtf:
        tx = gtf[tx_id]
        if len(tx.cdss) > 0:
            items = tx.format_region_bed12(tx.cdss, flank=extend)
            print('\t'.join(str(i) for i in items))
    return


def utr5_to_bed(gtf_file, extend=0):
    """
    print UTR5 of each transcript in bed12 format

    param: path to GTF file, gzipped format allowed.
    """
    gtf = parse_gtf(gtf_file)
    for tx_id in gtf:
        tx = gtf[tx_id]
        tx_utr5 = tx.five_prime_utrs
        if len(tx_utr5) > 0:
            items = tx.format_region_bed12(tx_utr5)
            print('\t'.join(str(i) for i in items))
    return


def utr3_to_bed(gtf_file, extend=0):
    """
    print UTR3 of each transcript in bed12 format

    param: path to GTF file, gzipped format allowed.
    """
    gtf = parse_gtf(gtf_file)
    for tx_id in gtf:
        tx = gtf[tx_id]
        tx_utr3 = tx.three_prime_utrs
        if len(tx_utr3) > 0:
            items = tx.format_region_bed12(tx_utr3, flank=extend)
            print('\t'.join(str(i) for i in items))
    return


def t2g(gtf_file, tfile):
    """
    convert transcript coordinates to genomic coordinates
    
    param: path to GTF file, gzipped format allowed.
    param tfile: tab-delimited file, 1st column=tx, 2nd column = tpos
    """
    gtf = parse_gtf(gtf_file)
    with open(tfile) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            try:
                tx = gtf[row[0]]
                gpos = tx.tpos_to_gpos(int(row[1]))
                row += [tx.gene.chrom, tx.gene.strand, str(gpos)]
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
                row += ['NA'] * 3
            print('\t'.join(row))
    return


def g2t(gtf_file, gfile):
    """
    convert genomic coordinates ot transcript coordinates

    param: path to GTF file, gzipped format allowed.
    param gfile: tab-delimited file, 1st column=tx, 2nd column = gpos
    """
    gtf = parse_gtf(gtf_file)
    with open(gfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                tpos, ptype = tx.gpos_to_tpos(int(row[1]))
                row += [str(tpos), ptype]
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
                row += ['NA'] * 2
            print('\t'.join(row))
    return


def tiv2giv(gtf_file, tivfile, append=False):
    """
    convert transcript intervals to genomic intervals

    param: path to GTF file, gzipped format allowed.
    param tivfile: tab-delimited, first three columns are tx_id, start, and end, 1-based
    """
    gtf = parse_gtf(gtf_file)
    with open(tivfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                givs = tx.tiv_to_giv(int(row[1]), int(row[2]))
                if append:
                    print('\t'.join([str(i) for i in tx.format_region_bed12(givs)] + row))
                else:
                    print('\t'.join(str(i) for i in tx.format_region_bed12(givs)))
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
            except IndexError:
                print('Cannot parse the entry:{}'.format(row), file=sys.stderr)
    return


def giv2tiv(gtf_file, givfile):
    """
    convert genomic intervals to transcript intervals

    param: path to GTF file, gzipped format allowed.
    param givfile: tab-delimited, first three columns are tx_id, start, and end, 1-based
    """
    gtf = parse_gtf(gtf_file)
    with open(givfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                if tx.gene.strand == '+':
                    tiv_l = list(tx.gpos_to_tpos(int(row[1])))
                    tiv_r = list(tx.gpos_to_tpos(int(row[2])))
                else:
                    tiv_l = list(tx.gpos_to_tpos(int(row[2])))
                    tiv_r = list(tx.gpos_to_tpos(int(row[1])))
                tiv = [str(tiv_l[0]), str(tiv_r[0]), tiv_l[1], tiv_r[1]]
                row += tiv
            except KeyError:
                row += ['NA'] * 4
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
            print('\t'.join(row))
    return


def tx_info(gtf_file, force_end_tags=False):
    """
    print summary information of each transcript

    param: path to GTF file, gzipped format allowed.
    note: stop codon is counted for CDS length, so that cds + utr5 + utr3 = transcript length
    """
    gtf, tx_meta, metavars= parse_gtf(gtf_file, parse_attrs=True)
    # support user provide list of attrs or tags?
    # attrs = attrs_user if attrs_user else list(metavars['attrs'].difference(['transcript_id', 'gene_id']))
    # tags = tags_user if tags_user else list(metavars['tags'])
    attrs = list(metavars['attrs'].difference(['transcript_id', 'gene_id']))
    tags = list(metavars['tags'])

    attrs.sort()
    tags.sort()
    if force_end_tags:
        tags += ['cds_start_NF', 'cds_end_NF', 'mRNA_start_NF', 'mRNA_end_NF']

    header = (['tx_name', 'gene_id', 'chrom', 'strand', 'nexon', 'tx_len', 'cds_len',
        'utr5_len', 'utr3_len'] + [i.lower() for i in attrs] + [i.lower() for i in tags])
    print('\t'.join(header))
    for tx_id in gtf:
        tx = gtf[tx_id]
        meta = tx_meta[tx_id]
        nexon = len(tx.exons)
        tx_len = len(tx)
        cds_len = sum(len(i) for i in tx.cdss) + sum(len(i) for i in tx.stop_codon)
        utr5_len = sum(len(i) for i in tx.five_prime_utrs)
        utr3_len = (tx_len - cds_len - utr5_len) if cds_len > 0 else 0
        out = ([tx.tx_id, tx.gene.gene_id, tx.gene.chrom, tx.gene.strand] +
            [str(i) for i in [nexon, tx_len, cds_len, utr5_len, utr3_len]] +
            [meta['attrs'].get(i, '') for i in attrs] +
            [str(meta['tags'].get(i, False)) for i in tags])
        print('\t'.join(out))
    return


def txinfo_basic(gtf_file):
    """
    print summary information of each transcript

    param: path to GTF file, gzipped format allowed.
    note: stop codon is counted for CDS length, so that cds + utr5 + utr3 = transcript length
    """
    gtf = parse_gtf(gtf_file, parse_attrs=False)

    header = ['tx_name', 'gene_id', 'chrom', 'strand', 'nexon', 'tx_len',
              'cds_len', 'utr5_len', 'utr3_len']
    print('\t'.join(header))
    for tx_id in gtf:
        if tx_id == '':
            continue
        tx = gtf[tx_id]
        nexon = len(tx.exons)
        tx_len = len(tx)
        cds_len = sum(len(i) for i in tx.cdss) + sum(len(i) for i in tx.stop_codon)
        utr5_len = sum(len(i) for i in tx.five_prime_utrs)
        utr3_len = (tx_len - cds_len - utr5_len) if cds_len > 0 else 0
        out = ([tx.tx_id, tx.gene.gene_id, tx.gene.chrom, tx.gene.strand] +
            [str(i) for i in [nexon, tx_len, cds_len, utr5_len, utr3_len]])
        print('\t'.join(out))
    return


def gtf_from_bed12(bed12_file):
    """
    parse transcript information represented in bed12 format
    NB: this function is write when adding `extract_thick`, but has not been tested yet
    """
    gtf = {}
    with open(bed12_file, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            ary = line.strip().split('\t')
            gene = Gene(gene_id=ary[3], chrom=ary[0], strand=ary[5])
            tx = Transcript(tx_id=ary[3], gene=gene)
            gstart = int(ary[1]) + 1
            gend = int(ary[2])
            block_sizes = [int(i) for i in ary[10].split(',') if i != '']
            block_starts = [int(i) for i in ary[11].split(',') if i != '']
            for b in range(int(ary[9])):
                b_start = gstart + block_starts[b]
                b_end   = gstart + block_starts[b] + block_sizes[b] - 1
                tx.add_region(region=Region(b_start, b_end), region_type='exon')
            gtf[ary[3]] = tx
    return gtf


def extract_thick(bed12_file):
    """
    extract bed12 format of thick regions in original bed12 format
    """
    with open(bed12_file, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            ary = line.strip().split('\t')
            # construct transcript model
            gene = Gene(gene_id=ary[3], chrom=ary[0], strand=ary[5])
            tx = Transcript(tx_id=ary[3], gene=gene)
            gstart = int(ary[1]) + 1
            gend = int(ary[2])
            block_sizes = [int(i) for i in ary[10].split(',') if i != '']
            block_starts = [int(i) for i in ary[11].split(',') if i != '']
            for b in range(int(ary[9])):
                b_start = gstart + block_starts[b]
                b_end   = gstart + block_starts[b] + block_sizes[b] - 1
                tx.add_region(region=Region(b_start, b_end), region_type='exon')
            # giv to tiv
            thick_gstart = int(ary[6]) + 1
            thick_gend   = int(ary[7])
            t1 = tx.gpos_to_tpos(thick_gstart)[0]
            t2 = tx.gpos_to_tpos(thick_gend)[0]
            if t1 < t2:
                thick_tstart = t1
                thick_tend   = t2
            else:
                thick_tstart = t2
                thick_tend   = t1
            # tiv to giv
            givs = tx.tiv_to_giv(thick_tstart, thick_tend)
            print('\t'.join(str(i) for i in tx.format_region_bed12(givs)))
    return


def main():
    # parent parser that holds common argument
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-g', '--gtf',
            type=str, default='-', help='input gtf file')
    
    
    # main parser with subparsers
    parser = argparse.ArgumentParser(prog='python gtf_flybase.py',
        description='GTF file manipulation')
    subparsers = parser.add_subparsers(title='GTF operations',
        help='supported operations', dest='subcmd')

    parser_txinfo = subparsers.add_parser('txinfo',
        help='summary information of each transcript',
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_txinfo.add_argument('-f', '--force_end_check', action='store_true',
        help='whether to include completeness checking tags')
    
    parser_txinfo_basic = subparsers.add_parser('txinfo_basic',
        help='summary information of each transcript (basic info only)',
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_tobed = subparsers.add_parser('convert2bed',
        help='convert GTF to bed12 format', parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_tobed.add_argument('-t', '--type',
        type=str, default='exon',
        choices=['exon', 'cds', 'utr5', 'utr3'],
        help='types of intervals to be converted to bed for each transcript')
    parser_tobed.add_argument('-e', '--extend',
        type=int, default=0,
        help='number of bases to extend at both sides')
    
    parser_t2g = subparsers.add_parser('t2g',
        help='convert tpos to gpos', parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_t2g.add_argument('-i', '--infile', type = str,
        help='tab-delimited file with the first two columns composed of'
        'tx_id and transcript coordinates')

    parser_g2t = subparsers.add_parser('g2t',
        help='convert gpos to tpos', parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_g2t.add_argument('-i', '--infile', type = str,
        help='tab-delimited file with the first two columns composed of '
        'tx_id and genomic coordinates')
    
    parser_tiv2giv = subparsers.add_parser('tiv2giv',
        help='convert tiv to giv', parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_tiv2giv.add_argument('-i', '--infile', type = str,
        help='tab-delimited file with the first three columns composed of '
        'tx_id, start and end coordinates')
    parser_tiv2giv.add_argument('-a', '--append', action='store_true',
        help='whether to append input at the end of the ouput')

    parser_giv2tiv = subparsers.add_parser('giv2tiv',
        help='convert giv to tiv', parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_giv2tiv.add_argument('-i', '--infile', type = str,
        help='tab-delimited file with the first three columns composed of '
        'tx_id, start and end coordinates')

    # bed12 input
    parser_extract_thick = subparsers.add_parser('extract_thick',
        help='Extract nested thick regions from bed12',
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    args = parser.parse_args()
    if args.subcmd == 'convert2bed':
        if args.type == 'exon':
            exon_to_bed(args.gtf, args.extend)
        elif args.type == 'cds':
            cds_to_bed(args.gtf, args.extend)
        elif args.type == 'utr5':
            utr5_to_bed(args.gtf, args.extend)
        else:
            utr3_to_bed(args.gtf, args.extend)
    elif args.subcmd == 'txinfo':
        tx_info(args.gtf, args.force_end_check)
    elif args.subcmd == 'txinfo_basic':
        txinfo_basic(args.gtf)
    elif args.subcmd == 't2g':
        t2g(gtf_file=args.gtf, tfile=args.infile)
    elif args.subcmd == 'g2t':
        g2t(gtf_file=args.gtf, gfile=args.infile)
    elif args.subcmd == 'tiv2giv':
        tiv2giv(gtf_file=args.gtf, tivfile=args.infile, append=args.append)
    elif args.subcmd == 'giv2tiv':
        giv2tiv(gtf_file=args.gtf, givfile=args.infile)
    elif args.subcmd == 'extract_thick':
        extract_thick(bed12_file=args.gtf)
    return


if __name__ == "__main__":
    main()
