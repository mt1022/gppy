# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

- Python: `/Users/mt1022/miniforge3/envs/main/bin/python`
- uv: `/Users/mt1022/miniforge3/envs/main/bin/uv`
- **Never modify files outside this project directory.**

## Build & Test

```bash
# Build distribution packages
/Users/mt1022/miniforge3/envs/main/bin/uv pip install build && /Users/mt1022/miniforge3/envs/main/bin/python -m build

# Install from local source (editable)
/Users/mt1022/miniforge3/envs/main/bin/uv pip install -e .

# There is no formal test suite. Manual verification:
# 1. Run the Jupyter notebook in test/ against test data
# 2. Cross-validate txinfo output with the R script test/check_txinfo.R
# 3. Run CLI against the bundled test GTF files:
/Users/mt1022/miniforge3/envs/main/bin/python -m gppy txinfo -g test/human.chrY.gtf | head
/Users/mt1022/miniforge3/envs/main/bin/python -m gppy convert2bed -g test/human.chrY.gtf -t cds | head
```

## Architecture

This is a **single-module** package. All logic lives in `gppy/gtf.py` (~750 lines). There are no third-party dependencies — only Python stdlib.

### Core data model (3 classes)

- **`Region`** — dataclass with `start`/`end` (both 1-based, inclusive).
- **`Gene`** — dataclass: `gene_id`, `gene_name`, `chrom`, `strand`.
- **`Transcript`** — the workhorse. Holds lists of `Region` for `exons`, `cdss`, and `stop_codon`. All coordinate conversion methods live here:
  - `tpos_to_gpos` / `gpos_to_tpos` — single-position conversion between transcript and genomic coordinates.
  - `tiv_to_giv` — interval conversion (handles splicing by splitting across exons).
  - `cpos_to_gpos` / `gpos_to_cpos` — CDS-coordinate conversion (delegates to tpos/gpos).
  - `format_region_bed12` — serializes a list of Regions into BED12 format.
  - Properties: `introns`, `five_prime_utrs`, `three_prime_utrs`, `cds_len`, `cds_start`/`cds_end`, `tx_start`/`tx_end`, `stop_codon_start`/`stop_codon_end`.

### Data flow

1. **`parse_gtf(gtf_file, parse_attrs=False)`** — reads a GTF file (plain or `.gz`, or `-` for stdin), builds `dict[str, Transcript]`. When `parse_attrs=True`, also extracts per-transcript metadata (gene_name, biotype, tags, etc.) from the attributes column. Calls `Transcript.update()` on each transcript (sorts exons/CDSs/stop_codons by start position).

2. **Converter functions** call `parse_gtf()`, then iterate over an input TSV file to perform the requested coordinate conversion, printing results to stdout.

3. **`main()`** — argparse entry point with subcommands bound to each converter function.

### Subcommands

| Command | Function | Input |
|---------|----------|-------|
| `txinfo` | `tx_info()` | GTF → TSV with full metadata |
| `txinfo_basic` | `txinfo_basic()` | GTF → TSV (basic columns only, no attrs parsing) |
| `convert2bed` | `exon_to_bed` / `cds_to_bed` / `utr5_to_bed` / `utr3_to_bed` | GTF → BED12 |
| `t2g` | `t2g()` | tx_id + tpos → chrom + strand + gpos |
| `g2t` | `g2t()` | tx_id + gpos → tpos + region_type |
| `tiv2giv` | `tiv2giv()` | tx_id + tiv(start,end) → BED12 giv |
| `giv2tiv` | `giv2tiv()` | tx_id + giv(start,end) → tiv(start,end) + types |
| `extract_thick` | `extract_thick()` | BED12 → BED12 (nested thick regions) |

### Coordinate conventions (critical)

- **gpos, giv, tpos, tiv**: all **1-based, inclusive** (e.g., nucleotides 5–10 = `(5, 10)`).
- **BED output**: **0-based, half-open** (e.g., nucleotides 5–10 of chr1 = `chr1  4  10`).
- Strand handling: on `-` strand, transcript position 1 maps to the last nucleotide of the last exon. Coordinate conversion methods handle this internally — callers just provide positions in transcript space (1 to `len(transcript)`).

### `scripts/gtf_flybase.py`

A variant of the main module for FlyBase GTF files, which use a different attribute format than Ensembl. It is a standalone script (not importable) with the same CLI interface. Keep it in sync with `gppy/gtf.py` for shared logic changes.

### `test/` directory

Contains GTF test fixtures (human chrY/X/22, Drosophila BDGP6) and expected outputs (BED12, TSV). The notebook `test.ipynb` exercises the full pipeline. Test data files are large (some >100MB) — verify `.gitignore` excludes new generated outputs.

### Version alignment

`gppy/__init__.py` hardcodes `__version__ = '0.1.3'`, but `setup.cfg` declares `version = 0.1.4`. When bumping version, update both.
