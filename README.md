# REBASE loader

## Description

This is a [redo](https://redo.readthedocs.io) based project for loading and
parsing of [REBASE](http://rebase.neb.com) data.
It is a set of files generated and updated through redo recipes (\*.do files).

## Structure

 *  `scripts/` – contains data parsers used in redo recipes.
 *  `raw_data/` – contains raw data loaded from REBASE ftp server.
 *  `proteins.fasta.gz` – protein sequences from `protein_org_seqs.txt`.
 *  `proteins.tsv` – tabular data derived from `protein_org_seqs.txt`.
 *  `genes.fasta.gz` – DNA sequences from `dna_seqs.txt`.
 *  `genes.tsv` – tabular data derived from `dna_seqs.txt`.
 *  `gold.fasta.gz` – protein sequences from `protein_gold_seqs.txt`.
 *  `gold.tsv` – tabular data derived from `protein_gold_seqs.txt`.
 *  `proteins.rnm.ids` – a list of all REBASE protein names.
 *  `proteins_np.rnm.ids` – a list of REBASE names of non-putative proteins.
 *  `proteins_gs.rnm.ids` – a list of all protein names from the Gold Standard.
 *  `rid-rnm.dct` – `REBASE_ID` to `REBASE_name` dict.
 *  `rebnr.clusters.tsv` – clusters of REBASE proteins with identical.
    sequences; non-putatives are preferred as representatives.
 *  `rebnr.fasta.gz` – REBNR representative sequences, contains all
    unique REBASE sequences.
 *  `rebnr.nrnm.ids` – a list of all REBNR names.
 *  `rebnr_np.nrnm.ids` – a list of REBNR representative names for clusters
     with non-putative proteins.
 *  `rebnr_gs.nrnm.ids` – a list of REBNR representative names for clusters
     that contain proteins from the Gold Standard.
 *  `names/` – contains data derived from interpretation of REBASE names.
     *  `name_data.tsv` – the main file that contains all data used for
        interpretation of REBASE component names and derived corrected
        annotations.
     *  `complexes.dct` – a list of REBASE complexes, `Complex_name` to
        `REBASE_name` dict.
     *  `rnm-rcnm.dct` – `REBASE_name` to `Complex_name` dict.
     *  `systems.dct` – a list of REBASE systems, `System_name` to
        `REBASE_name` dict.
     *  `rnm-rsnm.dct` – `REBASE_name` to `System_name` dict.

## Requirements

 *  [redo](https://redo.readthedocs.io)
 *  [lab-tools](https://github.com/isrusin/lab-tools) (for names/)
     -  `sieve_tsv.py` (for names/)
     -  `transform_tsv.py` (for names/)
 *  python 3
     -  [crc64iso](https://pypi.org/project/crc64iso)
     -  [etsv](https://github.com/isrusin/etsv) (for names/)
 *  coreutils
     -  `cut`
     -  `paste`
     -  `sort`
     -  `touch`
 *  other
     -  `curl`
     -  `gzip`, `zcat`
     -  `sed`
