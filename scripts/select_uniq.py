#! /usr/bin/env python3

"""Cluster identical sequences.

CD-HIT is surprisingly inefficient for this particular task. Besides, it
does not allow to prefer non-putative proteins to be cluster representatives.
"""

import argparse
import gzip
import sys
from os.path import basename
from textwrap import wrap


def add_seq(nm, seq, seqs):
    if nm is None:
        return False
    seqs.setdefault(seq, set()).add(nm)
    return True


def load_seqs(infasta):
    seqs = dict()
    nm = None
    seq = []
    for line in infasta:
        if line.startswith(">"):
            add_seq(nm, "".join(seq), seqs)
            nm = line.split()[0].lstrip(">")
            seq = []
        else:
            seq.append(line.strip())
    add_seq(nm, "".join(seq), seqs)
    return seqs


def write_seqs(oufasta, seqs, width=80):
    for name, seq in sorted(seqs.items()):
        print(f">{name}", file=oufasta)
        print(*wrap(seq, width), sep="\n", file=oufasta)


def select_repr(nms, putative_tag="P"):
    reprs = [nm for nm in nms if not nm.endswith(putative_tag)] or nms
    return sorted(reprs, key=len)[0]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Cluster identical sequences.")
    parser.add_argument(
        "infasta", metavar="FASTA", help="""input .fasta file, may be gzipped,
        .gz extension is required in this case"""
    )
    parser.add_argument(
        "-o", "--output-tag", dest="outag", metavar="TAG", help="""output name tag,
        default is 'nr-' + input name without .fa* extension; the script will
        create <TAG>.fasta.gz with representative sequences, and <TAG>.clusters.tsv
        with clusters"""
    )
    parser.add_argument(
        "-c", "--cluster-prefix", dest="prefix", metavar="TAG", default="nr-",
        help="cluster name prefix, default is 'nr-'"
    )
    parser.add_argument(
        "-w", "--sequence-width", dest="width", metavar="N", type=int, default=80,
        help="output sequence width, default is 80"
    )
    parser.add_argument(
        "-n", "--column-name", dest="coln", metavar="NAME", default="Sequence_ID",
        help="first column header, default is 'Sequence_ID'"
    )
    args = parser.parse_args(argv)
    if args.infasta == "-":
        infasta = sys.stdin
    elif args.infasta.endswith(".gz"):
        infasta = gzip.open(args.infasta, "rt")
    else:
        infasta = open(args.infasta)
    with infasta:
        seqs = load_seqs(infasta)
    name = basename(infasta.name.strip("<>")).rsplit(".fa", 1)[0]
    outag = args.outag or f"nr-{name}"
    reprs = dict()
    clusters = []
    for seq, nms in seqs.items():
        repr_nm = select_repr(nms)
        clusters.append((-len(nms), repr_nm, len(seq), nms))
        reprs[repr_nm] = seq
    with gzip.open(f"{outag}.fasta.gz", "wt") as oufasta:
        write_seqs(oufasta, reprs, width=args.width)
    with open(f"{outag}.clusters.tsv", "w") as outsv:
        print(
            f"#:{args.coln}\tSequence_length\tRepresentative\tCluster_ID",
            file=outsv
        )
        for num, cluster in enumerate(sorted(clusters)):
            cluster_nm = f"{args.prefix}{num}"
            _invlen, repr_nm, seq_len, nms = cluster
            for nm in sorted(nms):
                print(nm, seq_len, repr_nm, cluster_nm, sep="\t", file=outsv)


if __name__ == "__main__":
    sys.exit(main())

