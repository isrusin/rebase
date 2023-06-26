#! /usr/bin/env python3

"""Cluster identical sequences.

CD-HIT is surprisingly inefficient for this particular task. Besides, it
does not allow to select cluster representatives.
"""

import argparse
import gzip
import sys
from os.path import basename
from textwrap import wrap

from crc64iso.crc64iso import crc64


def load_nms(inlist):
    nms = set()
    for line in inlist:
        if line.startswith("#"):
            continue
        nms.add(line.strip())
    return nms


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


def select_repr(nms, preferred=set(), putative_tag="P"):
    reprs = set(nms).intersection(preferred)
    if not reprs:
        reprs = [nm for nm in nms if not nm.endswith(putative_tag)] or nms
    return select_shortest(reprs)


def select_shortest(nms):
    return sorted(nms, key=lambda nm: (len(nm), nm))[0]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Cluster identical sequences."
    )
    parser.add_argument(
        "infasta", metavar="FASTA", help="""input .fasta file, may be gzipped,
        .gz extension is required in this case"""
    )
    parser.add_argument(
        "-o", "--output-tag", dest="outag", metavar="TAG",
        help="""output name tag, default is 'nr-' + input name without .fa*
        extension; the script will create <TAG>.fasta.gz with representative
        sequences, and <TAG>.clusters.tsv with clusters"""
    )
    parser.add_argument(
        "-c", "--cluster-prefix", dest="prefix", metavar="TAG", default="nr-",
        help="cluster name prefix, default is 'nr-'"
    )
    parser.add_argument(
        "-w", "--sequence-width", dest="width", metavar="N", type=int,
        default=80, help="output sequence width, default is 80"
    )
    parser.add_argument(
        "-n", "--seq-title", dest="coln", metavar="STR", default="Sequence_ID",
        help="first column header, default is 'Sequence_ID'"
    )
    parser.add_argument(
        "-r", "--repr-title", dest="colr", default="Representative",
        metavar="STR", help="third column header, default is 'Representative'"
    )
    parser.add_argument(
        "--gold", dest="gold", metavar="FILE", type=argparse.FileType("r"),
        help="""list of protein which should be preferred as representative
        members"""
    )
    args = parser.parse_args(argv)
    preferred = set()
    if args.gold:
        with args.gold as inlist:
            preferred = load_nms(inlist)
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
        repr_nm = select_repr(nms, preferred)
        checksum = crc64(seq.upper())
        clusters.append((-len(nms), repr_nm, len(seq), nms, checksum))
        reprs[repr_nm] = seq
    with gzip.open(f"{outag}.fasta.gz", "wt") as oufasta:
        write_seqs(oufasta, reprs, width=args.width)
    with open(f"{outag}.clusters.tsv", "w") as outsv:
        print(
            f"#:{args.coln}\tSequence_length\t{args.colr}\tCluster_ID\tCRC64",
            file=outsv
        )
        for num, cluster in enumerate(sorted(clusters)):
            cluster_nm = f"{args.prefix}{num}"
            _invlen, repr_nm, seq_len, nms, checksum = cluster
            for nm in sorted(nms):
                print(
                    nm, seq_len, repr_nm, cluster_nm, checksum,
                    sep="\t", file=outsv
                )


if __name__ == "__main__":
    sys.exit(main())

