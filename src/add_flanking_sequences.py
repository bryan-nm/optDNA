#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add flanking cloning sequences to the ends of each sequence in a fasta. A
head sequence is prepended and a tail sequence is appended to every entry.
e.g. --head ATCAGATA --tail GAGATCAG

The flanks are written exactly as provided (case preserved), so lowercase can
be used to visually mark the non-coding cloning sequence if desired.
"""

import sys

def main(input_fasta, output_fasta, head, tail):

    entries = dt.read_fasta_multi(input_fasta)

    flanked_entries = []
    for seq_name, seq in entries:
        if len(seq) == 0:
            sys.stderr.write("Warning: sequence '%s' is empty, skipping\n" % seq_name)
            continue
        flanked_entries.append((seq_name, head + seq + tail))

    with open(output_fasta, "w+") as fastout:
        fastout.write("\n".join(">%s\n%s" % (seq_name, seq)
                                for seq_name, seq in flanked_entries))

if __name__ == "__main__":

    import argparse
    import DNA_tools as dt

    parser = argparse.ArgumentParser(
        description = "Add flanking cloning sequences to the ends of each sequence in a fasta")
    parser.add_argument('--input', '-i',
          dest = 'input_fasta',
          help = "input fasta")
    parser.add_argument('--output', '-o',
          dest = 'output_fasta',
          help = "output fasta with flanking sequences added")
    parser.add_argument('--head',
          dest = 'head',
          default = "",
          help = "sequence to prepend to the 5' end of every entry (default: none)")
    parser.add_argument('--tail',
          dest = 'tail',
          default = "",
          help = "sequence to append to the 3' end of every entry (default: none)")

    option = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.head, option.tail)
