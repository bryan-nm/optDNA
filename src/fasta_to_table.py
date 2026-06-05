#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert a fasta file of DNA or protein sequences into a table (.csv by
default). The first column is "Name" and the second column is "Sequence".
This is the inverse of table_to_fasta.py.
"""

import csv

def main(input_fasta, output_table, delimiter):

    #Allow the delimiter to be given as an escape sequence on the command line
    #(e.g. -d '\t'), and accept the convenience aliases "tab" and "comma"
    delimiter = {"tab": "\t", "comma": ","}.get(delimiter.lower(),
                                                delimiter.encode().decode("unicode_escape"))

    entries = dt.read_fasta_multi(input_fasta)

    with open(output_table, "w", newline = "") as table_out:
        writer = csv.writer(table_out, delimiter = delimiter)
        writer.writerow(["Name", "Sequence"])
        for seq_name, seq in entries:
            writer.writerow([seq_name, seq])

if __name__ == "__main__":

    import argparse
    import DNA_tools as dt

    parser = argparse.ArgumentParser(
        description = "Convert a fasta file into a table (.csv by default)")
    parser.add_argument('--input', '-i',
          dest = 'input_fasta',
          help = "input fasta file")
    parser.add_argument('--output', '-o',
          dest = 'output_table',
          help = "output table file (.csv)")
    parser.add_argument('--delimiter', '-d',
          dest = 'delimiter',
          default = ",",
          help = "column delimiter; accepts escapes like '\\t' and the aliases 'tab'/'comma' (default: ',')")

    option = parser.parse_args()

    main(option.input_fasta, option.output_table, option.delimiter)
