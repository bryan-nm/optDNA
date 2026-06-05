#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert a table of DNA or protein sequences (.csv, .tsv, or any delimited
text) into a fasta file. By default the first column holds names, the second
column holds sequences, and the first row is treated as a header. The output
fasta can be fed into the other tools in this project (e.g.
reverse_translate_proteins.py or optimize_DNA_for_assembly.py).
"""

import sys
import csv

def main(input_table, output_fasta, delimiter, name_col, seq_col, header):

    #Allow the delimiter to be given as an escape sequence on the command line
    #(e.g. -d '\t'), and accept the convenience aliases "tab" and "comma"
    delimiter = {"tab": "\t", "comma": ","}.get(delimiter.lower(),
                                                delimiter.encode().decode("unicode_escape"))

    records = []
    with open(input_table, newline = "") as table_in:
        reader = csv.reader(table_in, delimiter = delimiter)
        row_num = 0
        for row in reader:
            row_num += 1
            if header and row_num == 1:
                continue
            #Skip completely blank lines
            if len(row) == 0 or all(field.strip() == "" for field in row):
                continue
            #Make sure the requested columns actually exist in this row
            if len(row) < max(name_col, seq_col):
                sys.stderr.write("Warning: row %d has %d column(s), expected at least %d, skipping\n"
                                 % (row_num, len(row), max(name_col, seq_col)))
                continue
            name = row[name_col - 1].strip()
            seq = row[seq_col - 1].strip()
            if name == "" or seq == "":
                sys.stderr.write("Warning: row %d has an empty name or sequence, skipping\n" % row_num)
                continue
            records.append(">%s\n%s" % (name, seq))

    with open(output_fasta, "w+") as fastout:
        fastout.write("\n".join(records))

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--input',
          '-i',
          action = 'store',
          type = 'string',
          dest = 'input_table',
          help = "input table (.csv, .tsv, or other delimited text)")
    parser.add_option('--output',
          '-o',
          action = 'store',
          type = 'string',
          dest = 'output_fasta',
          help = "output fasta file")
    parser.add_option('--delimiter',
          '-d',
          action = 'store',
          type = 'string',
          dest = 'delimiter',
          help = "column delimiter; accepts escapes like '\\t' and the aliases 'tab'/'comma' (default: ',')",
          default = ",")
    parser.add_option('--name-col',
          '-n',
          action = 'store',
          type = 'int',
          dest = 'name_col',
          help = "1-based column number holding sequence names (default: 1)",
          default = 1)
    parser.add_option('--seq-col',
          '-s',
          action = 'store',
          type = 'int',
          dest = 'seq_col',
          help = "1-based column number holding sequences (default: 2)",
          default = 2)
    parser.add_option('--no-header',
          action = 'store_false',
          dest = 'header',
          help = "set if the table has no header row (default: assumes a header row)",
          default = True)

    (option, args) = parser.parse_args()

    main(option.input_table, option.output_fasta, option.delimiter,
         option.name_col, option.seq_col, option.header)
