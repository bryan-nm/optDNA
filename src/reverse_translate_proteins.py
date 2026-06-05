#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reverse translate a protein fasta into a DNA fasta using loaded codon
frequencies. No further optimization is performed here -- the intent is that
the output can be piped into optimize_DNA_for_assembly.py.
"""

def main(input_fasta, output_fasta, bias_file):

    entries = dt.read_fasta_multi(input_fasta)
    codon_freqs = dt.read_codon_freqs(bias_file)

    translated_entries = []
    for seq_name, pro_seq in entries:
        if len(pro_seq) == 0:
            sys.stderr.write("Warning: sequence '%s' is empty, skipping\n" % seq_name)
            continue
        try:
            DNA_seq = reverse_translate_protein(pro_seq, codon_freqs)
        except Exception as e:
            sys.stderr.write("Warning: failed to reverse translate sequence '%s' (%s), skipping\n"
                             % (seq_name, e))
            continue
        translated_entries.append((seq_name, DNA_seq))

    with open(output_fasta,'w+') as fastout:
        fastout.write('\n'.join(">%s\n%s" % (seq_name, DNA_seq)
                                for seq_name, DNA_seq in translated_entries))

def reverse_translate_protein(protein, codon_freqs):
    #The codon bias table only carries the 20 amino acids, so a trailing stop
    #(*) is handled separately by appending a stop codon, mirroring the
    #behavior in optimize_DNA_for_assembly.codon_optimize
    protein = protein.upper()
    if protein.endswith("*"):
        return dt.reverse_translate(protein.strip("*"), CB = codon_freqs) + "TAA"
    else:
        return dt.reverse_translate(protein, CB = codon_freqs)

if __name__ == "__main__":

    from optparse import OptionParser
    import DNA_tools as dt
    import sys

    parser = OptionParser()
    parser.add_option('--input',
          '-i',
          action = 'store',
          type = 'string',
          dest = 'input_fasta',
          help = "input fasta with your protein sequence(s) to be reverse translated")
    parser.add_option('--output',
          '-o',
          action = 'store',
          type = 'string',
          dest = 'output_fasta',
          help = "output fasta with your reverse translated DNA sequence(s)")
    parser.add_option('--bias',
          '-b',
          action = 'store',
          type = 'string',
          dest = 'bias_file',
          help = "codon bias file",
          default = "../codon_sets/ColiProteomeContent.tsv")

    (option, args) = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.bias_file)
