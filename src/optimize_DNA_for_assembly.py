#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:16:01 2021

@author: bryanandrews1
"""

def main(input_fasta, output_fasta, bias_file):

    entries = dt.read_fasta_multi(input_fasta)

    optimized_entries = []
    for seq_name, DNA_seq in entries:
        if len(DNA_seq) == 0:
            sys.stderr.write("Warning: sequence '%s' is empty, skipping\n" % seq_name)
            continue
        try:
            DNA_seq = codon_optimize(DNA_seq, codon_bias = bias_file)
            DNA_seq = fix_GC(DNA_seq)
            DNA_seq = strip_microhomology(DNA_seq)
            DNA_seq = strip_mononucleotide_tracts(DNA_seq)
            DNA_seq = strip_G_quadruplexes(DNA_seq)
        except Exception as e:
            sys.stderr.write("Warning: failed to optimize sequence '%s' (%s), skipping\n"
                             % (seq_name, e))
            continue
        optimized_entries.append((seq_name, DNA_seq))

    with open(output_fasta,'w+') as fastout:
        fastout.write('\n'.join(">%s\n%s" % (seq_name, DNA_seq)
                                for seq_name, DNA_seq in optimized_entries))

def strip_G_quadruplexes(DNA_seq):
    #Find any occurrence of 3 G's in a row and break it
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
    attempt_count = 0
    while True:
        attempt_count +=1
        if attempt_count > 100:
            break
        if DNA_seq.find("GGG") == -1:
            break
        else:
            pos = DNA_seq.find("GGG")
            if pos%3 == 0: #GGG
                codon_seq[pos//3] = random.choice(dt.codon_synonyms["GGG"])
            elif pos%3 == 1: #GGg
                codon1 = codon_seq[(pos-1)//3]
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos-1)//3] = random.choice(dt.codon_synonyms[codon1])
            elif pos%3 == 2: #Ggg
                codon1 = codon_seq[(pos-2)//3]
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos-2)//3] = random.choice(dt.codon_synonyms[codon1])
            #re-sync DNA_seq so the next .find() reflects the edits we just made
            DNA_seq = ''.join(codon_seq)
    return ''.join(codon_seq)
    
def strip_mononucleotide_tracts(DNA_seq):
    #Find any occurrence of 5+ of the same base in a row and break it
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
    attempt_count = 0
    while True:
        #keep going until you find all of them
        #if you find a 5-mer, you will have 2-3 codons to play with
        attempt_count +=1
        if attempt_count > 1000:
            break
        if   DNA_seq.find("AAAAA") != -1:
            pos = DNA_seq.find("AAAAA")
        elif DNA_seq.find("CCCCC") != -1:
            pos = DNA_seq.find("CCCCC")
        elif DNA_seq.find("GGGGG") != -1:
            pos = DNA_seq.find("GGGGG")
        elif DNA_seq.find("TTTTT") != -1:
            pos = DNA_seq.find("TTTTT")
        else:
            break
        
        #If you found a mononucleotide tract, find the codons and replace them
        #with randomly selected synonymous codons
        if pos%3 == 0: #aaaAA
            codon1 = codon_seq[pos//3]
            if len(dt.codon_synonyms[codon1]) > 0:
                codon_seq[pos//3] = random.choice(dt.codon_synonyms[codon1])
        elif pos%3 == 1: #aaAAA
            codon1 = codon_seq[(pos-1)//3]
            if len(dt.codon_synonyms[codon1]) > 0:
                codon_seq[(pos-1)//3] = random.choice(dt.codon_synonyms[codon1])
            codon2 = codon_seq[(pos+2)//3]
            if len(dt.codon_synonyms[codon2]) > 0:
                codon_seq[(pos+2)//3] = random.choice(dt.codon_synonyms[codon2])
        elif pos%3 == 2: #aAAAa
            codon1 = codon_seq[(pos-2)//3]
            if len(dt.codon_synonyms[codon1]) > 0:
                codon_seq[(pos-2)//3] = random.choice(dt.codon_synonyms[codon1])
            codon2 = codon_seq[(pos+1)//3]
            if len(dt.codon_synonyms[codon2]) > 0:
                    codon_seq[(pos+1)//3] = random.choice(dt.codon_synonyms[codon2])

        #re-sync DNA_seq so the next .find() reflects the edits we just made
        DNA_seq = ''.join(codon_seq)

    return ''.join(codon_seq)

def strip_microhomology(DNA_seq):
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
    #make a library of all 8mers in the sequence
    #An 8mer is the smallest kmer where you have at least two codons to play with
    #(Actually, you're guaranteed to have exactly two codons to play with)
    kmer_lib = {}
    for i in range(len(DNA_seq)-7):
        kmer = DNA_seq[i:i+8]
        if kmer_lib.get(kmer) == None:
            kmer_lib[kmer] = 0
        kmer_lib[kmer] +=1
        
    #Find all the cases where the 8mer occurs more than once
    for kmer in kmer_lib:
        if kmer_lib[kmer] == 1:
            continue
        elif kmer_lib[kmer] > 1:
            #go to position of first occurence
            pos = DNA_seq.find(kmer)
            if pos%3 == 0: #nnnNNNnn
                codon1 = codon_seq[pos//3]
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[pos//3] = random.choice(dt.codon_synonyms[codon1])
                codon2 = codon_seq[(pos//3)+1]
                if len(dt.codon_synonyms[codon2]) > 0:
                    codon_seq[(pos//3)+1] = random.choice(dt.codon_synonyms[codon2])
            elif pos%3 == 1: #nnNNNnnn
                codon0 = codon_seq[(pos-1)//3]
                if len(dt.codon_synonyms[codon0]) > 0:
                    codon_seq[(pos-1)//3] = random.choice(dt.codon_synonyms[codon0])
                codon1 = codon_seq[(pos+2)//3]
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos+2)//3] = random.choice(dt.codon_synonyms[codon1])
                codon2 = codon_seq[((pos+2)//3)+1]
                if len(dt.codon_synonyms[codon2]) > 0:
                    codon_seq[((pos+2)//3)+1] = random.choice(dt.codon_synonyms[codon2])
            elif pos%3 == 2: #nNNNnnnN
                codon0 = codon_seq[(pos-2)//3]
                if len(dt.codon_synonyms[codon0]) > 0:
                    codon_seq[(pos-2)//3] = random.choice(dt.codon_synonyms[codon0])
                codon1 = codon_seq[(pos+1)//3]
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos+1)//3] = random.choice(dt.codon_synonyms[codon1])
                codon2 = codon_seq[((pos+1)//3)+1]
                if len(dt.codon_synonyms[codon2]) > 0:
                    codon_seq[((pos+1)//3)+1] = random.choice(dt.codon_synonyms[codon2])

            kmer_lib[kmer] -=1
    return ''.join(codon_seq)

def codon_optimize(DNA_seq, codon_bias):
    codon_freqs = dt.read_codon_freqs(codon_bias)
    pro_seq = dt.translate(DNA_seq)
    if pro_seq.endswith("*"):
        syn_DNA_seq = dt.reverse_translate(pro_seq.strip("*"), CB = codon_freqs) + "TAA"
    else:
        syn_DNA_seq = dt.reverse_translate(pro_seq, CB = codon_freqs)
    return(syn_DNA_seq)

def fix_GC(DNA_seq, lim_high = 0.6, lim_low = 0.45):
    #Calculate GC content
    GC_init = (DNA_seq.count("G") + DNA_seq.count("C")) / len(DNA_seq)
    
    #If its within range, return the original DNA sequence back
    if GC_init <= lim_high and GC_init >= lim_low:
        return DNA_seq
    
    else:
        #split sequence into codons
        codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
    #If GC is too high
    if GC_init > lim_high:
        attempt_count = 0
        while True:
            attempt_count +=1
            #Abandon ship if it appears impossible to get GC low enough
            if attempt_count > 1000:
                sys.stderr.write("Warning: failed to fix GC to within specified limits, "
                                 "returning best effort\n")
                return(''.join(codon_seq))
            #Go to a random codon in the sequence
            pos = random.randrange(0, len(codon_seq))
            codon_init = codon_seq[pos]
            #Look for a random synonymous codon with lower GC content
            synonyms = random.sample(dt.codon_synonyms[codon_init], len(dt.codon_synonyms[codon_init]))
            for codon_new in synonyms:
                if (codon_new.count("G") + codon_new.count("C")) < (codon_init.count("G") + codon_init.count("C")):
                    #And replace the orignal codon with the new one
                    codon_seq[pos] = codon_new
            #Exit if you've gotten the GC low enough
            if (''.join(codon_seq).count("G") + ''.join(codon_seq).count("C")) / len(DNA_seq) < lim_high:
                return(''.join(codon_seq))

    elif GC_init < lim_low:
        attempt_count = 0
        while True:
            attempt_count +=1
            #Abandon ship if it appears impossible to get GC high enough
            if attempt_count > 1000:
                sys.stderr.write("Warning: failed to fix GC to within specified limits, "
                                 "returning best effort\n")
                return(''.join(codon_seq))
            #Go to a random codon in the sequence
            pos = random.randrange(0, len(codon_seq))
            codon_init = codon_seq[pos]
            #Look for a random synonymous codon with higher GC content
            synonyms = random.sample(dt.codon_synonyms[codon_init], len(dt.codon_synonyms[codon_init]))
            for codon_new in synonyms:
                if (codon_new.count("G") + codon_new.count("C")) > (codon_init.count("G") + codon_init.count("C")):
                    #And replace the orignal codon with the new one
                    codon_seq[pos] = codon_new
            #Exit if you've gotten the GC high enough
            if (''.join(codon_seq).count("G") + ''.join(codon_seq).count("C")) / len(DNA_seq) > lim_low:
                return(''.join(codon_seq))
    else:
        sys.stderr.write("Warning: something went wrong while trying to fix GC content, "
                         "returning sequence unchanged\n")
        return DNA_seq

if __name__ == "__main__":

    from optparse import OptionParser
    import DNA_tools as dt
    import random
    import sys

    parser = OptionParser()
    parser.add_option('--input',
          '-i',
          action = 'store',
          type = 'string',
          dest = 'input_fasta',
          help = "input fasta with your DNA sequence to be optimized")
    parser.add_option('--output',
          '-o',
          action = 'store',
          type = 'string',
          dest = 'output_fasta',
          help = "output fasta with your optimized DNA sequence")
    parser.add_option('--bias',
          '-b',
          action = 'store',
          type = 'string',
          dest = 'bias_file',
          help = "codon bias file",
          default = "../ColiProteomeContent.tsv")

    (option, args) = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.bias_file)
