#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:16:01 2021

@author: bryanandrews1
"""

def main(input_fasta, output_fasta, bias_file, banned = None, is_dna = False,
         seed = None):

    #Seed both random number generators so a run can be reproduced. Python's
    #random drives the cleanup steps (fix_GC, strip_*), while numpy drives the
    #weighted codon draw in dt.reverse_translate. A seed of None (the default)
    #leaves them at OS entropy, i.e. a different valid result each run.
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    entries = dt.read_fasta_multi(input_fasta)
    codon_freqs = dt.read_codon_freqs(bias_file)

    optimized_entries = []
    for seq_name, seq in entries:
        if len(seq) == 0:
            sys.stderr.write("Warning: sequence '%s' is empty, skipping\n" % seq_name)
            continue
        #Turn the input into a protein sequence before optimizing. By default the
        #input is assumed to be protein; with --DNA it is coding DNA that we
        #translate first. Validation happens here (outside the try) so that a
        #hard error can exit rather than being swallowed as a per-sequence skip.
        if is_dna:
            validate_dna_input(seq, seq_name)  #exits on non-DNA characters
            pro_seq = dt.translate(seq)
        else:
            warn_if_dna_input(seq, seq_name)   #warns if a protein looks like DNA
            pro_seq = seq.upper()
        try:
            DNA_seq = codon_optimize(pro_seq, codon_freqs)
            DNA_seq = fix_GC(DNA_seq)
            DNA_seq = strip_microhomology(DNA_seq)
            DNA_seq = strip_mononucleotide_tracts(DNA_seq)
            DNA_seq = strip_G_quadruplexes(DNA_seq)
            #Banned subsequence removal runs last so that no earlier step can
            #reintroduce a forbidden site (e.g. a restriction enzyme site)
            DNA_seq = strip_banned_subsequences(DNA_seq, banned, seq_name)
        except Exception as e:
            sys.stderr.write("Warning: failed to optimize sequence '%s' (%s), skipping\n"
                             % (seq_name, e))
            continue
        optimized_entries.append((seq_name, DNA_seq))

    with open(output_fasta,'w+') as fastout:
        fastout.write('\n'.join(">%s\n%s" % (seq_name, DNA_seq)
                                for seq_name, DNA_seq in optimized_entries))

def validate_dna_input(seq, seq_name = ""):
    #DNA input must contain only unambiguous DNA bases (A/C/G/T). Anything else
    #(other amino-acid letters, gaps, ambiguity codes, whitespace) means the
    #input isn't the coding DNA that --DNA promised, and dt.translate() can't
    #handle it, so we exit loudly rather than silently mangling the sequence.
    bad = sorted(set(seq.upper()) - set("ACGT"))
    if bad:
        label = (" '%s'" % seq_name) if seq_name else ""
        sys.stderr.write("Error: --DNA input sequence%s contains non-DNA "
                         "character(s): %s. Exiting.\n"
                         % (label, ", ".join(bad)))
        sys.exit(1)

def warn_if_dna_input(seq, seq_name = ""):
    #A/C/G/T are also valid amino acids, so a protein composed entirely of those
    #letters is ambiguous -- it may actually be DNA passed without the --DNA
    #flag. Warn but proceed, since it could legitimately be a short peptide.
    if set(seq.upper()) <= set("ACGT"):
        label = (" '%s'" % seq_name) if seq_name else ""
        sys.stderr.write("Warning: protein input sequence%s contains only DNA "
                         "characters (A/C/G/T); did you mean to pass --DNA?\n"
                         % label)

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

def strip_banned_subsequences(DNA_seq, banned, seq_name = "", max_attempts = 2000):
    #Remove arbitrary forbidden subsequences (e.g. restriction enzyme sites) by
    #swapping in synonymous codons, so the encoded protein is unchanged. Because
    #enzymes cut either strand, each banned motif's reverse complement is banned
    #too. This is meant to run last in the pipeline so nothing reintroduces a site.
    if not banned:
        return DNA_seq

    #Build the set of target sequences: each motif plus its reverse complement.
    valid_bases = set("ACGTN")
    targets = []
    seen = set()
    for motif in banned:
        for t in (motif.upper(), _safe_revcomp(motif.upper(), valid_bases)):
            if t and t not in seen:
                seen.add(t)
                targets.append(t)

    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
    blocked = set()  #(target, position) pairs we know can't be broken
    attempt_count = 0
    while True:
        seq_now = ''.join(codon_seq)
        #Find the earliest occurrence of any target that we haven't given up on.
        hit = None
        for t in targets:
            start = 0
            while True:
                pos = seq_now.find(t, start)
                if pos == -1:
                    break
                if (t, pos) not in blocked:
                    if hit is None or pos < hit[0]:
                        hit = (pos, t)
                    break
                start = pos + 1
        if hit is None:
            break
        attempt_count += 1
        if attempt_count > max_attempts:
            break
        pos, t = hit
        #Codons spanning the banned site; only those with synonyms are mutable.
        first_codon = pos // 3
        last_codon = (pos + len(t) - 1) // 3
        candidates = [c for c in range(first_codon, last_codon + 1)
                      if len(dt.codon_synonyms[codon_seq[c]]) > 0]
        if not candidates:
            #No covering codon can change without altering the protein -> give up
            #on this particular occurrence and move on to any others.
            blocked.add((t, pos))
            continue
        c = random.choice(candidates)
        codon_seq[c] = random.choice(dt.codon_synonyms[codon_seq[c]])

    final_seq = ''.join(codon_seq)
    remaining = sorted({t for t in targets if t in final_seq})
    if remaining:
        label = (" in '%s'" % seq_name) if seq_name else ""
        sys.stderr.write("Warning: unable to remove banned subsequence(s)%s: %s "
                         "(returning best effort)\n" % (label, ", ".join(remaining)))
    return final_seq

def _safe_revcomp(motif, valid_bases):
    #Only reverse-complement motifs made of recognized bases; dt.reverse_complement
    #would otherwise sys.exit() on an unexpected character.
    if set(motif) <= valid_bases:
        return dt.reverse_complement(motif)
    return None

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

def codon_optimize(pro_seq, codon_freqs):
    #codon_freqs is a codon-bias table already loaded by dt.read_codon_freqs().
    #A trailing stop (*) is handled separately by appending a stop codon, since
    #the bias table only carries the 20 amino acids.
    pro_seq = pro_seq.upper()
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

    import argparse
    import DNA_tools as dt
    import numpy as np
    import random
    import sys

    parser = argparse.ArgumentParser(
        description = "Optimize protein or DNA sequence(s) for synthesis/assembly")
    parser.add_argument('--input', '-i',
          dest = 'input_fasta',
          help = "input fasta; a protein fasta by default, or coding DNA if --DNA is set")
    parser.add_argument('--output', '-o',
          dest = 'output_fasta',
          help = "output fasta with your optimized DNA sequence")
    parser.add_argument('--bias', '-b',
          dest = 'bias_file',
          default = "../codon_sets/ColiProteomeContent.tsv",
          help = "codon bias file")
    parser.add_argument('--DNA',
          dest = 'is_dna',
          action = 'store_true',
          default = False,
          help = "treat the input as coding DNA to be translated first "
                 "(default: input is protein)")
    parser.add_argument('--ban',
          dest = 'banned',
          nargs = '*',
          default = [],
          metavar = 'SUBSEQ',
          help = "one or more subsequences to remove (e.g. restriction sites); "
                 "the reverse complement of each is removed as well")
    parser.add_argument('--seed',
          dest = 'seed',
          type = int,
          default = None,
          help = "integer seed for reproducible output; default is a random seed")

    option = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.bias_file,
         option.banned, option.is_dna, option.seed)
