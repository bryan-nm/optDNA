# optDNA

A small toolkit of command-line Python scripts for turning protein or DNA
sequences into **synthesis- and assembly-friendly DNA**. Given a protein (or an
already-coding DNA sequence), optDNA picks codons according to an organism's
codon-usage bias and then rewrites the sequence to remove features that cause
trouble during gene synthesis and DNA assembly (Gibson, Golden Gate, etc.) —
without ever changing the encoded protein.

Everything is plain Python (or `numpy`), and the scripts are designed to be chained together.

---

## What problem does this solve?

The same protein can be encoded by an astronomical number of DNA sequences
(most amino acids have 2–6 synonymous codons). Some of those DNA sequences are
much easier to synthesize and assemble than others. optDNA searches through
synonymous encodings to produce a sequence that:

- **uses codons preferred by your target host** (better expression),
- **has balanced GC content** (roughly 45–60%),
- **avoids repeated 8-mers** (microhomology that misroutes homology-based
  assembly),
- **avoids long single-base runs** (homopolymers that trip up synthesis and
  sequencing),
- **avoids G-quadruplex-forming G-runs** (`GGG`), and
- **avoids arbitrary forbidden motifs** you specify, such as restriction-enzyme
  recognition sites — on *both* strands.

Because every rewrite swaps one codon for a synonymous codon, the translated
protein is guaranteed to be identical to the input.

---

## Repository layout

```
optDNA/
├── src/
│   ├── DNA_tools.py                  # shared library: FASTA I/O, translation, codon tables
│   ├── table_to_fasta.py             # table (.csv/.tsv) -> FASTA
│   ├── reverse_translate_proteins.py # protein FASTA -> DNA FASTA (codon choice only)
│   ├── optimize_DNA_for_assembly.py  # the core optimizer (codon-optimize + clean up)
│   ├── add_flanking_sequences.py     # prepend/append cloning adapters
│   └── fasta_to_table.py             # FASTA -> table (.csv/.tsv)
└── codon_sets/
    ├── ColiProteomeContent.tsv       # E. coli codon bias (proteome-weighted)
    ├── ColiGenomeContent.tsv         # E. coli codon bias (genome-weighted)
    ├── HumanProteomeContent.tsv
    ├── HumanGenomeContent.tsv
    ├── YeastProteomeContent.tsv
    ├── YeastGenomeContent.tsv
    └── Neutral.tsv                   # equal weight to every synonymous codon
```

---

## Requirements

- Python 3 (developed and tested on 3.13; 3.8+ should be fine)
- `numpy`

```bash
pip install numpy
```

All scripts `import DNA_tools`, so run them **from inside the `src/`
directory** (or add `src/` to your `PYTHONPATH`). The default codon-bias paths
in the scripts are written relative to `src/` (`../codon_sets/...`), which is
another reason to run from there.

```bash
cd src
```

---

## The typical workflow

A common end-to-end run looks like this:

```
protein table ──▶ table_to_fasta ──▶ optimize_DNA_for_assembly ──▶ add_flanking_sequences ──▶ fasta_to_table
   (.csv)          (protein .fasta)      (clean DNA .fasta)            (order-ready .fasta)         (.csv)
```

`optimize_DNA_for_assembly.py` takes a **protein FASTA by default** and produces
clean, optimized DNA in one step — so most runs don't need a separate
reverse-translation step. You don't need every step, either:

- Already have a coding **DNA** FASTA? Feed it to the optimizer with the `--DNA`
  flag and it will translate first.
- Just want to reverse-translate a protein **without** the assembly cleanup? Use
  `reverse_translate_proteins.py` instead.

The table/FASTA converters are just conveniences for getting data in and out of
spreadsheets.

---

## The scripts, step by step

### 1. `table_to_fasta.py` — get sequences out of a spreadsheet

Converts a delimited table into a FASTA. By default it reads column 1 as the
name, column 2 as the sequence, and treats the first row as a header.

```bash
python table_to_fasta.py -i proteins.csv -o proteins.fasta
```

What it does:
- Parses the table with Python's `csv` module using the chosen delimiter.
- Skips the header row (unless `--no-header` is given), blank lines, and rows
  that are missing the requested columns or have an empty name/sequence
  (a warning is printed for each skip).
- Writes one `>name` / `sequence` record per row.

Key options:
- `-d, --delimiter` — column delimiter. Accepts escapes like `'\t'` and the
  aliases `tab` / `comma` (default `,`).
- `-n, --name-col` — 1-based column number for names (default `1`).
- `-s, --seq-col` — 1-based column number for sequences (default `2`).
- `--no-header` — set this if the table has no header row.

### 2. `reverse_translate_proteins.py` — protein → DNA (codon choice only)

Turns a protein FASTA into a DNA FASTA by choosing a codon for each amino acid.
This step performs **no** cleanup — it only picks codons. It's a standalone
convenience for when you want codon choice *without* the assembly optimization;
the optimizer in step 3 does its own reverse-translation, so you don't need to
run this first to feed it.

```bash
python reverse_translate_proteins.py -i proteins.fasta -o dna.fasta -b ../codon_sets/HumanProteomeContent.tsv
```

What it does:
- Loads the codon-bias table (see [Codon-bias files](#codon-bias-files)).
- For each protein, draws a codon per residue **randomly, weighted by the
  normalized frequency** in the bias table (`numpy.random.choice`). So a codon
  used 70% of the time in the host is chosen ~70% of the time here. This means
  output is **stochastic** — re-running produces a different (equally valid) DNA
  sequence.
- A trailing `*` (stop) in the protein is handled specially: it's stripped and a
  literal `TAA` stop codon is appended, because the bias tables only carry the
  20 amino acids.
- Empty sequences are skipped with a warning.

Options: `-i/--input`, `-o/--output`, `-b/--bias` (default
`../codon_sets/ColiProteomeContent.tsv`).

### 3. `optimize_DNA_for_assembly.py` — the core optimizer

This is the heart of the toolkit. By default it takes a **protein FASTA** and
returns a cleaned-up DNA FASTA. Pass `--DNA` to give it a **coding DNA FASTA**
instead, in which case it translates the DNA to protein first and then proceeds
identically. For each sequence it runs the following pipeline in order (`main`
in [optimize_DNA_for_assembly.py](src/optimize_DNA_for_assembly.py)):

```bash
# protein input (default)
python optimize_DNA_for_assembly.py -i proteins.fasta -o dna.optimized.fasta \
    -b ../codon_sets/ColiProteomeContent.tsv \
    --ban GAATTC GGATCC   # e.g. EcoRI and BamHI sites

# coding DNA input
python optimize_DNA_for_assembly.py -i dna.fasta -o dna.optimized.fasta --DNA \
    -b ../codon_sets/ColiProteomeContent.tsv
```

The per-sequence steps:

1. **`codon_optimize`** — Reverse-translate the protein into DNA using the
   codon-bias table, choosing codons by the same weighted-random draw as step 2.
   (With `--DNA`, the input DNA is translated to protein first, so the sequence
   is fully re-encoded in the host's preferred codons.) This gives the later
   steps synonymous room to work. A trailing stop is preserved as `TAA`.

2. **`fix_GC`** — Measure GC content. If it's already within **45%–60%**, leave
   it alone. Otherwise, repeatedly pick a random codon and swap it for a
   synonymous codon that moves GC in the needed direction, until the sequence is
   back in range (or after 1000 attempts, at which point it warns and returns
   the best effort).

3. **`strip_microhomology`** — Build a library of every 8-mer in the sequence
   (8 bp is the shortest window guaranteed to span two full codons). For any
   8-mer that appears more than once, mutate the codons covering the *first*
   occurrence to synonymous codons, breaking the repeat. Repeated 8-mers are a
   major cause of misassembly in homology-based methods.

4. **`strip_mononucleotide_tracts`** — Find any run of 5+ identical bases
   (`AAAAA`, `CCCCC`, `GGGGG`, `TTTTT`) and swap the codon(s) overlapping the run
   for synonyms to break it. Repeats until none remain (up to 1000 attempts).

5. **`strip_G_quadruplexes`** — Find any `GGG` and mutate a covering codon to a
   synonym to eliminate it, removing G-runs that can seed G-quadruplex secondary
   structure (up to 100 attempts).

6. **`strip_banned_subsequences`** — Remove arbitrary forbidden motifs supplied
   via `--ban`. For each motif it also bans the **reverse complement** (enzymes
   cut either strand). It finds the earliest remaining banned site and mutates a
   synonymous codon spanning it; if a site can't be broken without changing the
   protein it's marked "blocked" and skipped. **This step runs last on purpose**,
   so that no earlier cleanup step can reintroduce a forbidden site. If a banned
   motif truly cannot be removed (no synonymous codons available at that
   position), it warns and returns the best effort.

Notes on behavior:
- Every mutation is a synonymous codon swap (via the `codon_synonyms` table in
  `DNA_tools.py`), so **the protein never changes**.
- The process is **randomized**; different runs yield different valid sequences.
  There is no fixed seed.
- If a single sequence fails at any step, it's skipped with a warning and the
  rest continue.

Options:
- `-i/--input`, `-o/--output`
- `-b/--bias` — codon-bias file (default `../codon_sets/ColiProteomeContent.tsv`).
- `--DNA` — treat the input as coding DNA to be translated first (default: input
  is protein).
- `--ban SUBSEQ [SUBSEQ ...]` — one or more motifs to remove (e.g. restriction
  sites). The reverse complement of each is removed too. Optional; default is
  none.
- `--seed INT` — integer seed for reproducible output. With a seed, re-running on
  the same input produces the exact same DNA. Omit it (the default) for a fresh
  random result each run.

### 4. `add_flanking_sequences.py` — add cloning adapters

Prepends a head sequence and appends a tail sequence to every entry — handy for
adding Gibson/Golden Gate overhangs, primer-binding sites, or other constant
cloning adapters.

```bash
python add_flanking_sequences.py -i dna.optimized.fasta -o dna.flanked.fasta \
    --head atcagata --tail gagatcag
```

What it does:
- Writes `head + sequence + tail` for each entry, **preserving case exactly** as
  provided — so you can pass the flanks in lowercase to visually mark the
  non-coding adapter regions. Empty sequences are skipped with a warning.

Options: `-i/--input`, `-o/--output`, `--head` (5′ prepend, default empty),
`--tail` (3′ append, default empty).

> The flanks are added verbatim and are **not** re-checked by the optimizer, so
> add them *after* optimization and be aware they can reintroduce banned sites or
> GC/repeat issues at the junctions if you're not careful.

### 5. `fasta_to_table.py` — get sequences back into a spreadsheet

The inverse of `table_to_fasta.py`. Converts a FASTA into a two-column table
(`Name`, `Sequence`), which is convenient for ordering from a synthesis vendor.

```bash
python fasta_to_table.py -i dna.flanked.fasta -o order.csv
```

Options: `-i/--input`, `-o/--output`, `-d/--delimiter` (same delimiter handling
as above; default `,`).

---

## Codon-bias files

The `codon_sets/` directory holds tab-separated codon-usage tables. Each has
four columns:

| column       | meaning                                                         |
|--------------|-----------------------------------------------------------------|
| `codon`      | the three-letter codon (e.g. `GCA`)                             |
| `amino_acid` | the one-letter amino acid it encodes (e.g. `A`)                 |
| `raw_freq`   | overall frequency of the codon (informational)                 |
| `norm_freq`  | frequency **among synonymous codons for that amino acid** (sums to 1 per amino acid) |

`norm_freq` is what the reverse-translation actually samples from. For example,
in `ColiProteomeContent.tsv` the four alanine codons `GCA/GCC/GCG/GCT` have
`norm_freq` values that add up to 1.0, and a codon is chosen with that
probability.

Provided tables:

- **`ColiProteomeContent.tsv` / `ColiGenomeContent.tsv`** — *E. coli*
- **`HumanProteomeContent.tsv` / `HumanGenomeContent.tsv`** — human
- **`YeastProteomeContent.tsv` / `YeastGenomeContent.tsv`** — *S. cerevisiae*
- **`Neutral.tsv`** — every synonymous codon weighted equally (useful when you
  want maximal freedom for the cleanup steps and don't care about host bias)

"Proteome" tables weight codons by how often the protein they help encode is
actually expressed; "Genome" tables weight by raw gene content. To use your own
organism, create a TSV with the same four columns and point `--bias` at it.

---

## A complete example

Starting from a spreadsheet of proteins (`proteins.csv`, with a header row and
columns `name,protein`) and targeting *E. coli*, while removing EcoRI (`GAATTC`)
and BsaI (`GGTCTC`) sites:

```bash
cd src

# 1. spreadsheet -> protein FASTA
python table_to_fasta.py -i ../proteins.csv -o proteins.fasta

# 2. codon-optimize the proteins and clean up for assembly,
#    banning two restriction sites
python optimize_DNA_for_assembly.py -i proteins.fasta -o dna.opt.fasta \
    -b ../codon_sets/ColiProteomeContent.tsv --ban GAATTC GGTCTC

# 3. add cloning adapters (lowercase to mark them)
python add_flanking_sequences.py -i dna.opt.fasta -o dna.flanked.fasta \
    --head aatgata --tail gagatca

# 4. write an order-ready table
python fasta_to_table.py -i dna.flanked.fasta -o order.csv
```

If you were starting from a coding **DNA** FASTA instead of proteins, step 2
would be the same command with `--DNA` added.

---

## Good to know

- **Output is stochastic.** Reverse translation and every cleanup step use
  random synonymous choices. By default there is no fixed seed, so two runs on
  the same input give two different, equally valid sequences — re-run if you want
  a different candidate. Pass `--seed INT` to `optimize_DNA_for_assembly.py` to
  make a run reproducible (it seeds both Python's `random` and numpy's RNG).
- **The protein is always preserved.** Every edit is a synonymous codon swap.
- **Best-effort guarantees.** GC-fixing, tract/quadruplex removal, and banned-
  motif removal all have attempt caps. If a constraint can't be satisfied
  (e.g. a required codon has no synonym), the script prints a warning to
  `stderr` and returns the best sequence it found rather than failing.
- **DNA input must be in-frame coding DNA.** When you pass `--DNA` to
  `optimize_DNA_for_assembly.py`, it translates on frame 0 and assumes a length
  that is a multiple of 3. Protein input has no such requirement.
- **`DNA_tools.py`** is a shared utility library (FASTA readers, translation
  table, reverse complement, codon synonym/frequency tables, plus some extra
  helpers like a PCR simulator, fragment assembler, and Needleman–Wunsch
  aligner that the CLI scripts don't currently use). You generally won't call it
  directly.
