from Bio.SeqIO.FastaIO import SimpleFastaParser
import random

aa2dna = { 'F': 'TTT', 'L': 'CTT', 'I': 'ATT', 'M': 'ATG', 'V': 'GTT',
           'S': 'TCT', 'P': 'CCT', 'T': 'ACT', 'A': 'GCT', 'Y': 'TAT',
           'H': 'CAT', 'Q': 'CAA', 'N': 'AAT', 'K': 'AAA', 'D': 'GAT',
           'E': 'GAA', 'C': 'TGT', 'W': 'TGG', 'R': 'CGT', 'G': 'GGT'}

aa = [*aa2dna.keys()]

kmer_len = 8

# generate 1000 different k-mers
aa_kmers = []
dna_fwd_kmers = []
dna_rev_kmers = []

rev_cmp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

for i in range(100000):
    aa_seq = [random.choice(aa) for j in range(kmer_len)]
    dna_seq = [aa2dna[c] for c in aa_seq]

    aa_str = ''.join(aa_seq)
    dna_str = ''.join(dna_seq)
    dna_rev_str = ''.join([rev_cmp[c] for c in dna_str[::-1]])

    aa_kmers.append(aa_str)
    dna_fwd_kmers.append(dna_str)
    dna_rev_kmers.append(dna_rev_str)

dna_kmers = [dna_fwd_kmers, dna_rev_kmers]

# generate proteoms and genomes
n_samples = 100

for n_kmers_per_sample in [1, 1000]:
    aa_file = open(f'aa_{n_samples}x{n_kmers_per_sample}.fasta', 'w')
    dna_file = open(f'dna_{n_samples}x{n_kmers_per_sample}.fasta', 'w')

    for i in range(0, n_samples):

        indices = [random.choice(range(0,len(aa_kmers))) for j in range(n_kmers_per_sample)]
        fwd_or_rev = [random.choice([0, 1]) for j in range(n_kmers_per_sample)]

        aa_seq = [aa_kmers[j] for j in indices]
        dna_seq =[dna_kmers[fwd_or_rev[j]][indices[j]] for j in range(n_kmers_per_sample)]

        aa_str = '.'.join(aa_seq)
        dna_str = '.'.join(dna_seq)
        id = f'>seq_{i+1}\n'

        aa_file.write(id)
        aa_file.write(aa_str)
        aa_file.write('\n')

        dna_file.write(id)
        dna_file.write(dna_str)
        dna_file.write('\n')

    aa_file.close()
    dna_file.close()
