import re

mass_table = {
    "A": 89.047678,
    "C": 121.019749,
    "D": 133.037508,
    "E": 147.053158,
    "F": 165.078979,
    "G": 75.032028,
    "H": 155.069477,
    "I": 131.094629,
    "K": 146.105528,
    "L": 131.094629,
    "M": 149.051049,
    "N": 132.053492,
    "O": 255.158292,
    "P": 115.063329,
    "Q": 146.069142,
    "R": 174.111676,
    "S": 105.042593,
    "T": 119.058243,
    "U": 168.964203,
    "V": 117.078979,
    "W": 204.089878,
    "Y": 181.073893,
}

codons = {
    "UUU": "F",
    "CUU": "L",
    "AUU": "I",
    "GUU": "V",
    "UUC": "F",
    "CUC": "L",
    "AUC": "I",
    "GUC": "V",
    "UUA": "L",
    "CUA": "L",
    "AUA": "I",
    "GUA": "V",
    "UUG": "L",
    "CUG": "L",
    "AUG": "M",
    "GUG": "V",
    "UCU": "S",
    "CCU": "P",
    "ACU": "T",
    "GCU": "A",
    "UCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "UCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "UCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "UAU": "Y",
    "CAU": "H",
    "AAU": "N",
    "GAU": "D",
    "UAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "UAA": "Stop",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "UAG": "Stop",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "UGU": "C",
    "CGU": "R",
    "AGU": "S",
    "GGU": "G",
    "UGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "UGA": "Stop",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "UGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G"
}


def parse(path):
    """
    Function parses a fasta file and returns a dictionary with seq_ids and sequences
    :param path: path to fasta file
    :return: dictionary with seq_id as keys and sequence for seq_id as values
    """
    with open(path, 'r') as fasta:
        result = dict()
        for line in fasta:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
                result[seq_id] = ""
            else:
                result[seq_id] = f'{result[seq_id]}{line.strip()}'
    return result


def translate(rna_seq):
    """
    Function translates input RNA sequence to protein sequence
    :param rna_seq: string which represents target RNA sequence (it can start from AUG (M) or not,
    but it must have length multiplied by 3)
    :return: protein sequence for input RNA sequence
    """
    rna_seq = re.sub(" +", '', rna_seq.strip())
    assert not len(rna_seq) % 3, "Length of RNA sequences must multiple by 3"
    assert not 'T' in rna_seq, "RNA sequences can't contain Tymin"
    protein_seq = ""
    for codon in range(0, len(rna_seq), 3):
        protein_seq = f'{protein_seq}{codons[rna_seq[codon: codon + 3]]}'
    return protein_seq


def calc_mass(rna_seq):
    """
    Function receives RNA sequence as input and returns mass for encoded protein
    :param rna_seq: string which represents target RNA sequence (it can start from AUG (M) or not,
    but it must have length multiplied by 3)
    :return: summary mass for protein encoded by input RNA
    """
    rna_seq = re.sub(" +", '', rna_seq.strip()).upper()
    assert not len(rna_seq) % 3, "Length of RNA sequences must multiple by 3"
    assert not 'T' in rna_seq, "RNA sequences can't contain Tymin"
    protein_seq = translate(rna_seq)
    mass = float(0)
    for aa in protein_seq:
        mass += mass_table[aa]
    return mass


def orf(rna_seq):
    """
    Function receives RNA sequence as input and returns all possible ORFs
    If for one start codon there are several stops, the first one will be chosen
    :param rna_seq: string which represents target RNA sequence (it can start from AUG (M) or not,
    but it must have length multiplied by 3)
    :return: all possible ORFs within input RNA sequence (if for one start codon there are several stops,
    the first one will be chosen)
    """
    assert not len(rna_seq) % 3, "Length of RNA sequences must multiple by 3"
    assert not 'T' in rna_seq, "RNA sequences can't contain Tymin"
    rna_seq = re.sub(" +", '', rna_seq.strip())
    starts = list()
    stops = list()
    seqs = list()

    for codon in range(0, len(rna_seq) - 2):
        if codons[rna_seq[codon: codon + 3]] == "M":
            starts.append(codon)
        elif codons[rna_seq[codon: codon + 3]] == "Stop":
            stops.append(codon)
    starts.sort()
    stops.sort()
    for start in starts:
        for stop in stops:
            if start < stop and not ((stop - start) % 3) and (stop - start > 3):
                seq = ""
                for codon in range(start, stop, 3):
                    if codons[rna_seq[codon: codon + 3]] == "Stop":
                        break
                    seq += codons[rna_seq[codon: codon + 3]]
                seqs.append(seq)
                break
    return seqs
