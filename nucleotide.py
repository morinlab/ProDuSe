import itertools
import random
import math


IUPAC = {
    'A':['A'],
    'C':['C'],
    'G':['G'],
    'T':['T'],
    'W':['A','T'],
    'S':['C','G'],
    'M':['A','C'],
    'K':['G','T'],
    'R':['A','G'],
    'Y':['C','T'],
    'B':['C','G','T'],
    'D':['A','G','T'],
    'H':['A','C','T'],
    'V':['A','C','G'],
    'N':['A','C','G','T']
    }

CAPUI = {
    'A':'A', 
    'C':'C', 
    'G':'G', 
    'T':'T', 
    'AC':'M', 
    'AG':'R',
    'AT':'W',
    'CG':'S',
    'CT':'Y',
    'GT':'K', 
    'ACG':'V',
    'ACT':'H',
    'AGT':'D',
    'CGT':'B', 
    'ACGT':'N'
    }

DIST = {
    'A':[1,0,0,0],
    'C':[0,1,0,0],
    'G':[0,0,1,0],
    'T':[0,0,0,1],
    'W':[0.5,0,0,0.5],
    'S':[0,0.5,0.5,0],
    'M':[0.5,0.5,0,0],
    'K':[0,0,0.5,0.5],
    'R':[0.5,0,0.5,0],
    'Y':[0,0.5,0,0.5],
    # 'B':[0.0,1.0/3.0,1.0/3.0,1.0/3.0],
    # 'D':[1.0/3.0,0.0,1.0/3.0,1.0/3.0],
    # 'H':[1.0/3.0,1.0/3.0,0.0,1.0/3.0],
    # 'V':[1.0/3.0,1.0/3.0,1.0/3.0,0.0],
    'N':[0.25,0.25,0.25,0.25]
    }


COMPLEMENT = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    'N':'N'
    }

CAPITAL = {
    'a':'A',
    'A':'A',
    'c':'C',
    'C':'C',
    'g':'G',
    'G':'G',
    't':'T',
    'T':'T',
    'n':'N',
    'N':'N'
    }

VARIANT = {
    'A':['C','G','T'],
    'C':['A','G','T'],
    'G':['A','C','T'],
    'T':['A','C','G']
    }

BASE_TO_INDEX = {
    'A':0,
    'C':1,
    'G':2,
    'T':3,
    'N':4
    }

### Check if characters in seq are IUPAC nucleotide symbols
def is_iupac( seq ):
    is_seq_iupac = True
    for base in seq:
        is_seq_iupac = is_seq_iupac and base in IUPAC
    return is_seq_iupac


### Determines if seq is ambiguous (all symbols map to a single nucleotide)
def is_unambiguous( seq ):
    is_seq_unambiguous = True
    for base in seq:
        is_seq_unambiguous = is_seq_unambiguous and len(IUPAC[base]) == 1
    return is_seq_unambiguous


### Determines if seq is ambiguous (at least one symbol maps to two nucleotides)
def is_ambiguous( seq ):
    is_seq_ambiguous = False
    for base in seq:
        is_seq_ambiguous = is_seq_ambiguous or not len(IUPAC[base]) == 1
    return is_seq_ambiguous


### Determines /home/malbuque/projects/produse/produse2/ProDuSeical AND of seq being a subset of ref across all pos
def is_match( seq, ref ):
    is_seq_match = True
    for i in range(0, len(seq)):
        is_seq_match = is_seq_match and set(IUPAC[seq[i]]).issubset(IUPAC[ref[i]])
    return is_seq_match


### Produces all possible strings in order to make seq into ref. 
def make_reference( seq, ref ):
    # If seq[i] is a subset of ref[i], then seq is either more specific
    # or exactly the same as ref, so at position i an unambiguous
    # nucleotide should be from seq[i] 
    # Otherwise seq[i] is not a subset of ref[i], meaning that there
    # exists some mismatches between the unambiguous characters in seq
    # when compared to those in ref, so at position i an unambiguous 
    # nucleotide should be from ref[i]
    nucleotides_lists = []
    for i in range(0, len(seq)):
        if set(IUPAC[seq[i]]).issubset(IUPAC[ref[i]]):
            nucleotides_lists.append(IUPAC[seq[i]])
        else:
            nucleotides_lists.append(IUPAC[ref[i]])
    
    # Converts nucleotide lists into all possible char combinations
    char_lists = itertools.product(*nucleotides_lists)
    
    # Joins the character combinations into their appropriate strigs
    unambiguous_nucleotides =[]
    for chars in char_lists:
        unambiguous_nucleotides.append(''.join(chars))
    
    return unambiguous_nucleotides


### Produces all possible unambiguous nucleotide sequences from seq
def make_unambiguous( seq ):
    # Store all the unambiguous nucleotides in a list for each
    # particular base
    nucleotide_lists = []
    for base in seq:
        nucleotide_lists.append(IUPAC[base])
    # Convert this list into all possible combinations of those
    # characters
    char_lists = itertools.product(*nucleotide_lists)

    # Join the character lists into their appropriate strings
    unambiguous_nucleotides = []
    for chars in char_lists:
        unambiguous_nucleotides.append(''.join(chars))

    return unambiguous_nucleotides

def make_ambiguous(list_of_seq):

    ambiguous_list = []

    for i in range(len(list_of_seq[0])):
        ambiguous_list.append([])

    for seq in list_of_seq:
        for i in range(len(seq)):
            if not seq[i] in ambiguous_list[i]:
                ambiguous_list[i].append(seq[i])

    ambiguous_seq = ''

    for seq in ambiguous_list:
        seq.sort()
        ambiguous_seq = ''.join([ambiguous_seq, CAPUI[''.join(seq)]])

    return ambiguous_seq

def which_ambiguous( seq ):
    indexes = []
    for i in range(len(seq)):
        if len(IUPAC[seq[i]]) > 1:
            indexes.append(i)

    return indexes

def distance( seq, ref, indexes = [] ):
    distance = 0
    if indexes == []:
        indexes = range(len(seq))
    for i in indexes:
        if not seq[i] == 'N' and not ref[i] == 'N':
            if not set(IUPAC[seq[i]]).issubset(IUPAC[ref[i]]):
                distance = distance + 1
    return distance

def base_mismatch( seq, ref, previous = [] ):
    output = previous;
    if previous == []:
        output = [0] * len(seq)
    for i in range(len(seq)):
        if not seq[i] == 'N' and not ref[i] == 'N':
            if not set(IUPAC[seq[i]]).issubset(IUPAC[ref[i]]):
                output[i] += 1
    return output

def dist(first, second):
    tmp_dist = float(0)
    for i in range(len(first)):
        tmp_dist += (float(first[i]) - float(second[i])) ** 2
    return math.sqrt(tmp_dist)

class Counter:

    def __init__(self, list_of_bases, ref):
        self.counts = [0] * 5
        for base in list_of_bases:
            self.counts[BASE_TO_INDEX[base]] += 1
        self.alt = [4]
        self.ref = BASE_TO_INDEX[ref]
        for i in range(4):
            if i == self.ref:
                continue
            elif self.countsself.counts[self.alt[0]] and self.counts0:
                self.alt = [i]
            elif self.counts[i] == self.counts[self.alt[0]] and self.counts0:
                self.alt.append(i)
        self.alt = random.sample(self.alt, 1)[0]
        if self.alt == 4:
            self.alt = self.ref

    def get_count(self,base):
        return self.counts[BASE_TO_INDEX[base]]

    def get_alt(self):
        return INDEX_TO_BASE[self.alt]

    def get_ref(self):
        return INDEX_TO_BASE[self.ref]

def random_unambiguous(seq, n=1):
    all_output = [None] * n
    for j in range(n):
        output = list(seq)
        for i in range(len(seq)):
            output[i] = random.choice(IUPAC[seq[i]])
        all_output[j] = ''.join(output)
    return all_output

def random_mismatch(seq, prob):
    output = list(seq)
    for i in range(len(seq)):
        if random.uniform(0,1) <= prob:
            output[i] = random.choice(IUPAC["N"])
    return ''.join(output)

def dist_list(list_of_seq, indexes = []):
    n = len(list_of_seq)
    distances = [0] * (math.factorial(n) / math.factorial(2) / math.factorial(n-2))
    n_index = 0
    for i in range(n):
        for j in range(i+1, n):
            distances[n_index] = distance(list_of_seq[i], list_of_seq[j], indexes)
            n_index += 1
    return distances

def dist_pairwise(indexes=[], *lists_of_seqs):
    n = len(lists_of_seqs)
    distances = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(len(lists_of_seqs[i])):
                for l in range(len(lists_of_seqs[j])):
                    distances.append(distance(lists_of_seqs[i][k], lists_of_seqs[j][l], indexes))
    return distances

def complement( seq ):
    complementSeq = list(seq)
    for i in range(len(seq)):
        complementSeq[i] = COMPLEMENT[CAPITAL[seq[i]]]
    return ''.join(complementSeq)

def reverse( seq ):
    reverseSeq = list(seq)
    length = len(seq)
    for i in range(length):
        reverseSeq[i] = seq[length - i - 1]
    return ''.join( reverseSeq )

def reverseComplement( seq ):
    reverseComplementSeq = list(seq)
    length = len(seq)
    for i in range(length):
        reverseComplementSeq[i] = COMPLEMENT[CAPITAL[seq[length - i - 1]]]
    return ''.join(reverseComplementSeq)

def make_variant( ref ):
    return random.sample(VARIANT[ref], 1)[0]

def makeCapital ( base ):
    return CAPITAL[base];
