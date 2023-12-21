def reverse_complement(seq):
    ntComplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    revSeqList = list(reversed(seq))
    revComSeqList = [ntComplement[k] for k in revSeqList]
    revComSeq = ''.join(revComSeqList)
    return revComSeq
seq = ''
with open('C:/Users/lcc/Desktop/snp/primer/tol2-PEGFP-N-3_D10.seq') as f:
    for line in f:
        line = line.rstrip()
        seq += line.upper()
        print(reverse_complement(seq))


def DNA_complement(sequence):sequence = sequence.upper()sequence = sequence.replace('A', 't')sequence = sequence.replace('T', 'a')sequence = sequence.replace('C', 'g')sequence = sequence.replace('G', 'c')return sequence.upper() def DNA_reverse(sequence):sequence = sequence.upper()return sequence[::-1]












