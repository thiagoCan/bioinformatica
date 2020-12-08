ami = {
    'UUU' : 'F',      
    'CUU' : 'L',      
    'AUU' : 'I',      
    'GUU' : 'V',
    'UUC' : 'F',      
    'CUC' : 'L',      
    'AUC' : 'I',      
    'GUC' : 'V',
    'UUA' : 'L',      
    'CUA' : 'L',      
    'AUA' : 'I',      
    'GUA' : 'V',
    'UUG' : 'L',      
    'CUG' : 'L',      
    'AUG' : 'M',      
    'GUG' : 'V',
    'UCU' : 'S',      
    'CCU' : 'P',      
    'ACU' : 'T',      
    'GCU' : 'A',
    'UCC' : 'S',      
    'CCC' : 'P',      
    'ACC' : 'T',      
    'GCC' : 'A',
    'UCA' : 'S',      
    'CCA' : 'P',      
    'ACA' : 'T',      
    'GCA' : 'A',
    'UCG' : 'S',      
    'CCG' : 'P',      
    'ACG' : 'T',      
    'GCG' : 'A',
    'UAU' : 'Y',      
    'CAU' : 'H',      
    'AAU' : 'N',      
    'GAU' : 'D',
    'UAC' : 'Y',      
    'CAC' : 'H',      
    'AAC' : 'N',      
    'GAC' : 'D',
    'UAA' : 'Stop',   
    'CAA' : 'Q',      
    'AAA' : 'K',      
    'GAA' : 'E',
    'UAG' : 'Stop',   
    'CAG' : 'Q',      
    'AAG' : 'K',      
    'GAG' : 'E',
    'UGU' : 'C',      
    'CGU' : 'R',      
    'AGU' : 'S',      
    'GGU' : 'G',
    'UGC' : 'C',      
    'CGC' : 'R',      
    'AGC' : 'S',      
    'GGC' : 'G',
    'UGA' : 'Stop',   
    'CGA' : 'R',      
    'AGA' : 'R',      
    'GGA' : 'G',
    'UGG' : 'W',      
    'CGG' : 'R',      
    'AGG' : 'R',      
    'GGG' : 'G'
}

nucleotídios = ['A', 'T', 'C', 'G']

dnaCompl = {'A' : 'T',
             'T' : 'A',
             'C' : 'G',
             'G' : 'C'
             }

def trans(rna):
    """traduzir o RNAm para aminoácidos"""
    rna = rna.upper()
    protein = ''
    
    for i in range(0, len(rna), 3):
        if ami[rna[i:i+3]] != 'Stop':
            protein += ami[rna[i:i+3]]

    return protein

def RNA_DNA(a):
    """RNA para DNA"""
    a = a.upper()

    return a.replace('U', 'T')

def fitaCompl(a):
    """fita de DNA complementar"""
##    return ''.join([dnaCompl[i]   for i in a])

    mapping = a.maketrans('ATCG', 'TAGC')
    return a.translate(mapping)

def validarDNA(a):
    for i in a:
        if i not in nucleotídios:
            return False
    return True

def gc_content(b):
    return round((b.count('C') + b.count('G'))/len(b) * 100)

def gc_content_subseq(a, k):
    r = []
    for i in range(0, len(a), k):
        subseq = a[i:i+k]
        r.append(gc_content(subseq))

    return r
#==========================================================
FASTAdic = {}

f = open('gc_content.txt')

for i in f:
    i = i.replace('\n', '')
    if '>' in i:
        índice = i
        FASTAdic[índice] = ''
    else:
        FASTAdic[índice] += i

f.close()
print(FASTAdic)


bob = {key: gc_content(value)    for (key, value) in FASTAdic.items()}
print(bob)

MaxGCKey = max(bob, key=bob.get)
print ('MaxGCKey ', bob[MaxGCKey])

#====================================================================================

def aa_different(aa1, aa2):
    count = 0
    for i in range(len(aa1)):
        if aa1[i] != aa2[i]:
            count += 1

    return count
#===============================================================================

def pos_inicial(p, seq):

    comp_p = len(p)
    comp_seq = len(seq)

    if comp_p > comp_seq: return -1

    pulo = [comp_p  for k in range(256)]
    for u in range(comp_p - 1): pulo[ord(p[u])] = comp_p - u - 1
    pulo = tuple(pulo)

    s = comp_p - 1
    while s < comp_seq:
        j = comp_p - 1
        i = s

        while j >= 0 and seq[i] == p[j]:
            j -= 1
            i -= 1

        if j == -1: return i + 1
        s += pulo[ord(seq[s])]

    return -1


#===============================================================================
























    
