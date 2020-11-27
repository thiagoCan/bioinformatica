from bioferramentas import *

a = RNA_DNA('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGUUGA')
print(a)

b = fitaCompl(a)
print('fitaCompl')
print(b)

c = trans('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGUUGA')
print(c)

d = validarDNA('TACCGGTACCGCGGGTCTTGACTCTAGTTATCATGGGCATAATTGCCAACT')
print(d)

e = gc_content(c)
print('gc_content')
print(e)

