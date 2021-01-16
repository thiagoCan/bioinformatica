f = open('Tabela_Codon.txt')

dic = {}
for i in f:
    i = i.strip()
    codon, am = i.split('  ')
    dic[codon] = am

print(dic)
