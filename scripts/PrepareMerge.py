import re
import sys

Familyno = sys.argv[1]

file1 = open(('Trio_' + Familyno + '.denovo.Bayfilter.vcf'),'rt')
output1 = open(('Trio_' + Familyno + '.denovo.Bayfilter.content.txt'),'wt')
file2 = open(('Trio_' + Familyno + '.hg19_multianno.vcf'),'rt')
output2 = open(('Trio_' + Familyno + '.hg19_multianno.content.txt'),'wt')

for line1 in file1:
    if re.match('^[^#]',line1) is not None:
        output1.write(line1)

for line2 in file2:
    if re.match('^[^#]',line2) is not None:
        output2.write(line2)
        
file1.close()
file2.close()
output1.close()
output2.close()
