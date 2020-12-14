import re
import sys

Familyno = sys.argv[1]

file = open(('Trio_' + Familyno + '_reGT_step2_modified.vcf'),'rt')
removed = open(('Trio_' + Familyno + '_removed.txt'),'wt')
output = open(('Trio_' + Familyno + '_updated.vcf'),'wt')

for line in file:
    if re.match('^[^#]',line) is None:
        output.write(line)
    else:
        AC = (line.strip().split('\t')[7]).strip().split(';')[0]
        PL = (line.strip().split('\t')[8]).strip().split(':')[-1]
        Pro_GT = ((line.strip().split('\t')[9]).strip().split(';')[0]).strip().split(':')[0]
        Pro_PL = (((line.strip().split('\t')[9]).strip().split(';')[0]).strip().split(':')[-1]).strip().split(',')
        Mo_GT = ((line.strip().split('\t')[10]).strip().split(';')[0]).strip().split(':')[0]
        Mo_PL = (((line.strip().split('\t')[10]).strip().split(';')[0]).strip().split(':')[-1]).strip().split(',')
        Fa_GT = ((line.strip().split('\t')[11]).strip().split(';')[0]).strip().split(':')[0]
        Fa_PL = (((line.strip().split('\t')[11]).strip().split(';')[0]).strip().split(':')[-1]).strip().split(',')
        if AC == 'AC=0':
            removed.write(line)
            continue
        if PL != 'PL':
            removed.write(line)
            continue
        if ((Pro_GT == './.') or (Fa_GT == './.') or (Mo_GT == './.')):
            removed.write(line)
            continue
        if ((len(Pro_PL) != 3) or (len(Fa_PL) != 3) or (len(Mo_PL) != 3)):
            removed.write(line)
            continue
        else:
            output.write(line)


file.close()
removed.close()
output.close()

