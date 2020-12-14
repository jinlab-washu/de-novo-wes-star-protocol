# Author: Xue

import sys
file = sys.argv[1]
out = sys.argv[2]

file = open(file,'rt')		#'/ycga-gpfs/project/ysm/lifton/xz374/NDAR_BayOutput/Batch1/NDAR_Batch1.BayesianFilter.DenovoM'
out = open(out,'wt')		#'/ycga-gpfs/project/ysm/lifton/xz374/NDAR_BayOutput/Batch1/NDAR_Batch1.BayesianFilter.DenovoM_parsed'


for line in file:
    Denovo_Candidate, Proband_ID, chr, pos, ref, parents, alt, QS, PASS, AC, AF, AN, Func_refGene, Gene_refGene, ExonicFunc_refGene, AAChange_refGene, genomicSuperDups, snp138, esp6500v2si_all, Thousandg2015aug_all, ExAC_ALL, ExAC_AMR, ExAC_FIN, ExAC_NFE, ExAC_OTH, ExAC_AFR, ExAC_EAS, ExAC_SAS, gnomAD_exome_ALL, gnomAD_genome_ALL, SIFT_score, SIFT_converted_rankscore, SIFT_pred, MetaSVM_score, MetaSVM_rankscore, MetaSVM_pred, CADD_raw, CADD_raw_rankscore, CADD_phred, phyloP100way_vertebrate, phyloP100way_vertebrate_rankscore, SiPhy_29way_logOdds, SiPhy_29way_logOdds_rankscore, PROBAND_GT, PROBAND_ref_cov, PROBAND_nonref_cov, PROBAND_DP, PROBAND_GQ, PROBAND_DGQ, PROBAND_PL_hom1, PROBAND_PL_het, PROBAND_PL_hom2, FATHER_GT, FATHER_ref_cov, FATHER_nonref_cov, FATHER_DP, FATHER_GQ, FATHER_DGQ, FATHER_PL_hom1, FATHER_PL_het, FATHER_PL_hom2, MOTHER_GT, MOTHER_ref_cov, MOTHER_nonref_cov, MOTHER_DP, MOTHER_GQ, MOTHER_DGQ, MOTHER_PL_hom1, MOTHER_PL_het, MOTHER_PL_hom2, Bayes_Factor = line.strip().split('\t')
    if Denovo_Candidate == 'Denovo_Candidate':
        out.write(line)
    if Denovo_Candidate != 'Denovo_Candidate':
        if (((Func_refGene == 'exonic') or (Func_refGene == 'splicing') or (('splicing' in Func_refGene) and ('exonic' in Func_refGene))) \
        and float(PROBAND_DP) >= 10 \
        and float(AC) !=0 \
        and float(PROBAND_nonref_cov) >= 5 \
        and float(FATHER_DP) >= 10 \
        and float(MOTHER_DP) >=10 \
        and ((float(PROBAND_nonref_cov) >= 5 and float(PROBAND_nonref_cov) < 10 and (float(PROBAND_nonref_cov)/float(PROBAND_DP)) >= 0.28) or (float(PROBAND_nonref_cov) >= 10 and (float(PROBAND_nonref_cov)/float(PROBAND_DP)) >= 0.2)) \
        and float(FATHER_nonref_cov)/float(FATHER_DP) < 0.035):
            if (chr == 'Y' and (PROBAND_GT != FATHER_GT)):
                out.write(line)
            else:
                if (float(MOTHER_nonref_cov)/float(MOTHER_DP) < 0.035):
                    out.write(line)

out.close()
file.close()
