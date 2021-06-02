#For network permutation test

permutationnum=100
category=CS16-ac-hs_gain

#Final result for the permutation test is: results/cordev_all_3-wayOrtho_hs_enhancer_${category}_WC2_4_shuffled_qvalues.txt


mkdir shuffle intersect counts results
#shuffle the human gain features
for num in $(eval echo {1..$permutationnum})
do
	perl shuffle_humangain_140225.pl cordev_all_3-wayOrtho_hs_enhancer.bed cortex_annotation_140226_forcount_enhancer.txt shuffle/cordev_all_3-wayOrtho_hs_enhancer_s$num.bed
done

#intersect with human genes
cut -f 1,2,3,4 gencode.v10.wholeGene.exonTranscript_regDom.bed | intersectBed -a stdin -b cordev_all_3-wayOrtho_hs_enhancer.bed -wb | cut -f 4,8 > intersect/cordev_all_3-wayOrtho_hs_enhancer_intersect.bed;


for num in $(eval echo {1..$permutationnum})
do
	cut -f 1,2,3,4 gencode.v10.wholeGene.exonTranscript_regDom.bed | intersectBed -a stdin -b shuffle/cordev_all_3-wayOrtho_hs_enhancer_s$num.bed -wb | cut -f 4,8 > intersect/cordev_all_3-wayOrtho_hs_enhancer_s$num\_intersect.bed;
done

#count the shuffling

perl count_intersect_module_140227.pl intersect/cordev_all_3-wayOrtho_hs_enhancer_intersect.bed cortex_annotation_140226_forcount_enhancer.txt WC2_4_module_allinfo.txt counts/cordev_all_3-wayOrtho_hs_enhancer_intersect_${category}_WC2_4_count.txt ${category} all 1 1

for num in $(eval echo {1..$permutationnum})
	do
		perl count_intersect_module_140227.pl intersect/cordev_all_3-wayOrtho_hs_enhancer_s$num\_intersect.bed cortex_annotation_140226_forcount_enhancer.txt WC2_4_module_allinfo.txt counts/cordev_all_3-wayOrtho_hs_enhancer_s$num\_intersect_${category}_WC2_4_count.txt ${category} all 1 1;
done

#calculate enrichment test p-value and perform p-value correction
perl summarize_shuffle_results_module_131118.pl counts/cordev_all_3-wayOrtho_hs_enhancer_intersect_${category}_WC2_4_count.txt "counts/cordev_all_3-wayOrtho_hs_enhancer_s*_intersect_${category}_WC2_4_count.txt" results/cordev_all_3-wayOrtho_hs_enhancer_${category}_WC2_4_shuffled_results.txt;
Rscript cal_general_binom_permu_140312.R results/cordev_all_3-wayOrtho_hs_enhancer_${category}_WC2_4_shuffled_results.txt results/cordev_all_3-wayOrtho_hs_enhancer_${category}_WC2_4_shuffled_qvalues.txt
