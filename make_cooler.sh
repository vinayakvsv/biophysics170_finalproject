# cooler csort
cooler csort --nproc 4 -c1 1 -p1 2 -s1 3 -c2 4 -p2 5 -s2 6 -o test_pairwisemat_1000x1000.20kb.sorted test_pairwisemat_1000x1000.txt hg-data/hg19.ucsc.20kb.binsizes.bed
# cooler cload
cooler cload pairix hg-data/hg19.ucsc.20kb.binsizes.bed test_pairwisemat_1000x1000.20kb.sorted test_pairwisemat_1000x1000.20kb.sorted.cool
# cooler balance
cooler balance -p 6 test_pairwisemat_1000x1000.20kb.sorted.cool --force
#cooler balance -p 6 test_pairwisemat_1000x1000.10kb.sorted.cool --force
