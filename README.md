# SCGenoL

## Differences from monovar:


Note "amp_p_mat[g,h,x]", is the probability that x will be called by the sequencer given the underlying genotype is gh
(ignoring allelic dropout events). Therefore amp_p_mat[g,g,g] should be very close to 1, amp_p_mat[g,h,g] should be 
very close to 0.5 and amp_p_mat[g,h,i] should be just above or exactly zero. These could come from an empirical source
or be calculated by the program where the above values are adjusted by some small prior amplification error.

seq_p_mat implicitly does the 1/3 calculation in monovar for homozygous likelihoods, by definition.
### Heterozygous genotype maths:
Calculated for each cell at each locus

Likelihood(g=1) = p(d|g=1)
=p(d, g=ra|g=1) + p(d, g=rb|g=1) + p(d, g=rc|g=1)
=\sum_a p(d|g=ra) * p(g=ra|g=1)
Here we do not favour any heterozygous genotype, and all have likelihood (1/3). This may be changed to reflect empirical or dbSNP data
=1/3 \sum_a p(d|g=ra, ado) * p(ado) + p(d|g=ra, no ado)(1-p(ado))
p(d|g=ra, ado) = p(d|g=ra, drop r) * p(drop r) + p(d|g=ra, drop a) * p(drop a)
Here we assume either allele is equally likely to be dropped in an ado event and p(drop r) = p(drop a) = 0.5. This is unlikely to change
p(d|g=ra, ado, drop a) = \prod_i p(d_i|g=rr)
p(d|g=ra, ado, drop r) = \prod_i p(d_i|g=aa)

##TODO
multiprocessing
input pipeline script, creating pileup file and RECALIBRATING(?) base scores
Assumes phred33, maybe include option for phred64
handle indels
EM algorithm with neighbour joining tree(by hamming distance?): EM params: (ado rate, fp rate(?), SFS(?), amplification probabilities(incl perhaps cytosine deamination parameter)
for germline data, if homozygous update reference(?), if heterozygous include germline_SNV info
Q: build tree using posterior probabilities? update priors in a phylogeny aware way, i.e. use 'parent' as reference? perhaps do BFS where posterior of parent * p(child is descendant) = prior of child
handle CNVs
unify sites assumption: infinite or finite? need to remove g=2 handling.
