# SCGenoL

Note "amp_p_mat[g,h,x]", is the probability that x will be called by the sequencer given the underlying genotype is gh
(ignoring allelic dropout events). Therefore amp_p_mat[g,g,g] should be very close to 1, amp_p_mat[g,h,g] should be 
very close to 0.5 and amp_p_mat[g,h,i] should be just above or exactly zero. These could come from an empirical source
or be calculated by the program where the above values are adjusted by some small prior amplification error.

seq_p_mat implicitly does the 1/3 calculation in monovar for homozygous likelihoods, by definition.

##TODO
multiprocessing
input pipeline script, creating pileup file and RECALIBRATING(?) base scores
Assumes phred33, maybe include option for phred64
handle indels
handle CNVs
handle no reads for a cell-locus, defaults to favouring homozygous reference
For each site run a quick EM algorithm to obtain f_0 priors
For cells with no read depth use site-wise prior as posterior


EM algorithm with neighbour joining tree(by hamming distance?): EM params: (ado rate, fp rate(?), SFS(?), amplification probabilities(incl perhaps cytosine deamination parameter)
for germline data, if homozygous update reference(?), if heterozygous include germline_SNV info
Q: build tree using posterior probabilities? update priors in a phylogeny aware way, i.e. use 'parent' as reference? perhaps do BFS where posterior of parent * p(child is descendant) = prior of child
If using neighbour joining, all cells are considered leaves and all non-leaf nodes have unknown genotypes. These genotypes could be considered latent variables for an EM algorithm, and (once neighbour joined) maximised using the genotype likelihoods of the cells. Potentially problematic as not truly EM: the neighbour joining is not maximising likelihood(?).
