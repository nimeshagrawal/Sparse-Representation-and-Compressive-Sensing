# Sparse-Representation-and-Compressive-Sensing
MATLAB implementation of Orthogonal Matching Pursuit(OMP) Algorithm for Compressive Sensing (Signal Recovery)

Consider a sparse signal with sparsity level "K" generated with the non-zero support "S" 
determined by a Bernoulli sequence from the binary field {1,0}, such that s(n) = 1 indicates
that the nth value is non-zero, with P{s(n) = 1} = p and P{s(n) = 0} = 1 - p, and corresponding non-zero
amplitudes x(n) given as Gaussian random variables N(0,1) for the non-zero locations.

A recovery code based on OMP is implemented considering M by N sensing matrix P, whose entries are populated
by random values drawn from a Gaussian distribution, say N(0,1). Conclusions for recovery of signal based on
different number of Measurements (M) is being presented in the report.
