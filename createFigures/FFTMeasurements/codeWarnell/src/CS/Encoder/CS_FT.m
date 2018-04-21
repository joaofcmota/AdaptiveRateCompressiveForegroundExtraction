%CS_FT.m
%
%DESCRIPTION:
%   takes compressive measurements of a signal x using a measurement matrix
%   formed from selecting rows of the DFT matrix
%
%INPUTS:
%    x: the input signal (N-vector)
%
%    M: the number of rows to be used in this CS matrix
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%    flip: whether or not to flip the permutation vector (0 or 1)
%
%OUTPUTS:
%    y: the compressive measurement (M-vector)
%
%NOTES:
%    the 'flip' functionality is useful when using cross validation since
%    the matrix for one call should have completely different rows from a
%    matrix from a second call
%
%REFERENCES:

function y = CS_FT(x,M,rng_seed,flip)

%input handling
if nargin<4
    flip=0;
end

%length of the input
N = length(x);

%reset the state of the stream for repeatability
%rng(rng_seed,'twister');
seed = RandStream('mcg16807','Seed',rng_seed);
RandStream.setDefaultStream(seed);

%generate the M rows to 'sample' from the full Fourier matrix by using the
%provided stream and setting the substream to a constant (known) value
idx = randperm(N);
if( flip )
    idx = fliplr(idx);
end

%elementwise multiply x by values in 'norms' to ensure the sampled Fourier
%matrix has unit-normed columns
x = 1/sqrt(M)*x;

%take the DFT of x
X = fft(x);

%keep only the M entries corresponding to the first M indices in idx
y = X(idx(1:M));

%ensure column vector output
y = y(:);
