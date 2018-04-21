%CS_IFT.m
%
%DESCRIPTION:
%   uses inverse Fourier transform to "reconstruct" incomplete fourier
%   measurements (ie sets undefined values to 0)
%
%INPUTS:
%    y: the vector of incomplete Fourier measurements (M-vector)
%
%    N: the length of the output
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%    flip: whether or not to flip the permutation vector (0 or 1)
%
%OUTPUTS:
%    x: the "reconstructed" signal (N-vector)
%
%NOTES:
%    multiply by 'N' in the last step to account for the way Matlab
%    implements FFT/IFFT (no scaling on forward DFT matrix, 1/N scaling on
%    inverse DFT matrix)
%
%    the 'flip' functionality is useful when using cross validation since
%    the matrix for one call should have completely different rows from a
%    matrix from a second call
%
%REFERENCES:

function x = CS_IFT(y,N,rng_seed,flip)

%input handling
if nargin<4
    flip=0;
end

%get the length of the input
M = length(y);

%reset the state of the stream for repeatability
%rng(rng_seed,'twister');
seed = RandStream('mcg16807','Seed',rng_seed);
RandStream.setDefaultStream(seed);

%generate the indices of the M rows to 'sample' from the full Fourier
%matrix by using the provided stream and setting the substream to a
%constant (known) value
idx = randperm(N);
if(flip)
    idx = fliplr(idx);
end

%create the FT of x, and fill in with the measurement values
X = zeros(N,1);
X( idx(1:M) ) = y;

%take the IDFT of x, and scale output components by values in 'norms' to
%account for unit-normed columns in the row-sampled Fourier matrix
x = real( N/sqrt(M)*ifft(X)); %factor of 'N': see note above
