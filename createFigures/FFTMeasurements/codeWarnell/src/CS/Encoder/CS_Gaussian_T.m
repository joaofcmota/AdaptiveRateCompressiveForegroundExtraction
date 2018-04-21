%CS_Gaussian_T.m
%
%DESCRIPTION:
%   finds the matrix-vector product of the transpose of the Gaussian matrix
%   created by 'CS_Gaussian.m' with a measurement vector
%
%INPUTS:
%    y: the measurement vector (M-vector)
%
%    N: the length of the output
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%OUTPUTS:
%    x: the "reconstructed" signal (N-vector)
%
%NOTES:
%
%REFERENCES:

function x = CS_Gaussian_T(y,N,rng_seed)

%get the length of the input
y = y(:);
M = length(y);

%determine the block size for Gaussian sub-matrix generation
max_block_bytes = 5e6;
block_rows = floor(max_block_bytes/8/N);
full_blocks = floor(M/block_rows);

%loop through blocks
yslice = reshape(y(1:block_rows*full_blocks),block_rows,full_blocks);
x = zeros(N,1);
parfor i=1:full_blocks
    %set up the state of the random number generator for this block
    rng(i*rng_seed,'twister');
    
    %add the contribution from this block
    x = x + sqrt(1/M)*randn(N,block_rows) * yslice(:,i);
end

%add the contribution from the leftover block, if one exists
if(rem(M,block_rows)~=0)
    %set up the state of the random number generator for this block
    rng((full_blocks+1)*rng_seed,'twister');
    
    %add the contribution from this block
    x = x + sqrt(1/M)*randn(N,rem(M,block_rows)) * y((block_rows*full_blocks+1):end);
end