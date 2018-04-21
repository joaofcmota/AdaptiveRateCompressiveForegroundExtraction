%CS_Gaussian.m
%
%DESCRIPTION:
%   takes compressive measurements of a signal x using a measurement matrix
%   formed by drawing entries from iid Gaussian distributions
%
%INPUTS:
%    x: the input signal (N-vector)
%
%    M: the number of rows to be used in this CS matrix
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%OUTPUTS:
%    y: the compressive measurement (M-vector)
%
%NOTES:
%
%REFERENCES:

function y = CS_Gaussian(x,M,rng_seed)

%length of the input
x = x(:);
N = length(x);

%determine the block size for Gaussian sub-matrix generation
max_block_bytes = 5e6;
block_rows = floor(max_block_bytes/8/N);
full_blocks = floor(M/block_rows);

%loop through blocks
yslice = zeros(block_rows,full_blocks);
parfor i=1:full_blocks    
    %set up the state of the random number generator for this block
    rng(i*rng_seed,'twister');
    
    %calculate the output for this block
    yslice(:,i) = sqrt(1/M)*randn(N,block_rows).' * x;
end
y = yslice(:);

%calculate the contribution from the leftover block, if one exists
if(rem(M,block_rows)~=0)
    %set up the state of the random number generator for this block
    rng((full_blocks+1)*rng_seed,'twister');
    
    %calculate the output for this block
    y = [y; sqrt(1/M)*randn(N,rem(M,block_rows)).' * x];
end