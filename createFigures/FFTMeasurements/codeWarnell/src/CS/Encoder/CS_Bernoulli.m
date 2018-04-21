%CS_Bernoulli.m
%
%DESCRIPTION:
%
%INPUTS:
%    x: the input column vector with which to take the matrix-vector
%    product with the Bernoulli matrix
%
%    M: number of rows in the Bernoulli matrix
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%OUTPUTS:
%    y: the matrix-vector product of the generated Bernoulli matrix and the
%   input x
%
%NOTES:
%
%REFERENCES:

function y = CS_Bernoulli(x,M,rng_seed)

%length of the input
N = length(x);

%reset the state of the stream for repeatability
%rng(rng_seed,'twister');
seed = RandStream('mcg16807','Seed',rng_seed);
RandStream.setDefaultStream(seed);

%determine the block size for Gaussian sub-matrix generation
memory_block_bytes = 1e9;
memory_block_doubles = memory_block_bytes/8;
memory_block_rows = floor(memory_block_doubles/N);

%loop through blocks
y = zeros(M,1);
for i=1:ceil(M/memory_block_rows)
    row_idx_min = (i-1)*memory_block_rows+1;
    row_idx_max = min(M,i*memory_block_rows);
    y(row_idx_min:row_idx_max) = ...
        (binornd(1,1/2,[N (row_idx_max-row_idx_min+1)]).'*2/sqrt(M) - 1/sqrt(M))*x;
end
