%CS_Bernoulli_T.m
%
%DESCRIPTION:
%   finds the matrix-vector product of the transpose of the Gaussian matrix
%   created by 'CS_Bernoulli.m' with a measurement vector
%
%INPUTS:
%    y: the input "observation" vector
%
%    N: dimensionality of the output (number of columns in the Bernoulli
%    matrix)
%
%    rng_seed: seed for Matlab's 'rng' function to ensure reapeatability
%
%OUTPUTS:
%    x: the matrix-vector product of the transpose of the Bernoulli matrix
%    and y
%
%NOTES:
%
%REFERENCES:

function x = CS_Bernoulli_T(y,N,rng_seed)

%get the length of the input
M = length(y);

%reset the state of the stream for repeatability
%rng(rng_seed,'twister');
seed = RandStream('mcg16807','Seed',rng_seed);
RandStream.setDefaultStream(seed);

%determine the block size for Gaussian sub-matrix generation
memory_block_bytes = 1e9;
memory_block_doubles = memory_block_bytes/8;
memory_block_rows = floor(memory_block_doubles/N);

%loop through blocks
x = zeros(N,1);
for i=1:ceil(M/memory_block_rows)
    row_idx_min = (i-1)*memory_block_rows+1;
    row_idx_max = min(M,i*memory_block_rows);
    x = x + (binornd(1,1/2,[N (row_idx_max-row_idx_min+1)])*2/sqrt(M) - 1/sqrt(M)).'* ...
        y(row_idx_min:row_idx_max);
end
