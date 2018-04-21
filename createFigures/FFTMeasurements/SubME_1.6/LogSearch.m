function [MVy MVx] = LogSearch(Block, img_ref, xc, yc, SearchLimit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Search algorithm
%
% [MVy MVx] = LogSearch(Block, img_ref, xc, yc, SearchLimit)
% finds the motion vector of the (yc, xc)-th Block in the reference image
% img_ref.
% 
% The search steps are log scaled: coarse step in the beginning and finer
% steps later.
% 
% Input:
% Block         - the current block being searched
% img_ref       - the reference image
% xc, yc        - (xc, yc) is the center coordinate of Block
%
% Output:
% [MVy MVx]     - the motion vector of Block
%
% Stanley Chan
%  3 Jun, 2009
%  5 May, 2010 bug corrected, thanks to Ramsin Khoshabeh.
% 29 Jun, 2010 modified for color images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M N C]       = size(img_ref);
BlockSize   = size(Block,1);
L           = floor(BlockSize/2);
BlockRange  = -L:L-1;
SADmin      = 1e6;
y_min       = yc;
x_min       = xc;

% Rejection
if (yc<SearchLimit+L)||(yc>M-(SearchLimit+L))
    error('Can you set yc >%3g pixels from the boundary? \n',SearchLimit+L);
end

if (xc<SearchLimit+L)||(xc>N-(SearchLimit+L))
    error('Can you set xc >%3g pixels from the boundary? \n',SearchLimit+L);
end

x0 = xc;
y0 = yc;

% Main Loop
LevelMax   = 2;
LevelLimit = zeros(1,LevelMax+1);

for k = 1:LevelMax
    LevelLimit(k+1)             = max(2^(floor(log2(SearchLimit))-k+1),1);
    c                           = 2.^(0:log2(LevelLimit(k+1)));
    c(c+sum(LevelLimit(1:k))>SearchLimit) = [];
    x_range                     = zeros(1,2*length(c)+1);
    x_range(1)                  = 0;
    x_range(2:2:2*length(c))    = c;
    x_range(3:2:2*length(c)+1)  = -c;
    y_range = x_range;
    
    for i = 1:length(y_range)
        for j = 1:length(x_range)
            if SADmin>1e-3
                xt = x0 + x_range(j);
                yt = y0 + y_range(i);
                Block_ref  = img_ref(yt+BlockRange, xt+BlockRange, :);
                SAD        = sum(abs(Block(:) - Block_ref(:)))/(BlockSize^2);
                if SAD < SADmin
                    SADmin  = SAD;
                    x_min   = xt;
                    y_min   = yt;
                end
            else
                SADmin = 0;
                x_min  = xc;
                y_min  = yc;
            end
            
        end
    end
    x0 = x_min;
    y0 = y_min;
end

% Motion Vector (integer part)
MVx_int = xc - x_min;
MVy_int = yc - y_min;

% Taylor Refinement
Block_ref   = img_ref(y_min+BlockRange, x_min+BlockRange, :);
Taylor_sol  = Taylor_App(Block, Block_ref);

% Motion Vector (fractional part)
MVx_frac   = Taylor_sol(1);
MVy_frac   = Taylor_sol(2);

% Motion Vector (overall)
MVx = MVx_int + MVx_frac;
MVy = MVy_int + MVy_frac;
end





% Taylor Refinement
function x = Taylor_App(f, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor Refinement
% 
% This function computes the motion vector using Taylor series
% approximation.
% f(x + dx, y + dy) ~= f(x,y) + dx df/dx + dy df/dy
%
% Stanley Chan
% 3 Jun, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dfx dfy] = gradient(f);

a = sum(dfx(:).^2);
b = sum(dfx(:).*dfy(:));
d = sum(dfy(:).^2);

z = g-f;
p = sum(z(:).*dfx(:));
q = sum(z(:).*dfy(:));

A = [a b; b d];
rhs = [p;q];

if cond(A)>1e6
    x = [0 0]';
else
    x = A\rhs;
end
end