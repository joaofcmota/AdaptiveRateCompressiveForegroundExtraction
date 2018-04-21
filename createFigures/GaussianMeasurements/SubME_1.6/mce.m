function g = mce(img1, MVx, MVy, SearchLimit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integer pel motion compensated extrapolation
%
% g = mce(img1, MVx, MVy)
% constructs a motion compensated extrapolation frame from img1 according to the motion
% vectors specified by MVx and MVy (which are estimated based on img0 and img1)
%
% 
% Nikos Deligiannis
% 29 Oct, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BlockSize  = floor(size(img1,1)/size(MVx,1));
[m n C]    = size(img1);
M          = floor(m/BlockSize)*BlockSize;
N          = floor(n/BlockSize)*BlockSize;
f_prev     = img1(1:M, 1:N, 1:C);
g_mce      = zeros(M, N, C);
g_mce = g_mce + NaN;

% Enlarge the image boundaries by BlockSize/2 pixels
f_prev  = padarray(f_prev,  [BlockSize/2 BlockSize/2], 'replicate');
g_mce = padarray(g_mce, [BlockSize/2 BlockSize/2], 'replicate');

% Pad zeros to images to fit SearchLimit
f_prev  = padarray(f_prev,  [SearchLimit, SearchLimit]);
g_mce = padarray(g_mce, [SearchLimit, SearchLimit]);
map = zeros(size(g_mce,1), size(g_mce,2));

% Set parameters
[M N C]     = size(f_prev);
L           = floor(BlockSize/2);
BlockRange  = -L:L-1;
xc_range    = SearchLimit+L+1 : BlockSize : N-(SearchLimit+L);
yc_range    = SearchLimit+L+1 : BlockSize : M-(SearchLimit+L);


% Main Loop
for i = 1:length(yc_range)
    for j = 1:length(xc_range)
        xc = xc_range(j);
        yc = yc_range(i);
        
        Block = f_prev(yc + BlockRange, xc + BlockRange, :);
        MVx1 = MVx(i,j);
        MVy1 = MVy(i,j);
        % pass the block pixel values to the predicted frame
        % update the predicted frame
        [g_mce, map] = extrapBlockValues(Block, g_mce, map, xc, yc, MVx1, MVy1);
    end
end

g = g_mce( ((BlockSize/2)+SearchLimit+1):(((BlockSize/2)+SearchLimit+1)+m-1), ((BlockSize/2)+SearchLimit+1):(((BlockSize/2)+SearchLimit+1)+n-1));
g_map = map( ((BlockSize/2)+SearchLimit+1):(((BlockSize/2)+SearchLimit+1)+m-1), ((BlockSize/2)+SearchLimit+1):(((BlockSize/2)+SearchLimit+1)+n-1));
g = g./g_map;
for y = 1:m
    for x = 1:n
        if isnan(g(y,x))
            g(y,x) = img1(y,x);
        end
    end
end


