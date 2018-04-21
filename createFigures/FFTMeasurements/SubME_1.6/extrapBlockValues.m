function [g_mce, map] = extrapBlockValues(Block, g_mce, map, xc, yc, MVx1, MVy1)


%% Parameters
BlockSize   = size(Block,1);
L           = floor(BlockSize/2);
% BlockRange: -L:L-1;

x_mce = xc + MVx1;
y_mce = yc + MVy1;

%% Full Search Loop
% for i = -L:L-1
%     for j = -L:L-1
%         
%         x_mce_b = x_mce + j;
%         y_mce_b = y_mce + i;
%         
%         map(y_mce_b,x_mce_b) = map(y_mce_b,x_mce_b) + 1;
%         if isnan(g_mce(y_mce_b,x_mce_b))
%             g_mce(y_mce_b,x_mce_b) = Block(L+i+1,L+j+1);
%         else
%             g_mce(y_mce_b,x_mce_b) = g_mce(y_mce_b,x_mce_b) + Block(L+i+1,L+j+1);
%         end
%     end
% end

for i = 0:(BlockSize-1)
    for j = 0:(BlockSize-1)
        
        x_mce_b = x_mce + j;
        y_mce_b = y_mce + i;
        
        map(y_mce_b,x_mce_b) = map(y_mce_b,x_mce_b) + 1;
        if isnan(g_mce(y_mce_b,x_mce_b))
            g_mce(y_mce_b,x_mce_b) = Block(i+1,j+1);
        else
            g_mce(y_mce_b,x_mce_b) = g_mce(y_mce_b,x_mce_b) + Block(i+1,j+1);
        end
    end
end