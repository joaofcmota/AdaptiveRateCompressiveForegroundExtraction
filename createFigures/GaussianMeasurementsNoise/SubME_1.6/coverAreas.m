function g_covered = coverAreas(g, img1, simpleExtra);

[m n C]    = size(img1);
g_covered = padarray(g, [1, 1], NaN);

for y = 2:(m+1)
    for x = 2:(n+1)
        if isnan(g_covered(y,x))
            pred = 0;
            numOfpred = 0;
            if isnan( g_covered(y,x-1) ) 
            else
                pred = pred + g_covered(y,x-1);
                numOfpred = numOfpred + 1;
            end
            if isnan( g_covered(y-1,x-1) ) 
            else
                pred = pred + g_covered(y-1,x-1);
                numOfpred = numOfpred + 1;
            end
            if isnan( g_covered(y-1,x) ) 
            else
                pred = pred + g_covered(y-1,x);
                numOfpred = numOfpred + 1;
            end
            pred = pred + simpleExtra(y-1,x-1);
            numOfpred = numOfpred + 1;
            g_covered(y,x) = pred/numOfpred;
        end
    end
end

g_covered = g_covered( 2:(m+1), 2:(n+1));

% for y = 1:m
%     for x = 1:n
%         if isnan(g_covered(y,x))
%             g_covered(y,x) = img1(y,x);
%         end
%     end
% end

