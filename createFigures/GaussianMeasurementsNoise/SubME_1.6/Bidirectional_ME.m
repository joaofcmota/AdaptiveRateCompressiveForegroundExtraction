function [MVx, MVy] = Bidirectional_ME(img0, img1, opts)
% Forward ME
[MVx1 MVy1] = Motion_Est(img0, img1, opts);

% Backward ME
[MVx2 MVy2] = Motion_Est(img1, img0, opts);

% Motion Refinement
MVx(:,:,1) =  MVx1;
MVx(:,:,2) = -MVx2;

MVy(:,:,1) =  MVy1;
MVy(:,:,2) = -MVy2;

MVx = max(MVx, [], 3);
MVy = max(MVy, [], 3);