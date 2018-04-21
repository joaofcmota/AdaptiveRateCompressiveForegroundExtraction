% =========================================================================
% Hall sequence

addpath('../../../../datasets/');
addpath('../../results/VolkansAlg/');
addpath('../../results_img/VolkansAlg/');

load hall128x128.mat;
load ResultsHall128x128_allFrames.mat;

THRESH = 0.1;  % To remove small noise

% Original images
frames = [1 4 100 150 250];
len_fr = length(frames);

for i = 1 : len_fr
   
    fr = frames(i);
    
    figure(1);clf;
    imshow(uint8(imStack(:,:,fr)));
    drawnow;
    saveas(gcf, ['Hall_original_f', num2str(fr), '.pdf']);        
end

% estimated image
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['Hall_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(uint8(img_pred));
    drawnow;
    saveas(gcf, ['Hall_estimated_image_f', num2str(fr), '.pdf']);        
end

% reconstructed image
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['Hall_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(uint8(f_k));
    drawnow;
    saveas(gcf, ['Hall_reconstructed_image_f', num2str(fr), '.pdf']);        
end

% reconstructed foreground
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['Hall_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(abs(f_k-imStack(:,:,1))>THRESH);
    drawnow;
    saveas(gcf, ['Hall_reconstructed_foreground_f', num2str(fr), '.pdf']);        
end
% =========================================================================

%%
% =========================================================================
% PETS sequence
addpath('../../../../datasets/');
addpath('../../results/VolkansAlg/');
addpath('../../results_img/VolkansAlg/');

load PETS09_S2L1.mat;
load ResultsPETS_allFrames.mat;

THRESH = 0.1;  % To remove small noise

% Original images
frames = [1,5,75,100,170];
len_fr = length(frames);

for i = 1 : len_fr
   
    fr = frames(i);
    
    figure(1);clf;
    imshow(uint8(imStack(:,:,fr)));
    drawnow;
    saveas(gcf, ['PETS_original_f', num2str(fr), '.pdf']);        
end

% estimated image
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['PETS_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(uint8(img_pred));
    drawnow;
    saveas(gcf, ['PETS_estimated_image_f', num2str(fr), '.pdf']);        
end

% reconstructed image
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['PETS_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(uint8(f_k));
    drawnow;
    saveas(gcf, ['PETS_reconstructed_image_f', num2str(fr), '.pdf']);        
end

% reconstructed foreground
for i = 2 : len_fr     % not the background image
   
    fr = frames(i);
    
    load(['PETS_visualsframe_' , num2str(fr)])
        
    figure(2);clf;
    imshow(abs(f_k-imStack(:,:,1))>THRESH);
    drawnow;
    saveas(gcf, ['PETS_reconstructed_foreground_f', num2str(fr), '.pdf']);        
end


% =========================================================================
























