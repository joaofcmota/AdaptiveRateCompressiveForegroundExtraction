clear all
close all
clc

%% Input parameters
%Video sequence
%filename = 'videos/tennis_qcif.yuv';
%filename = 'videos/ice_qcif.yuv';
filename = 'videos/hall_qcif.yuv';
M1=176; % horizontal dimension of frame
N1=144; %original vertical dimention of frame
downsamplingFactor = 2;
numberOfFrames = 100; %number of frames to downsample
framesize=M1*N1 + (M1/2)*(N1/2) + (M1/2)*(N1/2);
% ME Parameters
opts.BlockSize   = 8;
opts.SearchLimit = 5;
%% Read two consecutive frames
load('videos/PETS09_S2L1.mat');
imStack1 = imStack(:,1:116,:);
numberOfFrames = 100;
M1 = siz(1); N1 = siz(1);
for frameIndex=3:numberOfFrames
    close all
%     fid = fopen(filename,'r');
%     fseek(fid,framesize*frameIndex,-1);
%     frame1=fread(fid,[M1,N1],'uchar');
%     frame1_u=fread(fid,[M1/2,N1/2],'uchar');
%     frame1_v=fread(fid,[M1/2,N1/2],'uchar');
%     
%     frame2=fread(fid,[M1,N1],'uchar');
%     frame2_u=fread(fid,[M1/2,N1/2],'uchar');
%     frame2_v=fread(fid,[M1/2,N1/2],'uchar');
%     
%     frame3=fread(fid,[M1,N1],'uchar');
%     frame3_u=fread(fid,[M1/2,N1/2],'uchar');
%     frame3_v=fread(fid,[M1/2,N1/2],'uchar');
    
% %     % Read image
%     img0 = im2double(frame1');
%     img1 = im2double(frame2');
%     img2 = im2double(frame3');

img0 = im2double(imStack1(:,:,frameIndex-2));
img1 = im2double(imStack1(:,:,frameIndex-1));
img2 = im2double(imStack1(:,:,frameIndex));

%     img0 = [ 0 0 255 255 255 255 0 0 0 0 0 0 0 0 0 0;  
%              0 0 255 255 255 255 0 0 0 0 0 0 0 0 0 0;
%              0 0 255 255 255 255 0 0 0 0 0 0 0 0 0 0;
%              0 0 255 255 255 255 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%          
%     img1 = [ 0 0 0 0 255 255 255 255 0 0 0 0 0 0 0 0;  
%              0 0 0 0 255 255 255 255 0 0 0 0 0 0 0 0;
%              0 0 0 0 255 255 255 255 0 0 0 0 0 0 0 0;
%              0 0 0 0 255 255 255 255 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%              
%     img2 = [ 0 0 0 0 0 0 255 255 255 255 0 0 0 0 0 0;  
%              0 0 0 0 0 0 255 255 255 255 0 0 0 0 0 0;
%              0 0 0 0 0 0 255 255 255 255 0 0 0 0 0 0;
%              0 0 0 0 0 0 255 255 255 255 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
%     A = 20+6.*randn(8,16);
%          
%     img0 = im2double(uint8(img0+A));
%     img1 = im2double(uint8(img1+A));
%     img2 = im2double(uint8(img2+A));
    % Motion estimation
    %[MVx, MVy] = Bidirectional_ME(img0, img1, opts);
    [MVx, MVy] = Motion_Estnew(img0, img1, opts);
    %[MVxOld, MVyOld] = Motion_Est(img0, img1, opts);
    
    % simple extrapolation via Motion Compensation
    %imgMC = reconstruct(img1, MVxOld, MVyOld, 1);
    %MVx = [0 2 0 0; 0 0 0 0]; MVy = [0 0 0 0; 0 0 0 0];
    % extrapolation using a linear motion model
    g_mce = mceNew(img1, img1, MVx, MVy, opts.SearchLimit);

    % Evaluation
    [M N C] = size(g_mce);
%     Res  = imgMC-img2(1:M, 1:N, 1:C);
%     MSE  = norm(Res(:), 'fro')^2/numel(imgMC);
%     PSNR(frameIndex) = 10*log10(255^2/MSE);
    
    Res1  = g_mce-img2(1:M, 1:N, 1:C);
    MSE1  = norm(Res1(:), 'fro')^2/numel(g_mce);
    PSNR1(frameIndex) = 10*log10(255^2/MSE1);

    %Show results
% figure(1);
% quiver(MVx(end:-1:1,:), MVy(end:-1:1,:));
% title('Motion Vector Field');
% 
% figure(2);
% subplot(221);
% imshow(uint8(img1)); title('img_1');
% 
% subplot(222);
% imshow(uint8(img2)); title('img_2');
% 
% subplot(223);
% imshow(uint8(g_mce)); title('g_mce');
% 
% subplot(224); 
% %T = sprintf('img_M - img_2, PSNR %3g dB', PSNR);
% imshow(uint8(Res1)); title('Res1');
% pause(3);
end

disp('Average PSNR: '); disp(mean(PSNR1(3:end)));

figure(3)
plot(1:numberOfFrames, PSNR1)%, 1:numberOfFrames, PSNR);