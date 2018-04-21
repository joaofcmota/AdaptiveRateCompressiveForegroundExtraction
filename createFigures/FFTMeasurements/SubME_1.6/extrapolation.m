clear all
close all
clc

%% Input parameters
%Video sequence
filename = 'videos/tennis_qcif.yuv';
M1=176; % horizontal dimension of frame
N1=144; %original vertical dimention of frame
downsamplingFactor = 2;
numberOfFrames = 10; %number of frames to downsample
framesize=M1*N1 + (M1/2)*(N1/2) + (M1/2)*(N1/2);
% ME Parameters
opts.BlockSize   = 4;
opts.SearchLimit = 16;
%% Read two consecutive frames

for frameIndex=1:1:numberOfFrames
    close all
    fid = fopen(filename,'r');
    fseek(fid,framesize*frameIndex,-1);
    frame1=fread(fid,[M1,N1],'uchar');
    frame1_u=fread(fid,[M1/2,N1/2],'uchar');
    frame1_v=fread(fid,[M1/2,N1/2],'uchar');
    
    frame2=fread(fid,[M1,N1],'uchar');
    frame2_u=fread(fid,[M1/2,N1/2],'uchar');
    frame2_v=fread(fid,[M1/2,N1/2],'uchar');
    
    frame3=fread(fid,[M1,N1],'uchar');
    frame3_u=fread(fid,[M1/2,N1/2],'uchar');
    frame3_v=fread(fid,[M1/2,N1/2],'uchar');
    downsampledFrame1 = downsample((downsample(frame1,downsamplingFactor))',downsamplingFactor);
    downsampledFrame2 = downsample((downsample(frame2,downsamplingFactor))',downsamplingFactor);
    downsampledFrame3 = downsample((downsample(frame3,downsamplingFactor))',downsamplingFactor);
   % Read image
    img0 = im2double(downsampledFrame1);
    img1 = im2double(downsampledFrame2);
    img2 = im2double(downsampledFrame3);
%  load('walkingGuy.mat');
% img0 = imStack(:,:,2);
% img1 = imStack(:,:,3);
% img2 = imStack(:,:,4);
    % Motion estimation
    %[MVx, MVy] = Bidirectional_ME(img0, img1, opts);
    [MVx, MVy] = Motion_Est(img0, img1, opts);

    % Motion Compensation
    imgMC = reconstruct(img1, MVx, MVy, 1);

    % Evaluation
    [M N C] = size(imgMC);
    Res  = imgMC-img2(1:M, 1:N, 1:C);
    MSE  = norm(Res(:), 'fro')^2/numel(imgMC);
    PSNR(frameIndex) = 10*log10(255^2/MSE);
    
    Res1  = img1-img2(1:M, 1:N, 1:C);
    MSE1  = norm(Res1(:), 'fro')^2/numel(img1);
    PSNR1(frameIndex) = 10*log10(255^2/MSE1);
% Show results
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
% imshow(uint8(imgMC)); title('img_M');
% 
% subplot(224); 
% T = sprintf('img_M - img_2, PSNR %3g dB', PSNR);
% imshow(uint8(Res)); title(T);

end