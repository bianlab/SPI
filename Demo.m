%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 22, 2017. Contact me: lihengbian@gmail.com.
% This demo does the simulation of single-pixel imaging, with different reconstruction methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
clc;
close all;
addpath(genpath(pwd));

%% Parameters
num_pixel = 64; % number of image pixels in each dimension
samplingRatio = 0.5;

para.tol = 1e-2;
para.min_iter = 30;
para.x0flag = 0; % initialization flag of the reconstructed image 0: all one; 1: pinv(A)*b.

%% Measurements
im = im2double(imread('cameraman.tif'));
im = imresize(im,[num_pixel,num_pixel]);
figure,imshow(im,[],'InitialMagnification',1000); title('Scene image');

num_pattern = round(samplingRatio * num_pixel * num_pixel); % number of illumination patterns
patterns =  rand(num_pixel,num_pixel,num_pattern);
measurements = sum(sum(repmat(im,[1,1,num_pattern]) .* patterns));
measurements = reshape(measurements,[],1);

%% initialization
[row, col, m] = size(patterns);
P = reshape(patterns, [row*col, m]);
P = P'; % each row represents a pattern

if para.x0flag == 1
    para.x0 = pinv(P)*measurements;
else
    para.x0 = ones(row * col,1);
end

%% 1. SPI differential ghost imaging (DGI) reconstruction
ind = 1;
fprintf('Begin reconstruction of DGI. \n');
tic
[im_r_DGI] = fun_SPI_R_DGI(patterns, measurements);
runTime(ind,2) = toc;
figure,imshow(im_r_DGI,[],'InitialMagnification',1000); title('Recovered im using DGI method');
im_r(:,:,ind) = im_r_DGI;

%% 2. SPI gradient descent (GD) reconstruction
ind = 2;
fprintf('Begin reconstruction of GD. \n');
tic
[im_r_GD, totaliter] = fun_SPI_R_GD(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_GD,[],'InitialMagnification',1000); title('Recovered im using GD method');
im_r(:,:,ind) = im_r_GD;

%% 3. SPI conjugate gradient descent (CGD) reconstruction
ind = 3;
fprintf('Begin reconstruction of CGD. \n');
tic
[im_r_CGD, totaliter] = fun_SPI_R_CGD(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_CGD,[],'InitialMagnification',1000); title('Recovered im using CGD method');
im_r(:,:,ind) = im_r_CGD;

%% 4. SPI Poisson maximum likelihood reconstruction
ind = 4;
fprintf('Begin reconstruction of Poisson. \n');
tic
[im_r_Poisson, totaliter] = fun_SPI_R_Poisson(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_Poisson,[],'InitialMagnification',1000); title('Recovered im using Poisson method');
im_r(:,:,ind) = im_r_Poisson;

%% 5. SPI alternating projection (AP) reconstruction
ind = 5;
fprintf('Begin reconstruction of AP. \n');
tic
[im_r_AP, totaliter] = fun_SPI_R_AP(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_AP,[],'InitialMagnification',1000); title('Recovered im using AP method');
im_r(:,:,ind) = im_r_AP;

%% 6. SPI Sparse representation (DCT) compressive sensing reconstruction
ind = 6;
fprintf('Begin reconstruction of sparse representation (DCT) compressive sensing. \n');
tic
[im_r_Sparse, totaliter] = fun_SPI_R_Sparse(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_Sparse,[],'InitialMagnification',1000); title('Recovered im using sparse representation (DCT) method');
im_r(:,:,ind) = im_r_Sparse;

%% 7. SPI Total variation (TV) compressive sensing reconstruction
ind = 7;
fprintf('Begin reconstruction of total variation (TV) compressive sensing. \n');
tic
[im_r_TV, totaliter] = fun_SPI_R_TV(patterns, measurements, para);
runTime(ind,2) = toc;
runTime(ind,1) = totaliter;
figure,imshow(im_r_TV,[],'InitialMagnification',1000); title('Recovered im using total variation (TV) method');
im_r(:,:,ind) = im_r_TV;

%% show results
for ind = 1:7
    im_r(:,:,ind) = im_r(:,:,ind) - min(min(im_r(:,:,ind)));
    im_r(:,:,ind) = im_r(:,:,ind)/max(max(im_r(:,:,ind)));
    imError(ind) = fun_error(im, im_r(:,:,ind));
end
imError = imError';

figure;
subplot(2,4,1), imshow(im,[],'InitialMagnification',1000); title('Groundtruth');
subplot(2,4,2), imshow(im_r_DGI,[],'InitialMagnification',1000); title('DGI');
subplot(2,4,3), imshow(im_r_GD,[],'InitialMagnification',1000); title('GD');
subplot(2,4,4), imshow(im_r_CGD,[],'InitialMagnification',1000); title('CGD');
subplot(2,4,5), imshow(im_r_Poisson,[],'InitialMagnification',1000); title('Poisson');
subplot(2,4,6), imshow(im_r_AP,[],'InitialMagnification',1000); title('AP');
subplot(2,4,7), imshow(im_r_Sparse,[],'InitialMagnification',1000); title('Sparse');
subplot(2,4,8), imshow(im_r_TV,[],'InitialMagnification',1000); title('TV');
