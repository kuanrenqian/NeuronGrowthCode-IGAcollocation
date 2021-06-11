%% Clip image data
close all
clear all
clc

img = imread('8787_orig.tif');
img1 = zeros(1120,1120,3);
img2 = img(1:1030,81:1200,:);
img1(26:1055,:,:) = img2;
figure
imagesc(img2)
set(gca,'position',[0 0 1 1],'units','normalized')