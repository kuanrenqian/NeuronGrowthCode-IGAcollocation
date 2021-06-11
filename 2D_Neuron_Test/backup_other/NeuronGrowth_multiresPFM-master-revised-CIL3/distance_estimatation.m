%% extract cell center and distance from image

% %%CIL3
img = imread('8786_orig.tif');
img = img(:,:,1);

x1 = [465,226];
x2 = [725,772];
dist = sqrt((x1(1,1)-x2(1,1))^2 + (x1(1,2)-x2(1,2))^2);
distx = x1(1,1)-x2(1,1);
disty = x1(1,2)-x2(1,2);

radius = 28;
theta = 0:0.1:2*pi;
xx = x1(1,1) + radius*cos(theta);
yy = x1(1,2) + radius*sin(theta);

figure
imagesc(img)
colormap gray
hold on
plot(xx,yy);

% %% CIL2
% img = imread('8783_orig.tif');
% img = img(:,:,1);
% 
% x1 = [538,503];
% x2 = [826,383];
% dist = sqrt((x1(1,1)-x2(1,1))^2 + (x1(1,2)-x2(1,2))^2);
% distx = x1(1,1)-x2(1,1);
% disty = x1(1,2)-x2(1,2);
% 
% radius = 28;
% theta = 0:0.1:2*pi;
% xx = x1(1,1) + radius*cos(theta);
% yy = x1(1,2) + radius*sin(theta);
% 
% figure
% imagesc(img)
% colormap gray
% hold on
% plot(xx,yy);

%% CIL4
% img = imread('8787_orig.tif');
% img = img(:,:,1);
% 
% x1 = [463,492];
% x2 = [1011,354];
% dist = sqrt((x1(1,1)-x2(1,1))^2 + (x1(1,2)-x2(1,2))^2);
% distx = x1(1,1)-x2(1,1);
% disty = x1(1,2)-x2(1,2);
% 
% radius = 28;
% theta = 0:0.1:2*pi;
% xx = x1(1,1) + radius*cos(theta);
% yy = x1(1,2) + radius*sin(theta);
% 
% figure
% imagesc(img)
% colormap gray
% hold on
% plot(xx,yy);

% %% CIl 1
% img = imread('8773_orig.tif');
% img = img(:,:,1);
% 
% x1 = [601,489];
% x2 = [823,400];
% dist = sqrt((x1(1,1)-x2(1,1))^2 + (x1(1,2)-x2(1,2))^2);
% distx = x1(1,1)-x2(1,1);
% disty = x1(1,2)-x2(1,2);
% 
% radius = 28;
% theta = 0:0.1:2*pi;
% xx = x2(1,1) + radius*cos(theta);
% yy = x2(1,2) + radius*sin(theta);
% 
% phi = zeros(size(img));
% for i =1:size(phi,1)
%     for j =1:size(phi,2)
%         if((i-x1(1,1))^2+ (j-x1(1,2))^2 <= radius^2)
%             phi(j,i) = 1;
%         end
%         
%         if((i-x2(1,1))^2+ (j-x2(1,2))^2 <= radius^2)
%             phi(j,i) = 1;
%         end
%     end
% end
% 
% figure
% imagesc(img)
% colormap gray
% hold on
% plot(xx,yy);