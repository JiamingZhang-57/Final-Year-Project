clc;
clear all;
close all;
tic;
c = 3e8;
f=20e9;
distance = 0.5;
lambda = c/f;%Maximum of f for multiple frequency
D_left = 0.3;
D_right = 0.295;
D_up = 0.295;
D_down = 0.30;
crossrange_resolution=((c/f)*distance)/(D_left+D_right);
aperture_xtr = -D_left:lambda/2:D_right;
aperture_ytr = -D_down:lambda/2:D_up;
[Xtr,Ytr] = ndgrid(aperture_xtr,aperture_ytr);
za = 0;
%-----------------Working Fundermental Part
%-----------------Try to define the target zone
scene_x = 0.3;
scene_y = 0.3;
pixel_x = -scene_x/2:crossrange_resolution/5:scene_x/2;%Cross-Range resolution for X and Y
pixel_y = -scene_y/2:crossrange_resolution/5:scene_y/2;
pixel_z = distance;
[X_pix,Y_pix,Z_pix] = ndgrid(pixel_x,pixel_y,pixel_z);%Target's coodinator
target=zeros(size(X_pix));
%Gun Target Upright
%target(round(length(pixel_x)/2)-10:round(length(pixel_x)/2)+10,round(length(pixel_y)/2)-5:round(length(pixel_y)/2)+5)=1;
%target(round(length(pixel_x)/2)-4:round(length(pixel_x)/2)+4,round(length(pixel_y)/2)-5:round(length(pixel_y)/2)-4)=1;
target = imread(['knife18.jpeg']);%remember to change save order!!!!!!
target = im2gray(target);
target= im2double(target);
[x,y] = size(X_pix);
target = imresize(target, [x,y]);
target=target./max(target(:));

%% 
k=2*pi*f/c;
s=zeros(size(Xtr));
for i=1:numel(Xtr)
         for j=1:numel(X_pix)
             s(i)=s(i)+target(j)*exp(-2*1i*k*sqrt((Xtr(i)-X_pix(j)).^2+(Ytr(i)-Y_pix(j)).^2+(0-Z_pix(j)).^2));
        end
end   

scene=zeros(size(X_pix));
for j=1:numel(X_pix)
         for i=1:numel(Xtr)
             scene(j)=scene(j)+s(i)*exp(+2*1i*k*sqrt((Xtr(i)-X_pix(j)).^2+(Ytr(i)-Y_pix(j)).^2+(0-Z_pix(j)).^2));
        end
end   






figure();
imagesc(target);
colorbar;
title('target')
hold
figure();
imagesc(abs(scene));
colorbar;
title('Reconstruction for  layers ')
hold
figure()
imagesc(abs(scene));


