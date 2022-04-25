clc;
clear all;
close all;
tic;
c = 3e8;
f=20e9;
Num_masks=1000;
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


H=zeros(Num_masks,numel(X_pix));%Measurement Matrix
E_txr=zeros(Num_masks,numel(X_pix));
E_trx=zeros(Num_masks,numel(X_pix));
k=2*pi*f/c;
X_pix_s=ones(numel(aperture_xtr),numel(aperture_xtr));% Make a matrix which has same size as Xtr
Y_pix_s=ones(numel(aperture_ytr),numel(aperture_ytr));% Make a matrix which has same size as Ytr
zero_matrix=zeros(size(Y_pix_s));%Make a zero matrix same size as Y_pix_s above
Z_pix_s=ones(size(zero_matrix));

for N=1:Num_masks
        real_m=hadamard(numel(aperture_xtr));
       imag_m=hadamard(numel(aperture_xtr));
       % M= normrnd(0,1,[numel(aperture_xtr),numel(aperture_ytr)])+(1j*normrnd(0,1,[numel(aperture_xtr),numel(aperture_ytr)]));
       a=randperm(numel(aperture_xtr),1);
       b=randperm(numel(aperture_xtr),1);
       c=randperm(numel(aperture_xtr),1);
       real_m=circshift(real_m,[a,b]);
       imag_m=circshift(imag_m,[b,c]);
       M=real_m+1j*imag_m;
       for j=1:numel(X_pix)   
                X_p_s=X_pix(j)*X_pix_s;
                Y_p_s=Y_pix(j)*X_pix_s;
                Z_p_s=Z_pix(j)*Z_pix_s;
                E_txr(N,j)=sum(M.*exp(-1*1i*k*sqrt((Xtr-X_p_s).^2+(Ytr-Y_p_s).^2+(0-Z_p_s).^2)),'all');%Using matrix array calculation
                E_trx(N,j)=sum(M.*exp(-1*1i*k*sqrt((Xtr-X_p_s).^2+(Ytr-Y_p_s).^2+(0-Z_p_s).^2)),'all');
                %for i=1:numel(Xtr)
                %   E_txr(N,j)=E_txr(N,j)+M(i)*exp(-1*1i*k*sqrt((Xtr(i)-X_pix(j)).^2+(Ytr(i)-Y_pix(j)).^2+(0-Z_pix(j)).^2));
                %end
                %for i=1:numel(Xtr)
                %    E_trx(N,j)=E_trx(N,j)+M(i)*exp(-1*1i*k*sqrt((Xtr(i)-X_pix(j)).^2+(Ytr(i)-Y_pix(j)).^2+(0-Z_pix(j)).^2));
                %end
           H(N,j)=E_trx(N,j).*E_txr(N,j);
        end
end    
g=zeros(Num_masks,1);
for t=1:Num_masks
    for j=1:numel(X_pix)
        g(t)=g(t)+(H(t,j).*target(j));
    end
    
end
g=awgn(g,20,'measured');
%% 

%-------------------------------Reconstruction Part(Need to modify)
scene_est=zeros(numel(target),1);
scene_est=H'*g;
f_est=gmres(H'*H,scene_est,1,1e-04,15);
f_est=f_est./max(abs(f_est(:)));

scene=zeros(size(X_pix));
for i=1:numel(f_est)%reshape f_est into scene matrix
    scene(i)=f_est(i);
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

%save('.mat','scene','-v7.3');                                        %save mat!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!
%load('Trainset.mat');
%save('.mat','Trainset','-v7.3'); 
%save('Trainoise.mat','Trainoise','-v7.3'); 
%save('Trainoisepred.mat','Trainoisepred','-v7.3'); 
%save('Trainnoise.mat','Trainnoise','-v7.3')
%Trainset=zeros(120,120,1,50);
%Trainset(:,:,1,1)=abs(scene);
%save('Trainset.mat','Trainset','-v7.3');
%Trainnoise(:,:,1,x)=Trainoise(:,:,1,y)[1,0,0,1,0,1,0,1,1,0];


s = svd(H);
s=s./max(s);
figure(); semilogy(s,'LineWidth',2); % We use semilogy for plotting singular values. This plots the vertical axis on a logarithmic grid (the singular values are still linear - just plotted on a log grid). 
xlabel('Number of measurement modes');
ylabel('Normalized Singular values');
set(gca,'FontSize',14)
set(gcf,'color','w');
grid on
toc
%%
%load('Trainset.mat')
%A(:,:,i)=scenei
%load('.mat')
%train=zeros(120,120,1,)
%train(:,:,1,)=abs(scene)
%Trainnoise(:,:,1,1)=abs(scene)