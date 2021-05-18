close all;
clear all;
clc;

%%

lambda=633*10^(-9);
d=(3.5+40)*10^(-2); %(4+3)/2=3.5
M=1024; N=1024;
delta_chi=10*10^-6; delta_eta=10*10^-6; 
theta_x = 100*1024*lambda;
theta_y = 0.05;

% Object-square

object=zeros(N,M); 
object(500:524,500:524)=1; 

%Object and its spectrum

FFT_of_the_object=fftshift(fft2(object)); % Following the instructions
figure;
suptitle('The object and his spectrum');
subplot(1,2,1);
imshow(object); colormap(gray); title('Object'); axis on;
subplot(1,2,2);
imshow(5.*mat2gray(abs(FFT_of_the_object))); 
colormap(gray); 
title('The spectrum of the object'); 
axis on;

% Creting CCD meshgrid M  
M1=1:M; %pixels
W_CCD=delta_chi*M;
m_CCD=delta_chi.*M1-W_CCD/2; 
deltaY=(lambda*d)/(M*delta_eta);

% Creting CCD meshgrid  N
M1=1:N; %pixels
W_CCD=delta_chi*M;
n_CCD=delta_eta.*M1-W_CCD/2;
deltaX=(lambda*d)/(N*delta_chi); 

% Calculation of Q_ccd
[X_ccd,Y_ccd]=meshgrid(m_CCD,n_CCD);
Q_CCD=exp((1i*pi)/(lambda*d)*(X_ccd.^2+Y_ccd.^2)); 

figure;
imagesc(angle(Q_CCD));
title('phase of the CCD plane'); %ribui

%--Recording

%   Creating b plane meshgrid 
Wb=deltaX*M; % Width of the mesgrid
n=deltaY*(1:N)-Wb/2; % shifting Y Axis to the center
m=deltaX*(1:M)-Wb/2; % shifting Y Axis to the center
[X,Y]=meshgrid(m,n); % Creating a matrix
Q_obj=exp((pi*1i)/(lambda*d)*(X.^2+Y.^2));%phase calculation

U_obj=Q_CCD.*fftshift(fft2(fftshift(object.*Q_obj)));
Ref_beam=20*exp(((2.*pi.*1i)/(lambda)).*(d+theta_x.*X_ccd+theta_y.*Y_ccd));

%Hologram Intensety 
H_tot = (abs(U_obj+Ref_beam)).^2;
figure;
imagesc(H_tot); 
title('Intensity of Hologram-Plane'); 
axis on;

%Reconstruction 
Q_CCD=exp((1i*pi)/(lambda*d)*(X_ccd.^2+Y_ccd.^2)); %phase calculation

U_rec = Q_obj.*fftshift(fft2(fftshift(H_tot.*conj(Ref_beam).*Q_CCD))); 
U_rec_abs=abs(U_rec);
figure;
imshow(10*mat2gray(U_rec_abs)); 
title('Reconstruction'); 
axis on;

