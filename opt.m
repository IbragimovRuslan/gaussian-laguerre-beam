%https://www.mathworks.com/matlabcentral/fileexchange/15459-basic-paraxial-optics-toolkit/content/transverse/LaguerreGaussianE.m

clear all;

w=[0.001; 0.001];
rseed=[0*max(w):max(w)/30:3*max(w)];
thetaseed=[0:360]*pi/180;
[r,theta]=meshgrid(rseed,thetaseed);
lambda = [1.064e-6 ; 1.064e-6];
R = [-30 ; -30];
q = (1./R - i* lambda./pi./w.^2).^(-1); 
a=[1;1];
p=[1;0]; m=[1;0];
E=LaguerreGaussianE([p,m,q,lambda,a],r,theta)+LaguerreGaussianE([p,-m,q,lambda,a],r,theta);
[x,y]=pol2cart(theta,r); colormap(bone);
GL = E(:,:,1);  %Init Gausse-Lagger

%hank = besselj(1,E(:,:,1));
%mesh(r,theta,reshape(abs(hank),size(r)));
%axis vis3d;


figure(1)
subplot(1,1,1); h1=pcolor(x,y,real(E(:,:,1)));
set(h1,'EdgeColor','none'); axis square; colormap(bone);
title('Gausse-Laguerr');
figure(2)
subplot(1,1,1); h1=pcolor(x,y,angle(E(:,:,1)));
set(h1,'EdgeColor','none'); axis square; colormap(bone);
title('Phase');
figure(3)
subplot(1,1,1); h1=pcolor(x,y,abs(E(:,:,1)));
set(h1,'EdgeColor','none'); axis square; colormap(bone);
title('Amplitude');
figure(4)
subplot(1,1,1); h1=pcolor(x,y,abs(E(:,:,1).^2));
set(h1,'EdgeColor','none'); axis square; colormap(bone);
title('Intensity');

%{
%subplot(2,1,1); h1=pcolor(x,y,real(E(:,:,1)));  set(h1,'EdgeColor','none'); axis square;
%set(h1,'EdgeColor','none'); axis square;



y1=fftshift(fft(GL)); 
N=length(y);         
n=-(N-1)/2:(N-1)/2;  
f=sqrt(y1.*conj(y1));  
fdouble = log(double(f));
figure(4);
subplot(1,1,1); h1=pcolor(x,y,(fdouble));
set(h1,'EdgeColor','none'); axis square;
title('Fourier transform of Gausse-Lagger beam');


figure(5);
plot(n,fdouble);
title('FFT GL Using plot');

fftOriginal = fft2(double(GL));
shiftedFFT = log(double(fftshift(fftOriginal)));
figure(6);
plot3(x,y,(shiftedFFT));
title('Real Part of Spectrum');

figure(7);
mesh(real(shiftedFFT));

figure(8);
plot3(x,y,(fdouble));
title('GL using plot3');


NFFT=1024; %NFFT-point DFT	 	 
X=fftshift(fft(GL,NFFT)); %compute DFT using FFT	 	 
fVals=(-NFFT/2:NFFT/2-1)/NFFT; %DFT Sample points	 	 
plot(fVals,abs(X));	  	 
title('Double Sided FFT');	 	 
xlabel('Sample points (N-point DFT)')	 	 
ylabel('DFT Values');

%}

%% check FFT transform of rectangle.
hw=10;

dt=.001;
t=[-60:dt:60];
           
%% ractangular pulse.
                                            
x1=(5/2)*(sign(t+hw)-sign(t-hw));
figure(5);
plot(t,x1);
title(['Rectangular pulse width-',num2str(2.*hw),'ms']);
xlabel('time(ms)');
ylabel('Amplitude(V)');
rge=40;
axis([-rge rge 0 6]);
%pause


%% rectangular pulse frequency content by fourier analysis.

y2=fftshift(fft(x1));  % moving the zero-frequency component to the center of the array
N2=length(y2);         %to take the frquecny axis of the hoarmonics.
n2=-(N2-1)/2:(N2-1)/2;  %divide the frequency compone
f2=sqrt(y2.*conj(y2)); % to take the amplitude of each hoarmony. 
figure(6);
plot(n2,f2);  
axis([-120 120 0 120000]);  
title(['Rectangular Pulse(width-',num2str(hw),'ms',')  frequency distribution']);
xlabel('frequency(Hz) ');
ylabel('Amplitude');


fs=80; %sampling frequency
L=length(GL);
NFFT = 1024;
X = fftshift(fft(GL,NFFT));
Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector
figure(7);
plot(f,abs(X)/(L),'r');
title('Magnitude of FFT');
xlabel('Frequency (Hz)')
ylabel('Magnitude |X(f)|');
xlim([-10 10])



L=length(GL);
res = BPF_2D(GL);
r = real((res));
figure(8)
subplot(1,1,1); h1=pcolor(x,y,r);
set(h1,'EdgeColor','none'); axis square; colormap(bone);
title('FFT');

figure(9)
colormap(bone);
mesh(abs(r));
title('FFT-MESH');

figure(10)
colormap(bone);
mesh(abs(GL));


fftOriginal = fft2(GL);
shiftedFFT = real(fftOriginal);
subplot(1,1,1); h1=pcolor(x,y,(shiftedFFT));
set(h1,'EdgeColor','none'); axis square;
title('Gausse-Laguerr');
