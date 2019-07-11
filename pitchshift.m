% pitchshift from dafx page 281

clear
clc
format short g

N = 2048;
Ha = round(N/8) % n2
Hs = round(N/4) % n1

Fs = 44100;
len = 1*Fs;
t = 0:len-1;
xin = 0.5*sin(2*pi*440*t/Fs);
%xin = randn(1, len);
%xin = ones(1,len);
%[xin, Fs] = audioread('batman.wav',[1, len]);
%xin = xin';

tstretch = Ha/Hs %(n2/n1)
resample =  1/tstretch
xin = xin / max(xin); % normalize to 1
w = hanning(N);
w = hamming(N);
w = 0.5 * (1 - cos(2*pi*[0:(N-1)]/(N-1)))';
omega = 2*pi*Hs*[0:N-1]/N;

x_ = zeros(1,len);
grain = zeros(1,len);
x1 = zeros(1,N);
Y1 = zeros(1,N);
Y1_ = zeros(1,N);
phi = zeros(1,N);
a = zeros(1,N);
deltaPhi = zeros(1,N);
phi0 = zeros(1,N);
psi = zeros(1,N);

% for linear interpolation of a grain of length N
% resampling factor Hs/Ha
lx = floor(N*resample);
x = 1+(0:lx-1)'*N/lx;
ix = floor(x) ;
ix1 = ix+1;
dx = x-ix;
dx1 = 1-dx;

pin = 0;
pout = 0;
pend = length(xin)-2*N;

while pin<pend
  x1 = xin(pin+1:pin+N); % frame with overlap of hopsize
  y1 = x1 .* w'; % windowing
  Y1 = fft(fftshift(y1), N);
  %----- phase processing / phase shift or time stretch
  phi = arg(Y1); %atan2(imag(Y1), real(Y1)); % 
  a = abs(Y1); % sqrt(real(Y1).*real(Y1) + imag(Y1).*imag(Y1)); % 
  deltaPhi = omega + princarg(phi-phi0-omega);
  phi0 = phi;
  psi = princarg(psi+deltaPhi*tstretch); 
  Y1_ = a .* exp(1i*(psi));
  
  if(tstretch == 1)
    grain = fftshift(real(ifft(Y1_, N)));
  else
    grain = fftshift(real(ifft(Y1_, N))) .* w';
  end  
  
  %----- interpolation (resampling)
  grain2 = [grain'; 0];
  grain3 = grain2(ix) .* dx1 + grain2(ix1) .* dx;
  x_(pout+1:pout+lx) = x_(pout+1:pout+lx) + grain3';
  pin = pin + Hs;
  pout = pout + Hs;
end
%figure(1),
%subplot(311), stem(xin(Ha+1:len-Ha));
%subplot(312), stem(x_(Ha+1:len-Ha));
%subplot(313), stem(abs(x(Ha+1:len-Ha)-x_(Ha+1:len-Ha)));
freq = -Fs/2:(length(xin))/Fs:Fs/2-(length(xin))/Fs;
freq = 0:(length(xin))/Fs:Fs-(length(xin))/Fs;
figure(2),
subplot(211), plot(freq, abs(fft(fftshift(xin))),'b'); 
hold on;
X_ = abs(fft(fftshift(x_)));
plot(freq, X_,'r');
hold off;
xlim([0 1e3])%length(X_)/4])
subplot(212), plot(phi);
hold on;
plot(psi,'r');
hold off;

%sound(xin, Fs);
pause;
sound(0.9*x_, Fs);%*Hs/Ha);

