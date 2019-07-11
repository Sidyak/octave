% stft

clear
clc
format short g

N = 2048;
Ha = N/2;     % TODO: fix for other hopsizes (e.g.) N/4
Hs = round(N/4);

Fs = 44100;
len = 7*Fs;
t = 0:len;
x = 0.5*sin(2*pi*440*t/Fs);
%x = randn(1, len);
%x = ones(1,len);
[x, Fs] = audioread('batman.wav',[1, len]);
x = x';

tstretch = Ha/Hs
Fsout = Fs * tstretch
x = x / max(x); % normalize to 1
w = hanning(N);
w = hamming(N);
w = 0.5 * (1 - cos(2*pi*[0:(N-1)]/(N-1)))';

ts = Ha/Hs;
omega = 2*pi*Hs*[0:N-1]/N;

x_ = zeros(1,len);
y1_ = zeros(1,len);
x1 = zeros(1,N);
Y1 = zeros(1,N);
Y1_ = zeros(1,N);
phi = zeros(1,N);
a = zeros(1,N);
deltaPhi = zeros(1,N);
phi0 = zeros(1,N);
psi = zeros(1,N);

for n=0:len/max(Ha,Hs)-2
  lo = n*Ha+1;
  up = n*Ha+N;
  x1 = x(lo:up); % frame with overlap of hopsize
  y1 = x1 .* w'; % windowing
  Y1 = fft(fftshift(y1), N);
  % phase processing / phase shift or time stretch
  phi = arg(Y1); %atan2(imag(Y1), real(Y1)); % 
  a = abs(Y1); % sqrt(real(Y1).*real(Y1) + imag(Y1).*imag(Y1)); % 
  deltaPhi = omega + princarg(phi-phi0-omega);
  phi0 = phi;
  psi = princarg(psi+deltaPhi*tstretch); 
  Y1_ = a .* exp(1i*(psi));
  % if tstretch is an integer, unwrapping is no issue and 
  % we can directly multiply phi with the stretch factor
  %Y1_ = a .* exp(1i*(phi*tstretch));
  % done phase processing
  if(tstretch == 1)
    y1_ = fftshift(real(ifft(Y1_, N)));
  else
    y1_ = fftshift(real(ifft(Y1_, N))) .* w';
  end  
  lo = n*Hs+1;
  up = n*Hs+N;
  x_(lo:up) = x_(lo:up) + y1_;
end
%figure(1),
%subplot(311), stem(x(Ha+1:len-Ha));
%subplot(312), stem(x_(Ha+1:len-Ha));
%subplot(313), stem(abs(x(Ha+1:len-Ha)-x_(Ha+1:len-Ha)));

figure(2),
subplot(211), plot(abs(fft(fftshift(x))),'b'); 
hold on;
plot(abs(fft(fftshift(x_))),'r');
hold off;
%xlim([0 4e3])
subplot(212), plot(phi);
hold on;
plot(psi,'r');
hold off;

%sound(x, Fs);
pause;
sound(x_, Fs);%*Hs/Ha);

