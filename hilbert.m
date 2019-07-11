% hilbert filter 
close all
clear
clc

pkg load signal % load signal processing package

N = 128;
n = -N/2:1:N/2;
h = 1./(pi*n) .* (1-cos(pi*n));
h(N/2+1) = 0; % div. by zero


w = hanning(N+1)';
w = hamming(N+1)';
%w = 0.5 * (1 - cos(2*pi*[0:(N)]/(N)));

h = h.*w;
figure(2)
subplot(211), stem(h);
[H, omega] = freqz(h, 1, 512);
%H = abs(fft(h, 512))
subplot(212), plot(omega, abs(H));


Fs = 44100;
dur = .1;
t = 0:1/Fs:dur-1/Fs;
x = sin(2*pi*200*t);
x = chirp(t, 1000, dur, 1500);
%x = ones(1,length(x));
%x = [1 x];
% something seems to be wrong with filter... no, input signal was low freq
y = filter(h, 1, x); % quadratur signal

%% fir filtering
%y_ = [];
%buffer = zeros(1,N);
%x = [buffer x];
%for i=1:length(x)-N-1
%  y_(i) = sum(h(N:-1:1).*buffer) + h(N+1)*x(N+i);
%  buffer = x(i:i+N-1);
%endfor

figure(1);
plot(x, 'b');
hold on
plot(y, 'r');
%plot(y_, 'm');
hold off

z = zeros(1,N/2); % in phase signal needs to be delayed by N/2 
x = [z x];
