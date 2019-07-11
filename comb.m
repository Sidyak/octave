% comb filter
close all
clear
clc

pkg load signal % load signal processing package

N = 3;
h = zeros(1,N+1);
h(1) = 1; 
%h(N/2) = 0.5; 
%h(N/2-2) = 0.25; 
h(N+1) = 1; 


w = 0.8*hanning(N+1)';
w = 0.8*hamming(N+1)';
w = 0.45 * (1 - cos(2*pi*[0:(N)]/(N)));

%h = h.*w;

%n = 0:N;
%c = cos(pi/N * (n +0.5+N/2) * (k+0.5);)

Fs = 44100;
dur = 1;
t = 0:1/Fs:dur-1/Fs;
fmod = 2e3;
%m = sin(2*pi*fmod*t);
x = 0.7*sin(2*pi*440*t);

%[x, Fs] = audioread('batman.wav',[1, length(t)]);
%x = x';

y = filter(h, 1, x); % quadratur signal

% modulation by e^j*pi/10
F=16;
t1 = zeros(F*N+1,N+1);
for f = 0:F*N
  t1(f+1,:) = h .* exp(+1i*2*pi*f/(F*N).*[0:N]/(N));
end
figure(2);
cmap = hsv(F*N+1);
plot(20*log10(abs(fft(t1(1,:), 1024))),'-','Color',cmap(1,:));
legend('1st');
hold on
li = ['1' ];
for f=1:F*N
  plot(20*log10(abs(fft(randn(1)*t1(f+1,:), 1024))),'-','Color',cmap(f+1,:));
  li = [li, sprintf('%d', f+1)];
end
legend();
hold off

