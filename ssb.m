% ssb modulation
close all
clear
clc

pkg load signal % load signal processing package

N = 128;
n = -N/2:1:N/2;
h = 1./(pi*n) .* (1-cos(pi*n));
h(N/2+1) = 0; % div. by zero


w = 0.8*hanning(N+1)';
w = 0.8*hamming(N+1)';
w = 0.45 * (1 - cos(2*pi*[0:(N)]/(N)));

h = h.*w;

Fs = 44100;
dur = 5;
t = 0:1/Fs:dur-1/Fs;
fmod = 440;
m1 = sin(2*pi*fmod*t);
[x, Fs] = audioread('batman.wav',[1, length(t)]);
x = x';

m2 = filter(h, 1, m1); % quadratur signal
m1 = [zeros(1,N/2) m1(1:length(m1)-N/2)]; % in phase signal needs to be delayed by N/2 

y1 = 0.7*(1+m1).*x;
y2 = 0.7*(1-m2).*x;

yL = y1-y2;
yR = y1+y2;

figure(1);
subplot(211), plot(x, 'm');
xlim([23e3 24e3]);
subplot(212),
plot(yL, 'b');
hold on
plot(yR, 'r--');
hold off
xlim([23e3 24e3]);

%sound(x, Fs);
pause;
y = zeros(2,length(yR));
y(1,:) = yL;
y(2,:) = yR;
sound(0.9*y, Fs);
