% comb iir filter
close all
clear
clc

pkg load signal % load signal processing package

x=zeros(100,1);x(1)=1; % unit impulse signaolf length 100
Fs = 44100;
dur = 5;
len = dur*Fs;
t = 0:len-1;
[x, Fs] = audioread('batman.wav',[1, length(t)]);
x = x';

K=ceil(0.003*Fs); % set to 3 ms
modu = K * sin(2*pi*2*t/(Fs*dur));
figure(1);
plot(modu);

g=0.7;
b_0=0.7;
b_1=0.7;
a_1=0.7;
xhold=0;
yhold=0;
k = K;
Delayline=zeros(K,1); % memory allocation for length K

for n=1:length(x);
  k=floor(modu(n)/(K/13));
  yh(n)=b_0*Delayline(k)+b_1*xhold-a_1*yhold;
  % 1st-order difference equation
  yhold=yh(n);
  xhold=Delayline(k);
  y(n)=x(n)+g*yh(n);
  Delayline=[y(n); Delayline(1:k-1)];
end ;

figure(2);
subplot(211), plot(x);
subplot(212), plot(y);

sound(0.9*x, Fs);
pause;
sound(0.9*y, Fs);
