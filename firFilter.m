% FIR Filter

clear
clc

Fs = 44100;
fp = 440;
fs = 2*fp;
W = fs; % cutoff
N = 64; % order
type = "low";
B = fir1(N, W, type);
stem(B);
plot(B)


