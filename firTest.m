% FIR Filter

clear

Fs = 44100;
fp = 440;
fs = Fs/2/4;
W = fs/Fs; % cutoff
N = 64; % order
type = "low";
B = fir1(N, W, type);
subplot(211), stem(B);
[H, W2] = freqz(B, 1, N);
subplot(212), plot(W2/(2*pi), 20*log10(abs(H)))

