clc;clear;close all
%%
file = load('ia.mat');   % load the dataset for Fig.7
ia = file.ia;
N = length(ia);
%%
ia = ia - mean(ia);
ia = ia ./ max(ia);
Fs = 1000;
%% FFT
ia1 = hilbert(ia);
f_ia = fft(ia1);
f_ia = abs(f_ia) / max( abs(f_ia) );
log_f = 10*log10( f_ia );
fs = ( 0:length(ia)-1 )*Fs / length(ia1);
%% ESPRIT
range = [40, 80];
G = 80;
[hat_f,hat_s]=ESPRIT_SelectOrder(ia,G,Fs,range);
hat_s = abs(hat_s) / max(abs(hat_s));
%% TOEO
i_toeo = TOEO(ia);
f_toeo = fft(i_toeo);
f_toeo = abs(f_toeo) / max(abs(f_toeo));
%% our method
Fc = 60;
baseline = cos( 1/Fs*2*pi*(1:N)*Fc );
ia = ia .* baseline.';
ia = hilbert(ia) - mean( hilbert(ia) );
ia = ia ./ max( abs(ia) );
r = 0.2;
f_sample = Fc + (40:r:80);
N_all = (0:1:length(ia)-1)';
[mu1,delta_invv1,f_sample_new1] = VBI_offgrid_CGDP(ia,N_all,f_sample,Fs);
mu1 = abs(mu1) / max( abs(mu1) );
log_mu = 10*log10( abs(mu1) );
%%
figure
subplot(2,2,1)
plot(fs, log_f + 40)
axis([40,80,5,50])
xlabel('Frequency [Hz]')
ylabel('Amp. [dB]')
title('FFT')

subplot(2,2,2)
stem( hat_f, 10*log10(hat_s) + 40 )
axis([40,80,5,50])
xlabel('Frequency [Hz]')
ylabel('Amp. [dB]')
title('ESPRIT')

subplot(2,2,3)
plot( fs(1:end-3), 10*log10(abs(f_toeo))+30 )
axis([0,20,5,50])
xlabel('Frequency [Hz]')
ylabel('Amp. [dB]')
title('TOEO')

subplot(2,2,4)
stem(f_sample_new1 - Fc, log_mu + 40)
axis([40,80,5,50])
xlabel('Frequency [Hz]')
ylabel('Amp. [dB]')
title('Our method')
