% Project part 4 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc

% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
% 
N = 3000;
p = 0.6;
Range = [-3,3];
Band = [0 bandpass_freq];
index_min = fix(N/2 - N/2*bandpass_freq);
index_plus = fix(N/2 + N/2*bandpass_freq);
r = prbs(N,max(Range),p);
stem(abs(fftshift(fft(r))))
hold on
xline(index_min,'LineWidth',2)
xline(index_plus,'LineWidth',2)
xlabel('Frequency index [-]')
ylabel('Amplitude [-]')
%% Nonparametric ID
close all
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
na = 4;
nb = 4;
nk = 3;

sys1 = arx(data,[na nb nk],arxOptions('Focus','prediction'));
sys2 = oe(data,[na nb nk],arxOptions('Focus','prediction'));
figure
resid([y u],sys1)
figure
resid([y u],sys2)
figure
bode(sys1)
hold on
bode(sys2)
legend('ARX','OE')


