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
na = 5;
nb = 5;
nk = 1;
% Most likely a 5th order OE model...
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

% Which model order to take:
% V = arxstruc(ze,zv,NN) !!
% Validation
%- check with respect to ETFE, but be careful since ETFE can also be wrong
%- Check pole zero map for maybe pole zero cancelations which will imply a
%too high model order.
%-residual test data

%% Step response of the system
close all
N = 3000;
r = zeros(N,1);
trials = 10000;
% z = zeros(N/2,1);
% o = ones(N/2,1);
r(1:N/2,1) = 0;
r(N/2:end,1) = 1;
y_sum = zeros(N,1);
for i = 1:trials
[u,y] = assignment_sys_20(r);
y_sum = y_sum + y;
end
y_avg = y_sum/trials;
figure
stem(y_avg)
% So the delay is nk = 1;
%% Box jenkins

nb = 5;
nc = 5;
nd = 7;
nf = 7;
nk = 1;

BJ = bj(data,[nb nc nd nf nk]);
figure
resid([y u],BJ)