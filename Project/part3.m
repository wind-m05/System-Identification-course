% Project part 3 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc

% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
%% PRBS1 vs PRBS2 vs Gaussian noise
N = 3000;
p = 0.6;
Range = [-3,3];
Band = [0 bandpass_freq];
index_min = fix(N/2 - N/2*bandpass_freq);
index_plus = fix(N/2 + N/2*bandpass_freq);
[r_1,freq] = idinput(N,'prbs',Band,Range);% Problem is that the spectrum does not account for the fact that the passband is at 0.667, so there is no extra power pushed into the lower frequencies.
r_2 = prbs(N,max(Range),p);
r_3 = idinput(N,'rgs',Band,Range);
stem(abs(fftshift(fft(r_1))))
hold on
stem(abs(fftshift(fft(r_2))))
stem(abs(fftshift(fft(r_3))))
xline(index_min,'LineWidth',2)
xline(index_plus,'LineWidth',2)
legend('R_1','R_2','R_3')
xlabel('Frequency index [-]')
ylabel('Amplitude [-]')

[u_1,y_1] = assignment_sys_20(r_1);
[u_2,y_2] = assignment_sys_20(r_2);
[u_3,y_3] = assignment_sys_20(r_3);
% Filtered gaussian noise
[phi_u_1,freq1] = pwelch(u_1);
[phi_u_2,freq2] = pwelch(u_2);
[phi_u_3,freq3] = pwelch(u_3);
sys_u1 = frd(phi_u_1,freq1);
sys_u2 = frd(phi_u_2,freq2);
sys_u3 = frd(phi_u_3,freq3);
figure
bodemag(sys_u1)
hold on
bodemag(sys_u2)
bodemag(sys_u3)
legend('U_1','U_2','U_3')
grid on
% Trapizoidal integration
u1_int = cumsum(phi_u_1);
u2_int = cumsum(phi_u_2);
u3_int = cumsum(phi_u_3);
figure
plot(u1_int,'LineWidth',2)
hold on
plot(u2_int,'LineWidth',2)
plot(u3_int,'LineWidth',2)
legend('U_1','U_2','U_3')
grid on
xlabel('Frequency index [-]')
ylabel('Summed Power [W]')
% Calculate power loss
[phi_r_1,freq1] = pwelch(r_1);
[phi_r_2,freq2] = pwelch(r_2);
[phi_r_3,freq3] = pwelch(r_3);
phi_dif_1 = phi_r_1 - phi_u_1;
phi_dif_2 = phi_r_2 - phi_u_2;
phi_dif_3 = phi_r_3 - phi_u_3;
dif_1_int = cumsum(phi_dif_1);
dif_2_int = cumsum(phi_dif_2);
dif_3_int = cumsum(phi_dif_3);
figure
plot(dif_1_int,'LineWidth',2)
hold on
plot(dif_2_int,'LineWidth',2)
plot(dif_3_int,'LineWidth',2)
grid on
legend('U_1','U_2','U_3')
xlabel('Frequency index [-]')
ylabel('Summed Power [W]')
% Time domain behaviour
figure
stem(r_1)
hold on
stem(r_2)
stem(r_3)

