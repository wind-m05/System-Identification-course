% Project part 2 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc
% Plot options
options.subplot = true; % true,false
options.xscale = 'lin'; % log,lin
options.yscale = 'mag'; % mag,db
options.plot = 'line'; %scatter,line

%% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;

%% ETFE white noise test
N = 1024;
trials = 10;
a = -2.8;
b = 2.8;
r = (b-a).*rand(N,1) + a;
[u,y] = assignment_sys_20(r);
data = iddata(y,u,1);
% ir = cra(data,length(data.u)-1);
M = 30;

% Many trials
figure()
for i = 1:trials
r = (b-a).*rand(N,1) + a;
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
respetfe = etfe(data,M);
% bp(respetfe,options,respetfe.frequency);
bode(respetfe)
hold on
end

%% Non-parametric identification 128 specific frequencies
T = 1024; % Period
Tnum = 1; % Number of periods
range = [-3 3]; % Range in time domain
band = [0 bandpass_freq]; % Frequency band
sinnum = 128; % Number of sinusoids
NumTrials = 10;
GridSkip = 1;
SineData = [sinnum,NumTrials,GridSkip];
[r,freq] = idinput([T 1 Tnum],'sine',band,range,SineData);
figure
stem(r)
fourier_grid = 2*pi/T*(1:GridSkip:fix(T/2));
fourier_grid_norm = fourier_grid/pi;
fourier_grid_red = fourier_grid_norm(fourier_grid_norm <= bandpass_freq);

rfft = fft(r);
figure
stem(abs(rfft))
%%
uabs = abs(rfft);
uabsn = uabs/max(uabs);
numsinu = sum(uabsn)/2;
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
% ir = cra(data,length(data.u)-1);
M = T; % No averaging

[u,y] = assignment_sys_20(r);
data = iddata(y,u);
respetfe = etfe(data);
respetfe_summed = respetfe*0;
clear respetfe
% Many trials
trials = 1;
figure(1)
for i = 1:trials
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
respetfe{i} = etfe(data,M); % MUST SPECIFY WHICH freq !
bode(respetfe{i},'r.')
respetfe_summed = respetfe_summed+respetfe{i};
hold on
    if i == trials
        hold on
        respetfe_expected = respetfe_summed/trials;
        bode(respetfe_expected,'b*')
    end
end
%% Just some test, appears to be unnecessesary
response = zeros(1,length(respetfe{1}.Frequency));
for i = 1:length(respetfe{1}.Frequency)
    for k = 1:length(respetfe)
        response(i) = response(i) + squeeze(respetfe{k}.ResponseData(i));  
    end
end
figure(1)
hold on
response = response/trials;
sys = frd(response,respetfe{1}.Frequency);
bode(sys,'b*')
%% Task 2.3 

%% Calculate estimate of Phi_v via E
respetfe_gem = respetfe{1}*0;
for i = 1:trials
respetfe_gem = respetfe_gem + respetfe{i};
end
respetfe_gem = respetfe_gem/trials;

respetfe_var = respetfe{1}*0;
for i = 1:trials
respetfe_var = respetfe_var + abs(respetfe{i} - respetfe_gem)^2;
end
respetfe_var = respetfe_var/trials;
figure
Phi_v_est = sum(abs(fft(u)).^2)/length(u) * respetfe_var;
bp(Phi_v_est,options)

%% Calculate estimate of Phi_v via spectral densities
% ufft = fft(u);
% yfft = fft(y);
% phi_u = (abs(ufft).^2)/length(ufft);
% phi_y = (abs(yfft).^2)/length(yfft);
phi_u = pwelch(u);
phi_y = pwelch(y);
phi_uy = cpsd(y,u);
phi_v_est_spect = phi_y - (abs(phi_uy).^2)./(phi_u);
figure
stem(phi_v_est_spect)

%% r = 0 
close all
N=1024;
r = zeros(1,N);
[u,y] = assignment_sys_20(r);
data = iddata(y,u,1);
% ir = cra(data,length(data.u)-1);
nfft = N/4;
window  = hann(nfft);
figure
pwelch(y,window,nfft/2,nfft)
figure
window  = hann(N);
nfft = N/64;
periodogram(y,window,nfft)
figure
stem(abs(fft(y)))
