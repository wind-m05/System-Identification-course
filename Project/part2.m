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

%% ETFE white noise test (Initial trial)
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

%% Reference definition
T = 1024; % Period
Tnum = 1; % Number of periods
range = [-3 3]; % Range in time domain
band = [0 bandpass_freq]; % Frequency band
sinnum = 128; % Number of sinusoids
NumTrials = 10;
GridSkip = 2;
SineData = [sinnum,NumTrials,GridSkip];
[r,freq] = idinput([T 1 Tnum],'sine',band,range,SineData);
fourier_grid = 2*pi/T*(1:GridSkip:fix(T/2));
fourier_grid_norm = fourier_grid/pi;
fourier_grid_red = fourier_grid_norm(fourier_grid_norm <= bandpass_freq);

[u,y] = assignment_sys_20(r);
data = iddata(y,u);
% ir = cra(data,length(data.u)-1);
M = T; % No averaging
% Just one trial to get the datastructure of the summed etfe right
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
respetfe = etfe(data);
respetfe_summed = respetfe*0;
clear respetfe


%% TF estimate vs ETFE vs SPA test
T = 1024; % Period
window = rectwin(T);
Tnum = 1; % Number of periods
range = [-3 3]; % Range in time domain
band = [0 bandpass_freq]; % Frequency band
sinnum = 128; % Number of sinusoids
NumTrials = 10;
GridSkip = 2;
SineData = [sinnum,NumTrials,GridSkip];
[r,freq] = idinput([T 1 Tnum],'sine',band,range,SineData);
% freq = [1:128]/128*pi;
[u,y] = assignment_sys_20(r);
M = T; % No averaging
data = iddata(y,u);
respetfe = etfe(data);
options.subplot = true; % true,false
options.xscale = 'lin'; % log,lin
options.yscale = 'mag'; % mag,db
options.plot = 'line'; %scatter,line
[resptfest,freq] = tfestimate(u,y,window,[],freq);
respetfe = etfe(data,M);
respspa = spa(data,T,freq);
sys = frd(resptfest,freq);
bode(sys,'b*')
hold on
bode(respetfe,'r*')
bode(respspa,'y*')

[pxx,f] = periodogram(u);
f_grid = f(pxx >=0.001); % this is the same as freq, which implies the same frequencies in u are present as in r.
% This also implies that the ETFE MUST be able to tell which frequencies it
% is getting, so why does it not only display and estimate the frequencies
% that I am actually supplying to the system?
%% ETFE only
trials = 10;
respetfe=[];
for i = 1:trials
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
respetfe{i} = etfe(data,M);
freqresp(respetfe{i},freq);
bode(respetfe{i},'r.');
respetfe_summed = respetfe_summed+respetfe{i};
hold on
%     if i == trials
%         hold on
%         respetfe_expected = respetfe_summed/trials;
%         bode(respetfe_expected,'b*')
%     end
end
grid on
%% Test to see if you can sum up the data structure idfrd
% response = zeros(1,length(respetfe{1}.Frequency));
% for i = 1:length(respetfe{1}.Frequency)
%     for k = 1:length(respetfe)
%         response(i) = response(i) + squeeze(respetfe{k}.ResponseData(i));  
%     end
% end
% figure(1)
% hold on
% response = response/trials;
% sys = frd(response,respetfe{1}.Frequency);
% bode(sys,'b*')

%% TFestimate only
trials = 10;
T = 1024; % Period
resptfest = cell(trials,1);
window = rectwin(T);
sys = [];
for i = 1:trials
[u,y] = assignment_sys_20(r);
[resptfest{i},freq] = tfestimate(u,y,window,0,freq);
sys{i} = frd(resptfest{i},freq);
bode(sys{i},'b.')
hold on
%     if i == trials
%         hold on
%         respetfe_expected = respetfe_summed/trials;
%         bode(respetfe_expected,'b*')
%     end
end
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
% Tried functions: xcorr, periodogram, pwelch, cpsd
% freq = [0:0.01:pi]
close all
N = length(u);
nfft = N/4;
window_ = hamming(nfft);
phi_u = pwelch(u,window_,nfft/2,freq);
phi_y = pwelch(y,window_,nfft/2,freq);
[phi_uy,freqcpsd] = cpsd(y,u,window_,nfft/2,freq);
% phi_uy = sqrt(phi_y.*conj(phi_u));
phi_v_est_spect = phi_y - (abs(phi_uy).^2)./(phi_u);
figure
stem(phi_v_est_spect)
% Now towards bode plot
semilogx(freqcpsd,mag2db(phi_v_est_spect),'b.','MarkerSize',10) 
grid on
xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')


%% r = 0 (not allowed I think)
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

%% Plot fourier grid
figure
plot(freq,'LineWidth',2)
xlabel('Index')
ylabel('Frequency [rad/s]')
grid on

%% Plot input signal r(t) u(t) and fourier transform
figure
subplot(211)
stem(r)
xlabel('Index')
ylabel('r(t) amplitude [-]')
grid on
subplot(212)
stem(abs(fftshift(fft(r))))
xlabel('Index')
ylabel('R(\omega) amplitude [-]')
grid on

%% Plot input signal u(t) and fourier transform
figure
subplot(211)
stem(u)
xlabel('Index')
ylabel('r(t) amplitude [-]')
grid on
subplot(212)
stem(abs(fftshift(fft(u))))
xlabel('Index')
ylabel('R(\omega) amplitude [-]')
grid on
