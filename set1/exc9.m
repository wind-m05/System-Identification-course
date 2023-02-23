clear all, close all, clc
options.subplot = false;
options.unwrap = false;
options.scale = 'db';
options.plot = 'line';
Ts = 1;
data = load('batch');
data = iddata(data.y,data.u,Ts);

% Analyze the correlation for window size
ir = cra(data,length(data.u)-1);
% ETFE
M = 140;
respetfe = etfe(data,M);

% SPA 
G = spa(data,M);

% TFestimate
nfft = length(data.u)/1;
window = hann(nfft);
data = load('batch');
[Txy,w] = tfestimate(data.u,data.y,window,nfft/2,nfft,1/Ts);
tfest = frd(Txy,w*2*pi);

%% Plot everything
figure()
bp(respetfe,options);
hold on
bp(G,options);
bp(tfest,options);
grid on
legend('etfe','spa','tfestimate')