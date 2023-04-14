% Project part 5 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc
% font = 12;
% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
% 
N = 1000;
p = 0.6; % probability of switching in prbs.
Range = [-3,3];
Band = [0 bandpass_freq];
% index_min = fix(N/2 - N/2*bandpass_freq);
% index_plus = fix(N/2 + N/2*bandpass_freq);

% Monte carlo runs 
nb = 5;
nf = 5;
nk = 1;
sims = 100;
model_runs = cell(sims,1);
for i = 1:sims
r = idinput(N,'rbs',Band,Range);
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
sys = oe(data,[nb nf nk],oeOptions('Focus','prediction')); % OE
model_runs{i} = sys;
end
save('OE_montecarlo','model_runs')
%% Further analysis
load('OE_montecarlo')
for i = 1:length(model_runs)
pzmap(model_runs{i})
hold on
end
