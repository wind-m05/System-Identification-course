%% Box jenkins
clear all, close all
load('median_F')
load('median_B')
BJ_old = load('BJ_working_53351');
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
N = 1000;
Range = [-3,3];
Band = [0 bandpass_freq];

r = idinput(N,'rbs',Band,Range);
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
M_init = idpoly([],median_B,BJ_old.BJ.C,BJ_old.BJ.D,median_F);
% Find median box jenkins for accurate initial noise model parameters
nb = 5; % From OE
nc = 2;
nd = 3;
nf = 5; % From OE
nk = 1;
M_init = idpoly([],median_B,[1 -3 3],[1 -3 3 1],median_F);
BJ = bj(data,M_init);
% BJ = bj(data,M_init);
% BJ = oe(data,[nb nf nk]);
figure
resid([y u],BJ)
% save('BJ_working_5325','BJ')
save('BJ_working_52351','BJ')