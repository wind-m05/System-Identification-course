clear all, close all, clc
options.subplot = false;
options.unwrap = false;
options.scale = 'db';
options.plot = 'line';

G0 = tf([1 -0.5 0.2 0.8 0.1], [1 -1.7 1.6 -0.8 0.25], 1);
H0 = tf(1, [1 -1.7 1.6 -0.8 0.25], 1);
N = 1024;
u = sign(randn(N, 1));
lambda = sqrt(0.1);
e = lambda * randn(N,1);
y = lsim(G0, u) + lsim(H0, e);
data = iddata(y,u,1);
% Analyze the correlation for window size
ir = cra(data,length(data.u)-1);
% ETFE
M = 40;
respetfe = etfe(data,M);
figure()
bp(respetfe,options,respetfe.frequency)
hold on
bp(G0,options,respetfe.frequency)
legend('ETFE','Actual plant G0')
