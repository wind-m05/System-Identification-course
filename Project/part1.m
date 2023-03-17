% Project part 1 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc
% Plot options
options.subplot = true; % true,false
options.xscale = 'lin'; % log,lin
options.yscale = 'mag'; % mag,db
options.plot = 'line'; %scatter,line

%% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.0001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
x = freq_query/pi;

figure()
[mag,phase_var,freq,l1,l2,p1,p2] = bp(F,options,freq);
xlabel(p1, 'Frequency [rad/s]');
ylabel(p1, 'Magnitude [dB]')
xlabel(p2, 'Frequency [rad/s]');
ylabel(p2, 'Degrees [^\circ]','Interpreter','tex')
l1.LineWidth = 2; l2.LineWidth = 2;

%% Test saturation function
x = linspace(0,100);
r = 1000*sin(x); % Iterate
[u,y] = assignment_sys_20(r);
grad_u = gradient(u);
sat_query = find(grad_u == 0);
[Mmin,Mmax] = bounds(u(sat_query));

