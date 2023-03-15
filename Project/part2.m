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
freq = (0:0.01:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
x = freq_query/pi;

% ETFE white noise test
N = 1024;
trials = 10;
a = -2.8;
b = 2.8;
r = (b-a).*rand(N,1) + a;
[u,y] = assignment_sys_20(r);
data = iddata(y,u,1);
ir = cra(data,length(data.u)-1);
M = 30;

%% Many trials
figure()
for i = 1:trials
r = (b-a).*rand(N,1) + a;
[u,y] = assignment_sys_20(r);
data = iddata(y,u,1);
respetfe = etfe(data,M);
% bp(respetfe,options,respetfe.frequency);
bode(respetfe)
hold on
end

%% Non-parametric identification 128 specific frequencies
% ETFE frequency input test
% N = 1024;
% bandpass = x;
% amp = 35;
% 
% kgrid = 0:N;
% fgrid = (-pi:2*pi/N:pi);
% for i = 1:N+1
%     if mod(kgrid(i),8) == 0
%         if abs(fgrid(i)) <= bandpass
%             u(i) = amp;
%         end
%     else
%     u(i) = 0;
%     end
% end
% r = ifft(u);
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u,1);
% ir = cra(data,length(data.u)-1);
% M = 40;

% Many trials
% figure()
% for i = 1:trials
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u,1);
% respetfe = etfe(data,M);
% % bp(respetfe,options,respetfe.frequency);
% bode(respetfe)
% hold on
% end


%% Plotting

% figure
% stem(u)
% title('u')
% figure
% plot(r)
% title('r')
% 
% figure
% u = fft(r);
% stem(u)
% title('Input spectrum')