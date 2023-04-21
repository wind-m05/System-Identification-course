clear all
close all

load sysOE.mat

z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
T = 1024; % Period
Tnum = 1; % Number of periods           %input length = T*Tnum=1024
range = [-3 3]; % Range in time domain
band = [0 bandpass_freq]; % Frequency band
sinnum = 128; % Number of sinusoids
NumTrials = 10;
GridSkip = 2;
SineData = [sinnum,NumTrials,GridSkip];
[r,freq] = idinput([T 1 Tnum],'sine',band,range,SineData);

[u,y] = assignment_sys_20(r);

simout=sim(sys1,u);

font = 16;

s = get(0, 'ScreenSize');
figure('Position', [10 50 600 425]);
plot(u,LineWidth=2)
xlim([1 100])
xlabel('Time [s]',FontSize=font,Interpreter='latex')
ylabel('Input $u(t)$',FontSize=font,Interpreter='latex')

figure('Position', [10 50 600 425]);
plot(y,LineWidth=2)             %Actual system
hold on;
plot(simout,LineWidth=2)        %OE model
xlim([1 100])
legend('Actual system','OE model',FontSize=font,Interpreter='latex')
xlabel('Time [s]',FontSize=font,Interpreter='latex')
ylabel('Output $y(t)$',FontSize=font,Interpreter='latex')

figure('Position', [10 50 600 425]);
plot(y-simout,LineWidth=2)             %error
xlim([1 100])
xlabel('Time [s]',FontSize=font,Interpreter='latex')
ylabel('Simulation error $e_{sim}(t)$',FontSize=font,Interpreter='latex')

%% TFestimate only
%close all
font = 16;
figure('Position', [10 50 600 425]);
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
%xlabel('Frequency [rad/s]',FontSize=font,Interpreter='latex')
hold on
grid on
%     if i == trials
%         hold on
%         respetfe_expected = respetfe_summed/trials;
%         bode(respetfe_expected,'b*')
%     end
end
bode(sys1,'r')
%bode(sys1,'Color',"#D95319")