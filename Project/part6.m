%% Box jenkins
clear all, close all
font = 12;
load('median_F')
load('median_B')
BJ_old = load('BJ_working_52351');
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
% M_init = idpoly([],median_B,[1 -3 3],[1 3 -3 1],median_F);
% M_init = idpoly([],median_B,[1 1 3],[1 -3 -10],median_F);
% BJ = bj(data,M_init);
load('BJ_report')
figure
resid([y u],BJ)
% save('BJ_working_5325','BJ')
% save('BJ_working_52351','BJ')
% save('BJ_working_52251','BJ')
% save('BJ_report','BJ')

%% Variance analysis for BJ [5,2,3,5,1]
load('OE_var_report')
load('BJ_report')
P_rho_BJ = diag(getcov(BJ))';
BJ_var_G = [P_rho_BJ(1:5),P_rho_BJ(11:end)];

bar_aux_init = [];
for i = 1:length(BJ_var_G)
    bar_aux_init = [bar_aux_init; [theoretical_var(i),BJ_var_G(i)]];
end
figure
bar(bar_aux_init)
legend('$\sigma_{OE}^2\;$ Output error variance','$\sigma_{BJ}^2\;$ Box Jenkins variance',FontSize=font,Interpreter='latex')
xlabel('Parameters $[b_0,b_1,b_2,b_3,b_4,f_1,f_2,f_3,f_4,f_5]$',FontSize=font,Interpreter='latex')
ylabel('$\sigma_{OE}^2$,$\;$ $\sigma_{BJ}^2$',FontSize=font+4,Interpreter='latex')
grid on
% legend('Output error variance','Box Jenkins variance')

%% Bode plot first 20 iniitlization mc runs