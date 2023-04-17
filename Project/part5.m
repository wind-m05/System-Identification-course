% Project part 5 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc
% font = 12;
% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
N = 1000;
Range = [-3,3];
Band = [0 bandpass_freq];

%% Monte carlo runs 
% nb = 5;
% nf = 5;
% nk = 1;
% sims = 100;
% model_runs = cell(sims,1);
% opt = oeOptions('focus','prediction');
% for i = 1:sims
% r = idinput(N,'rbs',Band,Range);
% [u,y] = assignment_sys_20(r); % Input every loop or not
% data = iddata(y,u);
% sys = oe(data,[nb nf nk],opt); % y,u instead of data
% model_runs{i} = sys;
% end
% save('OE_montecarlo','model_runs')

%% pzmap/bode plot for visualization
% for i = 1:length(model_runs)
% pzmap(model_runs{i})
% hold on
% end

%% Theoretical covariance (We choose model_runs{3}, because its most likely a global minimizer)
load('OE_montecarlo')
Prho = getcov(model_runs{3},'value');
imagesc(Prho);
colorbar;
title('Covariance Matrix');
xlabel('Variable');
ylabel('Variable');

% Monte Carlo variance and mean
aux_F = [];
aux_B = [];
for i = 1:length(model_runs)
    aux_F = [aux_F;model_runs{i}.F];
    aux_B = [aux_B;model_runs{i}.B];
end
mean_F = mean(aux_F);
median_F = median(aux_F);
var_F = var(aux_F);
mean_B = mean(aux_B);
var_B = var(aux_B);
median_B = median(aux_B);
save('median_B','median_B')
save('median_F','median_F')

% Compare the Monte Carlo variance with the theoretical variance
mean_vec = [mean_B(2:end) mean_F(2:end)];
var_vec = [var_B(2:end) var_F(2:end)];
theoretical_var = diag(Prho)';
bar_aux = [];
for i = 1:length(var_vec)
    bar_aux = [bar_aux; [var_vec(i),theoretical_var(i)]];
end
figure
bar(bar_aux)
legend('Monte Carlo variance','Theoretical variance')
Compare1 = [[var_B(2:end) var_F(2:end)];theoretical_var];

%% Monte Carlo with initialization of the points
% M_init = idpoly([],median_B,[],[],median_F);
% sims = 100;
% model_runs = cell(sims,1);
% opt.Focus = 'prediction';
% for i = 1:sims
% r = idinput(N,'rbs',Band,Range);
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u);
% sys = oe(data,M_init,opt); % OE
% model_runs{i} = sys;
% end
% save('OE_montecarlo_initialized','model_runs')

%% Variance analysis
%Monte carlo variance and mean

model_runs_init = importdata('OE_montecarlo_initialized.mat');

% Monte Carlo variance and mean
aux_F_init = [];
aux_B_init = [];
for i = 1:length(model_runs_init)
    aux_F_init = [aux_F_init;model_runs_init{i}.F];
    aux_B_init = [aux_B_init;model_runs_init{i}.B];
end

var_F_init = var(aux_F_init,0,1);
var_B_init = var(aux_B_init,0,1);

% Compare the Monte Carlo variance with the theoretical variance
var_vec_init = [var_B_init(2:end) var_F_init(2:end)];
theoretical_var = diag(Prho)';
bar_aux_init = [];
for i = 1:length(var_vec_init)
    bar_aux_init = [bar_aux_init; [var_vec_init(i),theoretical_var(i)]];
end
figure
bar(bar_aux_init)
legend('Monte Carlo variance','Theoretical variance')
Compare2 = [[var_B_init(2:end) var_F_init(2:end)];theoretical_var];

%% See difference bode plots
figure
for i = 1:length(model_runs)
bode(model_runs{i})
hold on
end
