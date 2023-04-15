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

%% Monte carlo runs 
% nb = 5;
% nf = 5;
% nk = 1;
% sims = 100;
% model_runs = cell(sims,1);
% for i = 1:sims
% r = idinput(N,'rbs',Band,Range);
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u);
% sys = oe(data,[nb nf nk],oeOptions('Focus','prediction')); % OE
% model_runs{i} = sys;
% end
% save('OE_montecarlo','model_runs')

%% pzmap/bode plot for visualization
% for i = 1:length(model_runs)
% pzmap(model_runs{i})
% hold on
% end

%% Theoretical covariance (We choose model_runs{3}, because its most likely a global minimizer)
% load('OE_montecarlo_original')
model_runs{3}.B = [model_runs{3}.B];
model_runs{3}.F = [model_runs{3}.F];
Prho = getcov(model_runs{3});
imagesc(Prho);
colorbar;
title('Covariance Matrix');
xlabel('Variable');
ylabel('Variable');
% Prho normalized
Prho_norm = Prho./max(max(Prho));
% Test parameter vector
sys = model_runs{3};
vector = [model_runs{3}.B(2:end) model_runs{3}.F(2:end)];
% Monte carlo variance and mean
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

% Compare the Monte Carlo variance with the theoretical variance
mean_vec = [mean_B(2:end) mean_F(2:end)];
var_vec = [var_B(2:end) var_F(2:end)];
var_vec_norm = var_vec./max(var_vec);
theoretical_var = diag(Prho_norm)';
bar_aux = [];
for i = 1:length(var_vec)
    bar_aux = [bar_aux; [var_vec_norm(i),theoretical_var(i)]];
end
bar(bar_aux)

%% Monte Carlo with initialization of the points
% M_init = idpoly([],median_B,[],[],median_F);
% nb = 5;
% nf = 5;
% nk = 1;
% sims = 100;
% model_runs = cell(sims,1);
% for i = 1:sims
% r = idinput(N,'rbs',Band,Range);
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u);
% sys = oe(data,M_init,oeOptions('Focus','prediction')); % OE
% model_runs{i} = sys;
% end
% save('OE_montecarlo_initialized','model_runs')

%% Variance analysis
% Monte carlo variance and mean
% model_runs_init = importdata('OE_montecarlo_initialized_original.mat');
% aux_F_init = [];
% aux_B_init = [];
% for i = 1:length(model_runs_init)
%     aux_F_init = [aux_F;model_runs_init{i}.F];
%     aux_B_init = [aux_B;model_runs_init{i}.B];
% end
% mean_F_init = mean(aux_F_init);
% var_F_init = var(aux_F_init);
% mean_B_init = mean(aux_B_init);
% var_B_init = var(aux_B_init);
% 
% % Compare the Monte Carlo variance with the theoretical variance
% mean_vec_init = [mean_F_init(2:end) mean_B_init(2:end)];
% var_vec_init = [var_F_init(2:end) var_B_init(2:end)];
% var_vec_norm_init = var_vec_init./max(var_vec_init);
% bar_aux_init = [];
% for i = 1:length(var_vec)
%     bar_aux_init = [bar_aux_init; [var_vec_norm_init(i),theoretical_var(i)]];
% end
% figure
% bar(bar_aux_init)
% 
% % Check if there are still outliers after initialization
% for i = 1:length(model_runs_init)
% bode(model_runs_init{i})
% hold on
% end
