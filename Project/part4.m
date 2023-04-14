% Project part 4 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc
font = 12;
% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
% 
N = 3000;
p = 0.6; % probability of switching in prbs.
Range = [-3,3];
Band = [0 bandpass_freq];
index_min = fix(N/2 - N/2*bandpass_freq);
index_plus = fix(N/2 + N/2*bandpass_freq);
r = prbs(N,max(Range),p);
% r = idinput(N,'prbs',Band,Range);
stem(abs(fftshift(fft(r))))
hold on
xline(index_min,'LineWidth',2)
xline(index_plus,'LineWidth',2)
xlabel('Frequency index [-]')
ylabel('Amplitude [-]')
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
%% parametric ID OE
nb = 5;
nf = 5;
nk = 1;
sys1 = oe(data,[nb nf nk],oeOptions('Focus','prediction')); % OE
figure
resid([y u],sys1);
%% parametric ID ARX
na = 5;
nb = 5;
nk = 1;
sys2 = arx(data,[na nb nk],arxOptions('Focus','prediction')); 
figure
resid([y u],sys2)
figure
bode(sys2)
%% parametric ID ARMAX
na = 5;
nb = 5;
nc = 1;
nk = 1;
sys3 = armax(data,[na nb nc nk],arxOptions('Focus','prediction')); 
figure
resid([y u],sys3)
figure
bode(sys3)

%% Find the delay with the use of the crosscorrelation
% Make white noise signal u
N = 1000000;
r = randn(1,N)*1;
[u,y] = assignment_sys_20(r);
xgrid = -(N-1):(N-1);

% Proof that u is white?
utest = xcorr(u);
figure
stem(xgrid,utest,'Linewidth',2)
xlim([-10 10])
title('Autocorrelation of the input signal u(t)')
xlabel('Lags',FontSize=font,Interpreter='latex')
ylabel('$Amplitude [-]$',FontSize=font,Interpreter='latex')
grid on

r = randn(1,N)*1;
[u,y] = assignment_sys_20(r);
crco = xcorr(y,u);
crco = crco./max(crco);
figure
stem(xgrid,crco,'Linewidth',2)
xlim([-10 10])
title('Normalized Cross correlation between output and input')
xlabel('Lags',FontSize=font,Interpreter='latex')
ylabel('$Amplitude [-]$',FontSize=font,Interpreter='latex')
grid on

%% Find the delay with step responses
close all
N = 3000;
r = zeros(N,1);
trials = 10000;
% z = zeros(N/2,1);
% o = ones(N/2,1);
r(1:N/2,1) = 0;
r(N/2:end,1) = 1;
y_sum = zeros(N,1);
for i = 1:trials
[u,y] = assignment_sys_20(r);
y_sum = y_sum + y;
end
y_avg = y_sum/trials;
figure
stem(y_avg)
% So now we know the delay is nk = 1; for sure

%% second validation test
yh = compare(data,sys1,1000);





%% Everything below this is backup, because the approach is wrong: for 
% an ARX structure it is not possible to consistently estimate G0 alone, 
% because G and H are not parameterized independently

% %% ARX model build
% clc, close all
% band = [0.4 bandpass_freq*pi];
% 
% [u,y] = assignment_sys_20(r);
% data = iddata(y,u);
% 
% na = 15;
% nb = 15;
% nk = 1;
% % Define the maximum model order to consider
% max_order = na+nb+nk;
% 
% M_arx = cell(na,nb); % ARX models
% err_arx = cell(na,nb); % ARX model errors 
% V_arx = cell(na,nb); % ARX models value functions
% opt = arxOptions('WeightingFilter',band,'Focus','prediction'); % Adjust
% for i = 1:na
%     for j = 1:nb
%         M_arx{i,j} = arx(data,[i j nk],opt);
%         err_arx{i,j} = pe(M_arx{i,j},data);
%         V_arx{i,j} = err_arx{i,j}.y'*err_arx{i,j}.y/(length(err_arx{i,j}.y));
%     end
% end
% 
% % Cross validation ARX models
% [u,y] = assignment_sys_20(r);
% crossdata = iddata(y,u);
% 
% err_arx_cross = cell(na,nb);
% V_arx_cross = cell(na,nb);
% 
% for i = 1:na
%     for j = 1:nb
%         err_arx_cross{i,j} = pe(M_arx{i,j},crossdata);
%         V_arx_cross{i,j} = sum(err_arx_cross{i,j}.y.^2)/(length(err_arx_cross{i,j}.y));
%     end
% end
% 
% 
% %% OE models build
% opt = oeOptions('WeightingFilter',band,'Focus','prediction'); % Adjust
% 
% nb = 10;
% nf = 10;
% nk = 1;
% 
% err_oe = cell(nf,nb);
% V_oe = cell(nf,nb);
% M_oe = cell(nf,nb);
% 
% for i = 1:nf
%     for j = 1:nb
%         M_oe{i,j} = oe(data,[j i nk],opt);
%         err_oe{i,j} = pe(M_oe{i,j},data);
%         V_oe{i,j} = (err_oe{i,j}.y'*err_oe{i,j}.y)/(length(err_oe{i,j}.y));
%     end
% end
% 
% % Cross validation OE models
% err_oe_cross = cell(nf,nb);
% V_oe_cross = cell(nf,nb);
% 
% for i = 1:nf
%     for j = 1:nb
%         err_oe_cross{i,j} = pe(M_oe{i,j},crossdata);
%         V_oe_cross{i,j} = (err_oe_cross{i,j}.y'*err_oe_cross{i,j}.y)/(length(err_oe_cross{i,j}.y));
%     end
% end
% 
% %% Plotting OE models
% figure()
% for i = 1:nf
%     for j = 1:nb
%     scatter(i+j,V_oe_cross{i,j})
%     hold on
%     end
% end
% xlim([0 i+j])
% ylim([min(min(cell2mat(V_oe_cross))) min(min(cell2mat(V_oe_cross)))+2])
% 
% %% Plotting ARX models
% figure()
% for i = 1:na
%     for j = 1:nb
%     scatter(i+j,V_arx_cross{i,j})
%     hold on
%     end
% end
% xlim([0 i+j])
% ylim([min(min(cell2mat(V_arx_cross))) min(min(cell2mat(V_arx_cross)))+1])
% 
% %% Plot bode of best and worst cost function values
% V_mat = cell2mat(V_oe_cross);
% [maxValue, maxIndex] = max(V_mat(:));
% [row, col] = ind2sub(size(V_mat), maxIndex);
% figure
% bode(M_oe{row,col})
% [minValue, minIndex] = min(V_mat(:));
% [row, col] = ind2sub(size(V_mat), minIndex);
% figure
% bode(M_oe{row,col})
% % Spa function for noise spectrum??
% 
% %% matlab functions test (from the lecture notes and then matlab function in the chapter model validations)
% [u,y] = assignment_sys_20(r);
% z = iddata(y,u);
% [u,y] = assignment_sys_20(r);
% nn = [30 30 5];
% zv = iddata(y,u);
% v = arxstruc(z,zv,nn);
% order = selstruc(v,0);
% 
% %% Aikake's information criterion (AIC) for OE models
% V_mat = cell2mat(V_oe_cross);
% V_AIC = zeros(size(V_mat));
% for i = 1:length(V_mat)
%     for j = 1:length(V_mat)
%        V_AIC(i,j) = 1/2*log(V_mat(i,j)) + (i+j+nk)/length(r);
%     end
% end
% [minvalue, minindex] = min(V_AIC(:));
% [row, col] = ind2sub(size(V_AIC), minindex);
% nn = [row col nk]; % This gives the right values... [5,5,1]
% 
% %% Aikake's information criterion (AIC) for ARX models
% V_mat = cell2mat(V_arx_cross);
% V_AIC = zeros(size(V_mat));
% for i = 1:length(V_mat)
%     for j = 1:length(V_mat)
%        V_AIC(i,j) = 1/2*log(V_mat(i,j)) + (i+j+nk)/length(r);
%     end
% end
% [minvalue, minindex] = min(V_AIC(:));
% [row, col] = ind2sub(size(V_AIC), minindex);
% nn = [row col nk];
% 
% %% Plot only the minimum of the scatter plot per order
% V_mat = cell2mat(V_oe_cross);
% figure()
% A = cell(nf+nb,1);
% A{1,1} = nan; % No possible permutation to get 1.
% for k = 2:nf+nb % Number of parameters
%     for i = 1:nf
%         for j = 1:nb
%             if i+j == k
%             A{k} = [A{k},V_mat(i,j)]; 
%         %         scatter(i+j,V_oe_cross{i,j})
%         %         hold on
%             end
%         end
%     end
%     A{k} = min(A{k});
% end
% A = cell2mat(A);
% scatter(1:length(A),A,'*')
