% Project part 4 system identification Michiel Wind, Jelle Cruijsen
clear all, close all, clc

% Determine frequency response butterworth filter
z = tf('z',1);
F = (0.505+1.01*z^(-1)+0.505*z^(-2))/(1+0.7478*z^(-1)+0.2722*z^(-2));
freq = (0:0.001:pi);
freq_query = getGainCrossover(F,db2mag(-3)); % Crossover frequency at -3dB
bandpass_freq = freq_query/pi;
% 
N = 3000;
p = 0.6;
Range = [-3,3];
Band = [0 bandpass_freq];
index_min = fix(N/2 - N/2*bandpass_freq);
index_plus = fix(N/2 + N/2*bandpass_freq);
% r = prbs(N,max(Range),p);
r = idinput(N,'prbs',Band,Range);
stem(abs(fftshift(fft(r))))
hold on
xline(index_min,'LineWidth',2)
xline(index_plus,'LineWidth',2)
xlabel('Frequency index [-]')
ylabel('Amplitude [-]')
[u,y] = assignment_sys_20(r);
data = iddata(y,u);
%% parametric ID
na = 5;
nb = 5;
nk = 1;
% Most likely a 5th order OE model...
sys1 = arx(data,[na nb nk],arxOptions('Focus','prediction'));
sys2 = oe(data,[na nb nk],arxOptions('Focus','prediction'));
figure
resid([y u],sys1)
figure
resid([y u],sys2)
figure
bode(sys1)
hold on
bode(sys2)
legend('ARX','OE')

% Which model order to take:
% V = arxstruc(ze,zv,NN) !!
% Validation
%- check with respect to ETFE, but be careful since ETFE can also be wrong
%- Check pole zero map for maybe pole zero cancelations which will imply a
%too high model order.
%-residual test data

%% Step response of the system to find the delay
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
% So the delay is nk = 1;

%% Cross validation test
close all
na = 15;
nb = 15;
nk = 1;
% Define the maximum model order to consider
max_order = na+nb+nk;

M_arx = cell(na,nb);
err_arx = cell(na,nb);
V_arx = cell(na,nb);


% Cross validation tests ARX

for i = 1:na
    for j = 1:nb
        M_arx{i,j} = arx(data,[i j nk]);
        err_arx{i,j} = pe(M_arx{i,j},data);
        V_arx{i,j} = sum(err_arx{i,j}.y.^2)/(length(err_arx{i,j}.y));
    end
end

% Cross validation
[u,y] = assignment_sys_20(r);
crossdata = iddata(y,u);

err_arx_cross = cell(na,nb);
V_arx_cross = cell(na,nb);

for i = 1:na
    for j = 1:nb
        err_arx_cross{i,j} = pe(M_arx{i,j},crossdata);
        V_arx_cross{i,j} = sum(err_arx_cross{i,j}.y.^2)/(length(err_arx_cross{i,j}.y));
    end
end


% OE model
opt = oeOptions('WeightingFilter',[0 0.5*pi],'Focus','prediction');

nb = 15;
nf = 15;
nk = 1;

err_oe = cell(nf,nb);
V_oe = cell(nf,nb);
M_oe = cell(nf,nb);

for i = 1:nf
    for j = 1:nb
        M_oe{i,j} = oe(data,[j i nk],opt);
        err_oe{i,j} = pe(M_oe{i,j},data);
        V_oe{i,j} = sum(err_oe{i,j}.y.^2)/(length(err_oe{i,j}.y));
    end
end

% Cross validation OE
err_oe_cross = cell(nf,nb);
V_oe_cross = cell(nf,nb);

for i = 1:nf
    for j = 1:nb
        err_oe_cross{i,j} = pe(M_oe{i,j},crossdata);
        V_oe_cross{i,j} = sum(err_oe_cross{i,j}.y.^2)/(length(err_oe_cross{i,j}.y));
    end
end

%% Plotting
figure()
for i = 1:na
    for j = 1:nb
    scatter(i+j,V_oe_cross{i,j})
    hold on
    end
end
xlim([0 i+j])
V_mat = cell2mat(V_oe_cross);
[maxValue, maxIndex] = max(V_mat(:));
[row, col] = ind2sub(size(V_mat), maxIndex);
figure
bode(M_oe{row,col})
[minValue, minIndex] = min(V_mat(:));
[row, col] = ind2sub(size(V_mat), minIndex);
figure
bode(M_oe{row,col})
% Spa function??