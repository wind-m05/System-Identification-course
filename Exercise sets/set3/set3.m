close all, clc, clear all
load('dataG0arx.mat')
z = tf('z',1);
b1 = 0.10276;
b2 = 0.18123;
f1 = 1.99185;
f2 = 2.20265;
f3 = -1.84083;
f4 = 0.89413;
G0 = (z^(-3)*b1+z^(-4)*b2)/(1-f1*z^(-1)+f2*z^(-2)+f3*z^(-3)+f4*z^(-4));

%% ARX
G_ARX = arx([y u],[4 2 3]);
resid([y u],G_ARX)
%% OE
G_OE = oe([y u],[2 4 3]);
resid([y u],G_OE)

%% BODE plots
bode(G0)
hold on
bode(G_ARX)
bode(G_OE)

%% Question 4
clear all, close all, clc
N = 500; % Number of data points
B = [1 -.5 .2 .8 .1]; % Numerator coefficients of G 0
A = [1 -1.7 1.6 -.8 .25]; % Denominator coefficients of G 0
G_0 = tf(B,A,-1,'Variable','z^-1'); % TF of the true system G 0
H_0 = tf(1,A,-1,'Variable','z^-1'); % TF of the true noise model H 0
u = randn(N,1); % Gaussian white noise w/ unit variance
e = randn(N,1); % Gaussian white noise w/ unit variance
y = lsim(G_0,u)+lsim(H_0,e); % Simulated output

G_ARX = arx([y u], [4,5,0], arxOptions('Focus','prediction')); % Alter the weighting filter such that the model is great predictor.
resid(G_ARX,[y u])
pzmap(G_ARX)
num_poles = length(pole(G_ARX))
num_zeros = length(zero(G_ARX))

bode(G_0)
hold on
bode(G_ARX)

%% FIR model
clear all,  close all, clc
N = 500; % Number of data points
B = [1 -.5 .2 .8 .1]; % Numerator coefficients of G 0
A = [1 -1.7 1.6 -.8 .25]; % Denominator coefficients of G 0
G_0 = tf(B,A,-1,'Variable','z^-1'); % TF of the true system G 0
H_0 = tf(1,A,-1,'Variable','z^-1'); % TF of the true noise model H 0
u = randn(N,1); % Gaussian white noise w/ unit variance
e = randn(N,1); % Gaussian white noise w/ unit variance
y = lsim(G_0,u)+lsim(H_0,e); % Simulated output
G_FIR = arx([y u], [0,12,0], arxOptions('Focus','prediction'));
bode(G_0)
hold on
bode(G_FIR)
legend('G0','FIR')
%% 3th order ARX, OE
rng(42)
clear all,  close all, clc
N = 500; % Number of data points
B = [1 -.5 .2 .8 .1]; % Numerator coefficients of G 0
A = [1 -1.7 1.6 -.8 .25]; % Denominator coefficients of G 0
G_0 = tf(B,A,-1,'Variable','z^-1'); % TF of the true system G 0
H_0 = tf(1,A,-1,'Variable','z^-1'); % TF of the true noise model H 0
u = randn(N,1); % Gaussian white noise w/ unit variance
e = randn(N,1); % Gaussian white noise w/ unit variance
y = lsim(G_0,u)+lsim(H_0,e); % Simulated output
G_ARX3 = arx([y u], [3,3,0], arxOptions('Focus','prediction'));
G_OE3 =   oe([y u], [3,3,0], arxOptions('Focus','prediction'));

bode(G_0)
hold on
bode(G_ARX3)
bode(G_OE3)
legend('G0','G arx','G OE')

num_poles = length(pole(G_OE3))
num_zeros = length(zero(G_OE3))
%% Prefiltering for ARX fitting
clear all,  close all, clc
rng(41)
N = 500; % Number of data points
B = [1 -.5 .2 .8 .1]; % Numerator coefficients of G 0
A = [1 -1.7 1.6 -.8 .25]; % Denominator coefficients of G 0
G_0 = tf(B,A,-1,'Variable','z^-1'); % TF of the true system G 0
H_0 = tf(1,A,-1,'Variable','z^-1'); % TF of the true noise model H 0
u = randn(N,1); % Gaussian white noise w/ unit variance
e = randn(N,1); % Gaussian white noise w/ unit variance
y = lsim(G_0,u)+lsim(H_0,e); % Simulated output
omega_p = 1;
fdata = idfilt([y u],[0,omega_p]);
G_ARX3 = arx([y u], [3,3,0], arxOptions('Focus','prediction'));
G_ARX3f = arx(fdata, [3,3,0], arxOptions('Focus','prediction'));
G_OE3 =   oe([y u], [3,3,0]);

bode(G_0)
hold on
bode(G_ARX3)
bode(G_ARX3f)
bode(G_OE3)
legend('Original','ARX','ARX prefiltered','OE')