clear all, close all, clc
load('batch.mat')
data = iddata(y,u);
na=5;
nb=5;
nk=1;
OE = oe(data,[na nb nk],arxOptions('Focus','prediction'));
figure
resid([y u],OE)
nb = 5;
nc = 5;
nd = 5;
nf = 5;
nk = 1;

BJ = bj(data,[nb nc nd nf nk],arxOptions('Focus','prediction'));
figure
resid([y u],BJ)
%%

[A,B,C,D,F] = polydata(BJ); % Retrieve polynomial coeff.
H = tf(C,D,1); G = tf(B,F,1); % Estimated TFs
[Phiu,w] = pwelch(u); % Estimate input power spectrum
[Gmag,~] = bode(G,w);
[Hmag,~] = bode(H,w); % TF magnitudes
sigma_ehat = BJ.Report.Fit.MSE; % Estimate variance of e(t)
Phiv = squeeze(Hmag).^2*sigma_ehat; % Estimate noise power spectrum
Phiy = squeeze(Gmag).^2.*Phiu+Phiv; % Estimate output power spectrum
figure; % Plot power spectra (ratios)
loglog(w,Phiv,'linewidth',2); hold on;
loglog(w,Phiv./Phiu,'linewidth',2);
loglog(w,Phiv./Phiy,'linewidth',2);
grid on