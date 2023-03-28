%% Box jenkins
% 
nb = 5;
nc = 5;
nd = 7;
nf = 7;
nk = 1;

BJ = bj(data,[nb nc nd nf nk]);
figure
resid([y u],BJ)