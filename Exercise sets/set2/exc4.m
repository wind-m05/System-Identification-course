close all, clc
load('dataG0oe.mat')
z = tf('z',1);
b1 = 0.10276;
b2 = 0.18123;
f1 = 1.99185;
f2 = 2.20265;
f3 = -1.84083;
f4 = 0.89413;
G0 = (z^(-3)*b1+z^(-4)*b2)/(1-f1*z^(-1)+f2*z^(-2)+f3*z^(-3)+f4*z^(-4));

bode(G0)
hold on
bode(oe243)
bode(arx423)
legend('G0','OE','ARX')