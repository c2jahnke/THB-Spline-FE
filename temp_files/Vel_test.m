clc, clear all

N = 1000;
bas = bSplBas(0,10,3,N,0.001);


tic
D = bas.generBasis2;
disp('Time: ');
toc


tic
C = bas.generBasis;
disp('Time: ');
toc