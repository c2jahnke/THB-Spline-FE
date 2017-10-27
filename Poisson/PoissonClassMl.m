clc, clear all, close all
% f = parameters.f;
% g = parameters.g;
a = -1;
b = 1;
p = 10;
N = 80;
resol = 0.01;
lvl = 5;% not more than 4 levels work!
obj = thbSplBasML(a,b,p,N,resol,lvl);
fBas = obj.levelBas{2};
cBas = obj.levelBas{1};
sigma = 0.05;
mu = 0;

f = @(x) - 1/sqrt(2*pi*sigma^2) * ((x-mu).*(x-mu)/sigma.^4 - 1/sigma^2).*exp(-1/2*(x-mu).*(x-mu)/sigma.^2);%pi.^2/4*sin(pi/2*x);%source function 
g = @(x) 1/sqrt(2*pi*sigma^2).*exp(-1/2*(x-mu).*(x-mu)/sigma^2);%sin(pi/2*x);

dBC = boundCond('Dirichlet','Dirichlet',0,0);% Dichlet BC
% refArea = [];
% refArea = [-6 6]; % refArea = [4/3 8/3];
% obj.ThbRefinement1DML(1,refArea);
% refArea2 = [-4 4];%[3*pi/8 7*pi/8];
% obj.ThbRefinement1DML(2,refArea2);
% refArea3 = [-2 2];%[4*pi/8 6*pi/8 ];
% obj.ThbRefinement1DML(3,refArea3);
% refArea4 = [-1 1];%[4*pi/8 6*pi/8 ];
% obj.ThbRefinement1DML(4,refArea4);
ps = PoissSolv(obj,f);

obj.plotBas
[Stiffn, rhs, iLvl,iBasisFctInd] = ps.assembleMl();

y = ps.solveSyst(Stiffn,rhs,dBC);

fplot(g,[obj.levelBas{1}.a obj.levelBas{1}.b],'b')
hold on;
uh =  ps.generSolThb(y);%
hold on;
plot(obj.getAllKnots,0.01,'k*', 'markers',4)
plot(obj.levelBas{1}.plotVector,uh,'r')
errorEst(uh,g,obj);
[H1err, L2err, Enerr] = ps.errCalc(g,y)