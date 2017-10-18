clc, clear all, close all
% f = parameters.f;
% g = parameters.g;
a = -1;
b = 1;
p = 1;
N = 10;
resol = 0.001;
lvl = 5;% not more than 3 levels work!
obj = thbSplBasML(a,b,p,N,resol,lvl);
sigma = 0.1;
mu = 0;

f = @(x) - 1/sqrt(2*pi*sigma^2) * ((x-mu).*(x-mu)/sigma.^4 - 1/sigma^2).*exp(-1/2*(x-mu).*(x-mu)/sigma.^2);%pi.^2/4*sin(pi/2*x);%source function 
g = @(x) 1/sqrt(2*pi*sigma^2).*exp(-1/2*(x-mu).*(x-mu)/sigma^2);%sin(pi/2*x);

dBC = boundCond('Dirichlet','Dirichlet',0,0);% Dichlet BC
refArea = [];
refArea = [-0.4 0.4]; % refArea = [4/3 8/3];
obj.ThbRefinement1DML(1,refArea,f);
% refArea2 = [-0.4 0.2];%[3*pi/8 7*pi/8];
% obj.ThbRefinement1DML(2,refArea2,f);
% refArea3 = [-0.4 0.4];%[4*pi/8 6*pi/8 ];
% obj.ThbRefinement1DML(3,refArea3,f);
% refArea4 = [-0.2 0.2];%[4*pi/8 6*pi/8 ];
% obj.ThbRefinement1DML(4,refArea3,f);

obj.plotBas
[Stiffn, rhs, iLvl,iBasisFctInd] = assemblePoissThbMl(obj,refArea,f);

y = solveSyst(obj,Stiffn,rhs,iLvl,iBasisFctInd,dBC);

fplot(g,[obj.levelBas{1}.a obj.levelBas{1}.b],'b')
hold on;
uh =  generSolThb(obj,y);
hold on;
plot(obj.getAllKnots,0.01,'k*', 'markers',4)
errorEst(uh,g,obj);
[H1err, L2err, Enerr] = errCalc(obj,g,y,iBasisFctInd,iLvl)