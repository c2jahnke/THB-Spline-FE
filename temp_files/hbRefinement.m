%close all, clear all, clc
a = 0;
b = 10;
p = 1;
N = 5;
resol = 0.01;
lvl = 2;
obj = hbSplBasML(a,b,p,N,resol,lvl);

f = @(x) sin(x); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); % not needed
Points(:,1) = linspace(0,10, size(Points,1));
Points(:,2) = f(2*pi*Points(:,1));


cBas = obj.levelBas{1};
fBas = obj.levelBas{2};
refArea = [4 10];
%% spline curve
% figure(1)
% ImpointCurvePlot(cBas.n,cBas.sP,cBas.p,cBas.knotVector,cBas.plotVector,Points)
% movegui(1,'west')
% 
 [U, Ubar, Points, Qw] = HbRefinement1D(cBas,fBas,refArea,Points);
% 
% nOF = length(cBas.activeIndex) + length(fBas.activeIndex);
% 
% figure(3)
% 
% ImpointCurvePlot(nOF-1,cBas.sP,cBas.p,Ubar,cBas.plotVector,Qw)
% movegui(3,'east')

     cBas.plotBasisStruct(cBas.generBasisRed);
     hold on;
     fBas.plotBasisStruct(fBas.generBasisRed);

    
    
% [THB0, THB1, trunq, q, trP] = ThbRefinement1D(cBas,fBas,refArea,Qw); 
%   figure %
% for ll = 1:cBas.n
%         plot(fBas.plotVector,THB0(:,ll)); hold all
% end
% fBas.plotBasisStruct(THB1) % inner basis

%% Multiple refinement
% 
% cBas =  obj.levelBas{2};
% fBas =  obj.levelBas{3};
% refArea = [6 8];
% [U, Ubar, Points, Qw] = HbRefinement1D(cBas,fBas,refArea);
% [THB1_, THB2, trunq2, q2] = ThbRefinement1D(cBas,fBas,refArea);
% 
% for ll = 5+cBas.p:cBas.n
%         plot(fBas.plotVector,THB1_(:,ll));
%         hold all
% end
% 
% fBas.plotBasisStruct(THB2) % inner basis
