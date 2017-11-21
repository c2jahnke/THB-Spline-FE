% scratch file for HB-splines
close all, clear all, clc
a = 0;
b = 10;
p = 1;
N = 5;
resol = 0.1;
lvl = 2;
obj = hbSplBasML(a,b,p,N,resol,lvl);

f = @(x) sin(x); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); 

cBas = obj.levelBas{1};
fBas = obj.levelBas{2};
refArea = [4 8];
%% spline curve
% figure(1)
% ImpointCurvePlot(cBas.n,cBas.sP,cBas.p,cBas.knotVector,cBas.plotVector,Points)
% movegui(1,'west')
% 
 [~, ~, ~, Qw] = HbRefinement1D(cBas,fBas,refArea,Points);


     cBas.plotBasisStruct(cBas.generBasisRed);
     hold on;
     fBas.plotBasisStruct(fBas.generBasisRed);

    %THB0 coarse functions truncated
    
[THB0, THB1, trunq, q, trP] = ThbRefinement1D(cBas,fBas,refArea,Qw); 
  figure %
for ll = cBas.activeIndex
        plot(fBas.plotVector,THB0(:,ll+1)); hold all
end
fBas.plotBasisStruct(fBas.generBasisRed);
hold off;
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
