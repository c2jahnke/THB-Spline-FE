close all, clear all, clc
a = 0;
b = 10;
p = 2;
N = 10;
resol = 0.001;
lvl = 2;
obj = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x) sin(x); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); 

cBas = obj.levelBas{1};
fBas = obj.levelBas{2};
refArea = [2 8]; 
%THB0 coarse functions truncated
tic  
[THB0, THB1, trunq, q, trP] = obj.ThbRefinement1DML(1,refArea,f); 
disp('Old refinement: ')
toc
%   figure %
%   hold all
% for ll = cBas.activeIndex
%         plot(fBas.plotVector,THB0(:,ll+1)); 
% end
% hold off;
tic
D = cBas.generBasisRed(fBas);
disp('Time: ')
toc


cBas.plotBasisStruct(D);
hold all
fBas.plotBasisStruct(fBas.generBasisRed(fBas));
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
