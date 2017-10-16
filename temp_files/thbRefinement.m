close all, clear all, clc
a = -1;
b = 1;
p = 3;
N = 4;
resol = 0.1;
lvl = 3;
obj = thbSplBasML(a,b,p,N,resol,lvl);

f = @(x) sin(x); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); 
refArea = [-0.5 0.5]; 

[THB0, THB1, trunq, q, trP] = obj.ThbRefinement1DML(1,refArea,f); 
refArea2 = [-0.25 0.25];
[THB1_, THB2, trunq2, q2, trP2] = obj.ThbRefinement1DML(2,refArea2,f);

% u = 0.55
% test = obj.evalBasis(u)
% [lvl, Ind] = obj.getActiveFctIndU(u)

obj2 = thbSplBasML(a,b,p,N,resol,2);
for k = obj.levelBas{1}.activeIndex
    figure
plot2DTHB(obj,1,k)
end
%plotOne2Dbasis(p,m,knotVector,iU,iV,plotVector,sP,resol)
a 
%THB0 coarse functions truncated
% tic  
% [THB0, THB1, trunq, q, trP] = obj.ThbRefinement1DML(1,refArea,f); 
% disp('Old refinement: ')
% toc
% 
% tic
% D = cBas.generBasisRed(fBas);
% disp('Class refinement: ')
% toc

% Multiple refinement
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
