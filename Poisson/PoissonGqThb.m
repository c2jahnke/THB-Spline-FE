clc, clear all, close all
a = 0;
b = 9;
p = 2;
N = 9;
resol = 0.002;
lvl = 2;
obj = thbSplBasML(a,b,p,N,resol,lvl);
f = @(x) pi.^2/4*sin(pi/2*x);%sin(x); 
dBC = 0;% Dichlet BC
    % change generation of Points!
%Points = zeros(obj.levelBas{1}.n,2); % needed, just for Refinement


refArea = [6 8];
obj.ThbRefinement1DML(1,refArea,f);
refArea2 = [2 4];
obj.ThbRefinement1DML(1,refArea2,f);

if(isempty(refArea))
    nOE = N;
    nOF = obj.levelBas{1}.n;
else

nOF = length(obj.levelBas{1}.activeIndex) +length(obj.levelBas{2}.activeIndex); % number o functions
nOE =nOF - p; % number of elements
end
Stiffn = zeros(nOF); %basis.n
ngp = max(obj.levelBas{1}.p+1,sqrt(obj.levelBas{1}.p^2 -2*obj.levelBas{1}.p+1));
% first part to integrate right hand side sufficiently accurate, second
% part to generate stiffness matrix exactly
rhs = zeros(nOF,1); % number of basis functions!
allKnots = unique(obj.getAllKnots);

for el = 1 : nOE % loop over elements
    [s,w]=lgwt(ngp,allKnots(el),allKnots(el+1));
    bVal = []; % basis evaluation
    gradVal = []; % derivative evaluation
    for j = length(s):-1:1 
        temp = obj.evalDersBasis(s(j)); % method of class to evaluate derivatives!
        [lvl, Ind] = obj.getActiveFctIndU(s(j));
        lIndex = [lvl ; Ind];
        
        bVal(j,:) = temp(1,:);
        gradVal(j,:) = temp(2,:); 
    end
    elRhs = zeros(size(bVal,2),1); 
    elStiff = zeros(size(bVal,2));
    elSInd = cell(size(bVal,2));
    elRInd =cell(size(bVal,2),1);
    for ii0 = 1 : size(bVal,2) % cBas.p+1
        elRhs(ii0) = sum(w.*f(s).*bVal(:,ii0));
        elRInd{ii0} = lIndex(:,ii0);
        elStiff(ii0,ii0) = sum(w.*gradVal(:,ii0).^2);
        elSInd{ii0,ii0} = [lIndex(:,ii0) lIndex(:,ii0)];
        for jj = ii0+1 : size(bVal,2) 
            elStiff(ii0,jj) = sum(w.*gradVal(:,ii0).*gradVal(:,jj));
            elSInd{ii0,jj} = [lIndex(:,ii0) lIndex(:,jj)];
            elStiff(jj,ii0) = elStiff(ii0,jj);
            elSInd{jj,ii0} = elSInd{ii0,jj};
        end      
    end
    % generation of element stiffness matrix and of index matrix done!
    

    for l = 1 : obj.level
        for k = 0:length(obj.levelBas{l}.activeIndex)-1
            if(l >1)
            ind_1 = length(obj.levelBas{l-1}.activeIndex) +k+1;
            else
                ind_1 = k+1; % index for basis functions
            end
            for ii1 = 1 : size(bVal,2)
                if(elRInd{ii1} == [l;obj.levelBas{l}.activeIndex(k+1)])
                    
                   rhs(ind_1) = rhs(ind_1) + elRhs(ii1);
                   iLvl(ind_1) = l;
                   iBasisFctInd(ind_1) = obj.levelBas{l}.activeIndex(k+1);
                end
             end
            
            for ll = l : obj.level
                for kk = 0:length(obj.levelBas{ll}.activeIndex)-1 
                    if(ll >1)
                    ind_2 = length(obj.levelBas{ll-1}.activeIndex) +kk+1;
                    else
                    ind_2 = kk+1;
                    end
                    for  ii = 1 : size(bVal,2)
                        for jj = ii : size(bVal,2)
                        if(elSInd{ii,jj} == [l ll; obj.levelBas{l}.activeIndex(k+1)...
                                obj.levelBas{ll}.activeIndex(kk+1)])
               
                            Stiffn(ind_1,ind_2) = Stiffn(ind_1,ind_2) + elStiff(ii,jj);
                            Stiffn(ind_2,ind_1) = Stiffn(ind_1,ind_2);
                        end
                        end
                    end
                end 
            end
        end
    end
    
end
% reorder stiffness matrix and rhs to enforce boundary conditions
% only if refArea(1) == b

%% BC with Lagrange multipliers
A = zeros(nOE+obj.levelBas{1}.p +1);
b = zeros(nOE+obj.levelBas{1}.p +1,1);
A(1,2) = 1;
A(2,1) = 1;
A(2:end,2:end) = Stiffn;
b(1,1) = dBC; % Dichlet BC
b(2:end) = rhs;

u = A\b; 
%u_ord = reorderU(obj,u,nOF,iLvl,iBasisFctInd); % reorder u
%y =  u_ord(2:end);
%% reordering is not necessary!
y = u(2:end);
uh =  generSolThb(obj,y);
hold on;
g = @(x) sin(pi/2*x);
fplot(g,[obj.levelBas{1}.a obj.levelBas{1}.b],'b')

errorEst(uh,g,obj);
%%
figure
obj.levelBas{1}.plotBasisStruct(obj.levelBas{1}.generBasisRed(obj.levelBas{2}))
hold on
obj.levelBas{2}.plotBasisStruct(obj.levelBas{2}.generBasisRed(obj.levelBas{2}))
