function Stiffn = PoissonErrorFinding()
a = -1;
b = 1;
p = 1;
N = 2;
resol = 0.001;
lvl = 2;
obj = hbSplBasML(a,b,p,N,resol,lvl);
f = @(x) 3; % pi.^2/4*sin(pi/2*x);%3;%.^(1/31); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); % needed, just for Refinement


refArea = [-1 1];
[U, Ubar, ~, Qw] = HbRefinement1D(obj.levelBas{1},obj.levelBas{2},refArea,Points);

if(isempty(refArea))
    nOE = N;
    nOF = obj.levelBas{1}.n;
else
    %% not correct for p > 1
nOF = length(obj.levelBas{1}.activeIndex) +length(obj.levelBas{2}.activeIndex); % number o functions
nOE =nOF - p; % number of elements
end
Stiffn = zeros(length(obj.levelBas{1}.activeIndex) +length(obj.levelBas{2}.activeIndex)); %basis.n % length(cBas.activeIndex) +length(fBas.activeIndex) 
% Ax = b
% Lapl(u) = f
ngp = max(obj.levelBas{1}.p+1,sqrt(obj.levelBas{1}.p^2 -2*obj.levelBas{1}.p+1));
% first part to integrate right hand side sufficiently accurate, second
% part to generate stiffness matrix exactly
rhs = zeros(nOF,1); % basis.n number of basis functions!
allKnots = unique(obj.getAllKnots);

% number of elements: N+(refArea(2) - refArea(1))/cBas.knotspan
for el = 1 : nOE % loop over elements
    [s,w]=lgwt(ngp,allKnots(el),allKnots(el+1));
    bVal = []; % zeros(ngp, basis.p+1); % basis evaluation
    gradVal = []; % zeros(ngp,basis.p+1); % derivative evaluation
    for j = length(s):-1:1 %basis.evalDersBasis(s(j))% not yet fully tested!
        temp = obj.evalDersBasis(s(j)); %DersBasisFuns(k,s(j),basis.p,basis.knotVector,1) % replace by method of class!
        [lvl, Ind] = obj.getActiveFctIndU(s(j));
        lIndex = [lvl ; Ind];
        
        bVal(j,:) = temp(1,:);
        gradVal(j,:) = temp(2,:); % 
    end
    elRhs = zeros(size(bVal,2),1); % dynamic
    elStiff = zeros(size(bVal,2));
    elSInd = cell(size(bVal,2));
    elRInd =cell(size(bVal,2),1);
    for ii0 = 1 : size(bVal,2) % cBas.p+1
        elRhs(ii0) = sum(w.*f(s).*bVal(:,ii0));
        elRInd{ii0} = lIndex(:,ii0);
        lIndex(:,ii0);
        elStiff(ii0,ii0) = sum(w.*gradVal(:,ii0).^2);
        elSInd{ii0,ii0} = [lIndex(:,ii0) lIndex(:,ii0)];
        [lIndex(:,ii0) lIndex(:,ii0)];
        for jj = ii0+1 : size(bVal,2) % cBas.p +1
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
                    [l;obj.levelBas{l}.activeIndex(k+1)];
                   rhs(ind_1) = rhs(ind_1) + elRhs(ii1);
                end
             end
            
            for ll = l : obj.level
                for kk = 0:length(obj.levelBas{ll}.activeIndex)-1 % true?
                    if(ll >1)
                    ind_2 = length(obj.levelBas{ll-1}.activeIndex) +kk+1;
                    else
                    ind_2 = kk+1;
                    end
                    for  ii = 1 : size(bVal,2)
                        for jj = ii : size(bVal,2)
  

                        if(elSInd{ii,jj} == [l ll; obj.levelBas{l}.activeIndex(k+1)...
                                obj.levelBas{ll}.activeIndex(kk+1)])
                            elStiff(ii,jj) = elStiff(ii,jj);
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
end
