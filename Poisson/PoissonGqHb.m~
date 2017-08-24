%clc, clear all, close all
a = 0;
b = 10;
p = 1;
N = 5;
resol = 0.001;
lvl = 2;
obj = hbSplBasML(a,b,p,N,resol,lvl);
f = @(x) (x);%.^(1/31); %(x).^(1/3);
    % change generation of Points!
Points = zeros(obj.levelBas{1}.n,2); % not needed
Points(:,1) = linspace(0,10, size(Points,1) );
Points(:,2) = f(2*pi*Points(:,1));


cBas = obj.levelBas{1};
fBas = obj.levelBas{2};
refArea = [4 10];
[U, Ubar, Points, Qw] = HbRefinement1D(cBas,fBas,refArea,Points);
cBas = obj.levelBas{1};
fBas = obj.levelBas{2};

if(isempty(refArea))
    nOE = N;
    nOF = cBas.n;
else
nOE =N+(refArea(2) - refArea(1))/cBas.knotspan; % number of elements
nOF = length(cBas.activeIndex) +length(fBas.activeIndex); % number o functions
end
Stiffn = zeros(length(cBas.activeIndex) +length(fBas.activeIndex)); %basis.n % length(cBas.activeIndex) +length(fBas.activeIndex) 
elStiff = zeros(cBas.p +2); % 
% Ax = b
% Lapl(u) = f
ngp = obj.levelBas{1}.p+3;

rhs = zeros(nOF,1); % basis.n number of basis functions!
elRhs =  zeros(cBas.p+2,1);
%allKnots = basis.getAllKnots;
allKnots = unique(obj.getAllKnots);
flag = true;

% number of elements: N+(refArea(2) - refArea(1))/cBas.knotspan
for el = 1 : nOE % loop over elements
    [s,w]=lgwt(ngp,allKnots(el),allKnots(el+1));
    bVal = []; % zeros(ngp, basis.p+1); % basis evaluation
    gradVal = []; % zeros(ngp,basis.p+1); % derivative evaluation
    for j = length(s):-1:1 %basis.evalDersBasis(s(j))% not yet fully tested!
        temp = obj.evalDersBasis(s(j))%DersBasisFuns(k,s(j),basis.p,basis.knotVector,1) % replace by method of class!
        [lvl, Ind] = obj.getActiveFctIndU(s(j));
        lIndex = [lvl ; Ind]
        % make sure, the order of the basis functions is correct!

        bVal(j,:) = temp(1,:);
        gradVal(j,:) = temp(2,:); % 
    end
    
    elRhs = zeros(size(bVal,2),1); % to large!
    elStiff = zeros(size(bVal,2));
    elSInd = cell(size(bVal,2));
    elRInd =cell(size(bVal,2),1);
    for ii = 1 : size(bVal,2) % cBas.p+1
        elRhs(ii) = sum(w.*f(s).*bVal(:,ii));
        elRInd{ii} = lIndex(:,ii);
        elStiff(ii,ii) = sum(w.*gradVal(:,ii).^2);
        elSInd{ii,ii} = [lIndex(:,ii) lIndex(:,ii)];
        for jj = ii+1 : size(bVal,2) % cBas.p +1
            elStiff(ii,jj) = sum(w.*gradVal(:,ii).*gradVal(:,jj));
            elSInd{ii,jj} = [lIndex(:,ii) lIndex(:,jj)];
            elStiff(jj,ii) = elStiff(ii,jj);
            elSInd{jj,ii} = elSInd{ii,jj};
        end      
    end
    
    

    for l = 1 : obj.level
        for k = 0:length(obj.levelBas{l}.activeIndex)-1
            if(l >1)
            ind_1 = length(obj.levelBas{l-1}.activeIndex) +k+1;
            else
                ind_1 = k+1;
            end
            for ii = 1 : size(bVal,2)
                if(elRInd{ii} == [l;obj.levelBas{l}.activeIndex(k+1)])
                   rhs(ind_1) = elRhs(ii);
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
                         [l ll; obj.levelBas{l}.activeIndex(k+1)...
                                obj.levelBas{ll}.activeIndex(kk+1)]

                        if(elSInd{ii,jj} == [l ll; obj.levelBas{l}.activeIndex(k+1)...
                                obj.levelBas{ll}.activeIndex(kk+1)])
                            elSInd{ii,jj} 
                            elStiff(ii,jj) = elStiff(ii,jj)
                            Stiffn(ind_1,ind_2) = Stiffn(ind_1,ind_2) + elStiff(ii,jj);
                            Stiffn(ind_2,ind_1) = Stiffn(ind_1,ind_2)
                        end
                        end
                    end
                end 
            end
        end
    end
    
end
Stiffn

%% BC with Lagrange multipliers
A = zeros(nOE+obj.levelBas{1}.p +1);
b = zeros(nOE+obj.levelBas{1}.p +1,1);
A(1,2) = 1;
A(2,1) = 1;
A(2:end,2:end) = Stiffn;
b(1,1) = 0; % Dichlet BC
b(2:end) = rhs;

u = A\b; % use sum u * basis.functions as a curve,
y = u(2:end);
x = linspace(cBas.a,cBas.b,size(y,1));
Points = [x' y];
% define as method of the class bSplBas, input parameter Points
hbKnotVector = allKnots;
for k = 1 : cBas.p
    hbKnotVector = [cBas.a hbKnotVector cBas.b];
end % cBas.n
%% change so that hb-basis is used
ImpointCurvePlot(nOF,cBas.sP,cBas.p,hbKnotVector,cBas.plotVector,Points)



