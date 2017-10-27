
function [Stiffn, rhs, iLvl,iBasisFctInd] = assemblePoissThbMl(obj,f)

allKnots = obj.getAllKnots;
   
nOF = obj.nOF; % number of functions
nOE =length(allKnots) -1; % number of elements
Stiffn = zeros(nOF); %basis.n
ngp = max(obj.levelBas{1}.p+1,sqrt(obj.levelBas{1}.p^2 -2*obj.levelBas{1}.p+1));
% first part to integrate right hand side sufficiently accurate, second
% part to generate stiffness matrix exactly
rhs = zeros(nOF,1); % number of basis functions!


for el = 1 : nOE % loop over elements
    %% export evaluation to seperate function
    [bVal,gradVal,lIndex,s,w] = evalEl(obj,el);
    % assemble element stiffness matrix
    elRhs = zeros(size(bVal,2),1); 
    elStiff = zeros(size(bVal,2));
    elSInd = cell(size(bVal,2));
    elRInd =cell(size(bVal,2),1);
    for ii0 = 1 : size(bVal,2) 
        elRhs(ii0) = sum(w.*arrayfun(f,(s)).*bVal(:,ii0)); % changed to arrayfun for non vectorized functions
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
    
    % rewrite assembly for arbitrary levels
    for l = 1 : obj.level
        for k = 1:length(obj.levelBas{l}.activeIndex)
            if(l >1) % raise index according to lower levels
                ind_1 = 0;
                for tInd = 1 : l-1
                    ind_1 = ind_1 + length(obj.levelBas{tInd}.activeIndex);
                end
                ind_1 = ind_1 + k;
            else
                ind_1 = k; % index for basis functions
            end
            for ii1 = 1 : size(bVal,2)
                if(elRInd{ii1} == [l;obj.levelBas{l}.activeIndex(k)])
                   rhs(ind_1) = rhs(ind_1) + elRhs(ii1);
                   iLvl(ind_1) = l; % keep track of all indices
                   iBasisFctInd(ind_1) = obj.levelBas{l}.activeIndex(k);
                end
             end
            
            for ll = l : obj.level
                for kk = 1:length(obj.levelBas{ll}.activeIndex) 
                    if(ll >1) % raise index according to lower levels
                        ind_2 = 0;
                    for tInd = 1 : ll-1
                        ind_2 = ind_2 + length(obj.levelBas{tInd}.activeIndex);
                    end
                    ind_2 = ind_2 + kk;
                    else
                    ind_2 = kk;
                    end
                    for  ii = 1 : size(bVal,2)
                        for jj = ii : size(bVal,2)
                        if(elSInd{ii,jj} == [l ll; obj.levelBas{l}.activeIndex(k)...
                                obj.levelBas{ll}.activeIndex(kk)])
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
