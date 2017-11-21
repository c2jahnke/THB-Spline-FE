function [bVal,gradVal,lIndex,s,w] = evalEl(obj,el)
% used in assemblePoissThbMl, errEst
    ngp = max(obj.levelBas{1}.p+1,sqrt(obj.levelBas{1}.p^2 -2*obj.levelBas{1}.p+1));
    allKnots = obj.getAllKnots;
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

end