function solVal = evalSol(obj,u,iBasisFctInd,iLvl,y)
basVal= obj.evalBasis(u);
[uLvl,uInd] = obj.getActiveFctIndU(u);
sULvl = length(uLvl);
%u_ord = reorderU(obj,u,iLvl,iBasisFctInd)
yS = [];
for k = 1 : obj.nOF
    for l = 1 : sULvl
        if ( isequal(iBasisFctInd(k),uInd(l)) && isequal(iLvl(k),uLvl(l)) )
        yS = [yS y(k)];
        end
    end
%    disp('Warning: evaluation of solution not possible.');
end
solVal = basVal*yS';
end