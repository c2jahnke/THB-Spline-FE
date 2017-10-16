function y = solveSyst(obj,Stiffn,rhs,iLvl,iBasisFctInd,dBC)
% maybe put this to other file/function /merge with assembly
[sdLvl, sdInd] = obj.getActiveFctIndU(obj.levelBas{1}.a);
sEval = obj.evalBasis(obj.levelBas{1}.a);
sIndEval = find(sEval);
sFLvl = find( sdLvl(sIndEval(1)) == iLvl);
sFInd = find(sdInd(sIndEval(1)) == iBasisFctInd);
sInd = intersect(sFLvl,sFInd); % start index

% not nice, but works
[edLvl, edInd] = obj.getActiveFctIndU(obj.levelBas{1}.b-obj.levelBas{1}.resol);
eEval = obj.evalBasis(obj.levelBas{1}.b-obj.levelBas{1}.resol);
eIndEval = find(eEval);
eFLvl = find( edLvl(eIndEval(end)) == iLvl); % 
eFInd = find(edInd(eIndEval(end)) == iBasisFctInd);
eInd = intersect(eFLvl,eFInd); % end index

if(dBC.mixed == true)
    A = zeros(obj.nOF +1);
    b = zeros(obj.nOF +1,1);
    if(strcmp(dBC.east,'Dirichlet') && strcmp(dBC.west,'Neumann'))
        % does not work this way!
        A(1,eInd+1) = 1;
        A(eInd+1,1) = 1;
        b(1,1) = dBC.eVal; % Dichlet BC
       A(2:end,2:end) = Stiffn;
        b(2:end) = rhs;
        u = A\b; 
        y = u(2:end);   
    elseif(strcmp(dBC.west,'Dirichlet') && strcmp(dBC.east,'Neumann'))
        A(1,sInd+1) = 1;
        A(sInd+1,1) = 1;
        b(1,1) = dBC.wVal; % Dichlet BC
        A(2:end,2:end) = Stiffn;
        b(2:end) = rhs;
        u = A\b; 
        y = u(2:end);
    end
elseif(strcmp(dBC.east,'Dirichlet'))
    A = zeros(obj.nOF +2);
    b = zeros(obj.nOF +2,1);
    A(1,sInd+2) = 1;
    A(sInd+2,1) = 1;
    b(1,1) = dBC.wVal;
    A(2,eInd+2) = 1;
    A(eInd+2,2) = 1;
    b(2,1) = dBC.eVal;
    
    A(3:end,3:end) = Stiffn;
    b(3:end) = rhs;
    u = A\b; 
    y = u(3:end);
elseif(strcmp(dBC.west,'Neumann'))
    y = Stiffn\b;

end

end