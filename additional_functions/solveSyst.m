function y = solveSyst(obj,Stiffn,rhs,iLvl,iBasisFctInd,dBC)
% maybe put this to other file/function /merge with assembly
[edLvl, edInd] = obj.getActiveFctIndU(obj.levelBas{1}.a);
eEval = obj.evalBasis(obj.levelBas{1}.a);
eIndEval = find(eEval);
eFLvl = find( edLvl(eIndEval(1)) == iLvl);
eFInd = find(edInd(eIndEval(1)) == iBasisFctInd);
sInd = intersect(eFLvl,eFInd); % start index

[edLvl, edInd] = obj.getActiveFctIndU(obj.levelBas{1}.b);
eEval = obj.evalBasis(obj.levelBas{1}.a);
eIndEval = find(eEval);
eFLvl = find( edLvl(eIndEval(1)) == iLvl);
eFInd = find(edInd(eIndEval(1)) == iBasisFctInd);
eInd = intersect(eFLvl,eFInd); % end index

if(dBC.mixed == true)
    A = zeros(obj.nOF +1);
    b = zeros(obj.nOF +1,1);
    if(strcmp(dBC.east,'Dirichlet') && strcmp(dBC.west,'Neumann'))
        % does not work this way!
        A(end,eInd-1) = 1;
        A(eInd-1,end) = 1;
        b(end,1) = dBC.eVal; % Dichlet BC
        A(1:end-1,1:end-1) = Stiffn;
        b(1:end-1) = rhs;
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

    
elseif(strcmp(dBC.west,'Neumann'))
    y = Stiffn\b;

end

end