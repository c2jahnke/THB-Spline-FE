function [THB0, THB1 trunq, q,trP] = ThbRefinement1D(cBas,fBas,refArea,Qw)
    % function evaluation at a small set of points
    % solve basis at p points per element
    % element boundaries plus subdivide element into p subelements
    % we need p(p+1) +1 points per B-Spline
    h = fBas.knotspan/fBas.p; 
    solveVector = fBas.a:h:fBas.b-h; % -h, so that no zero line arises for p = 1
    % 
    trP = Qw;
    V0 = zeros(length(solveVector),cBas.n);
    for ll = 1:cBas.n %cBas.activeIndex % 
        for ks = 1 : length(solveVector)
            V0(ks,ll) = OneBasisFun(cBas.p,cBas.m,cBas.knotVector,ll-1,solveVector(ks));
        end
    end
    
    V1 = zeros(length(solveVector),fBas.n); % adapt size, error 
    for ll = 1:fBas.n %fBas.activeIndex %
        for ks = 1:length(solveVector)
            V1(ks,ll) = OneBasisFun(fBas.p,fBas.m,fBas.knotVector,ll-1,solveVector(ks));
        end
    end

    q = zeros(size(V1,2),size(V0,2)); % adapt size: for p = 1, error on right bound!
    for k = cBas.activeIndex %1 : size(V0,2) %
        q(:,k+1) = V1\V0(:,k+1);
    end
    if(cBas.p == 1 && refArea(2) ~= cBas.b) % temporary fix for p = 1, right border
        q(end,end) = 1;
    end
    % stabilize q
    % Konditionszahlabhängig
    % choosing only the necessary amount of points helps in     
    TOL = 10e-10;
    q(q(:) < TOL) = 0;
    % charactaristic matrix Char0 Char1
    % coefficient matrix Points Qw
    % refinement matrix Q
    % For THB Splines
    trunq = q;% do truncation
    lInd = fBas.getIndexU(refArea(1))+1;
    hInd = fBas.getIndexU(refArea(2))-fBas.p;
    
    trunq(lInd:hInd,:) = 0;

    THB0 = fBas.generBasis*(trunq);
    THB0(THB0(:) < TOL) = 0;
    THB1 = fBas.generBasisRed;
    % check if partition of unity property holds
    for k = 1 : length(fBas.plotVector)
        if( abs(sum(THB0(k,:)) + sum(THB1(k,:))-1) > 10e-3 )
            disp('Warning, sum(THB0(k,:)) + sum(THB1(k,:)) = ');
            disp(sum(THB0(k,:)) + sum(THB1(k,:)) );
            disp(' for k = ');
            disp(k);
        end
    end
    
end
