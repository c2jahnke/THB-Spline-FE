function ThbRef1D(obj,rLevel)
    % function evaluation at a small set of points
    % solve basis at p points per element
    % element boundaries plus subdivide element into p subelements
    % we need p(p+1) +1 points per B-Spline
    cBas = obj.levelBas{rLevel};
    fBas = obj.levelBas{rLevel+1};
    refArea = obj.levelBas{rLevel}.refArea;
    h = fBas.knotspan/fBas.p; 
    solveVector = fBas.a:h:fBas.b-h; % -h, so that no zero line arises for p = 1
    
    A0 = zeros(length(solveVector),cBas.n);
    for ll = 1:cBas.n %cBas.activeIndex % 
        for ks = 1 : length(solveVector)
            A0(ks,ll) = OneBasisFun(cBas.p,cBas.m,cBas.knotVector,ll-1,solveVector(ks));
        end
    end
    
    A1 = zeros(length(solveVector),fBas.n); % adapt size, error 
    for ll = 1:fBas.n %fBas.activeIndex %
        for ks = 1:length(solveVector)
            A1(ks,ll) = OneBasisFun(fBas.p,fBas.m,fBas.knotVector,ll-1,solveVector(ks));
        end
    end
    cLind = cBas.getIndexU(refArea(1))+1
    cHind = cBas.getIndexU(refArea(2))-cBas.p
    q = zeros(size(A1,2),size(A0,2)); % adapt size: for p = 1, error on right bound!
    for k = cBas.activeIndex %1 : size(A0,2) %
        q(:,k+1) = A1\A0(:,k+1);
    end
    if(cBas.p == 1 && refArea(2) ~= cBas.b) % fix for p = 1, right border
        q(end,end) = 1;
    end
    % stabilize q
    % Konditionszahlabhängig
    % choosing only the necessary amount of points helps in     
    TOL = 10e-10;
    q(q(:) < TOL) = 0;
    % For THB Splines
    trunq = q;
    lInd = fBas.getIndexU(refArea(1))+1;
    hInd = fBas.getIndexU(refArea(2))-fBas.p;

    trunq(lInd:hInd,:) = 0; % do truncation

    redTrunqDiff = abs(sum(trunq -q));
    truncIndex = find(redTrunqDiff)-1;
%     truncIndex = [];
%     for k = 1 : cBas.n
%         if(redTrunqDiff(k) > 0 )
%             truncIndex = [truncIndex k-1];
%         end
%     end
    cBas.truncIndex = truncIndex;
    cBas.trunc = trunq;
end
