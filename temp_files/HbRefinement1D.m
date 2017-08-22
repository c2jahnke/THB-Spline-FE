    function [U, Ubar, Points, Qw] = HbRefinement1D(cBas,fBas,refArea, Points)
    % TESTING
    % carefull: cBas, fBas are references which are changed within this
    % function

    U = cBas.knotVector;
    if(isempty(refArea))
        Ubar = cBas.knotVector;
        Qw = Points;
        return;
    end
    cBas.refArea = refArea;

    cFctStart = cBas.getIndexU(refArea(1)) -1;
    cFctEnd = cBas.getIndexU(refArea(2)) -cBas.p;
    % private attribute is changed, use setter instead!
    % Problem: 1st knot corresponds p-th index
    % Problem: index corresponds to knots, not to basis functions
    % maybe add concept of activeKnots! 
     if(cFctStart == cBas.p-1 && cFctEnd==cBas.n - cBas.p)
         cBas.activeIndex = [];
         cBas.activeKnots = [];
         cBas.setCharM();
     elseif(cFctStart==cBas.p-1) % handle right border
         cBas.activeIndex = [cFctEnd:cBas.n-1];
         cBas.activeKnots = [cBas.knotVector(cFctEnd+cBas.p+1):cBas.knotspan:cBas.knotVector(end-cBas.p)];
         cBas.setCharM();
     elseif(cFctEnd == cBas.n-cBas.p)
         cBas.activeIndex =[0:cFctStart];
         cBas.activeKnots = [cBas.knotVector(cBas.p):cBas.knotspan:cBas.knotVector(cFctStart+2)];
         cBas.setCharM();
     else
         cBas.activeIndex =[0:cFctStart cFctEnd:cBas.n-1]; % active Index w.r.t basis functions
         % active knots with respect to knots!
         cBas.activeKnots = [cBas.knotVector(cBas.p):cBas.knotspan:cBas.knotVector(cFctStart+2) ...
             cBas.knotVector(cFctEnd+cBas.p+1):cBas.knotspan:cBas.knotVector(end-cBas.p)];
         cBas.setCharM();
     end
    fFctStart = fBas.getIndexU(refArea(1)) - 1;
    fFctEnd = fBas.getIndexU(refArea(2)) - fBas.p;
    
     if(fFctStart == cBas.p-1 && fFctEnd == fBas.n-cBas.p)
     elseif(fFctStart == cBas.p-1) %Test
        fBas.activeIndex =0:(fFctEnd-1);
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:fBas.knotVector(fFctEnd+fBas.p+1)];
        fBas.setCharM();
     elseif(fFctEnd == fBas.n-cBas.p) % Test!!
        fBas.activeIndex =(fFctStart+1):fBas.n-1;
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:fBas.knotVector(fFctEnd+fBas.p+1)];
        fBas.setCharM();
     else
        fBas.activeIndex =(fFctStart+1):(fFctEnd-1);
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:...
        fBas.knotVector(fFctEnd+fBas.p+1)];
        fBas.setCharM();
     end
    h0 = cBas.knotspan;
    h1 = fBas.knotspan;
    
    assert(~(refArea(end) - refArea(1) < cBas.p*h1),'Error: Omega1 to small for one basis function.');


    % box aligned refinement
    X = refArea(1)+h1:h0:refArea(2) -h1;%cBas.b-h1; %% not really general, work this out later
    r= size(X,2)-1; % choose only refinement points in refArea
    [Ubar, Qw] = RefineKnotVectCurve(cBas.n,cBas.p,U,Points,X,r); % the curve points Qw are not needed
    
    end
