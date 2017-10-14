function [U, Ubar, Points, Qw] = HbRef1D(obj,rLevel, Points)
    % carefull: cBas, fBas are references which are changed within this
    % function
    cBas = obj.levelBas{rLevel};
    fBas = obj.levelBas{rLevel+1};
    refArea = obj.levelBas{rLevel}.refArea;
    if(rLevel > 1)
        gBas = obj.levelBas{rLevel -1}; % ground base
    end
    U = cBas.knotVector;
    if(isempty(refArea))
        Ubar = cBas.knotVector;
        Qw = Points;
        fBas.activeIndex = [];
        return;
    end
    cBas.refArea = refArea;

    cFctStart = cBas.getIndexU(refArea(1)) -1;
    cFctEnd = cBas.getIndexU(refArea(2)) -cBas.p;
    % private attribute is changed, use setter instead!
    % Problem: 1st knot corresponds p-th index
    % Problem: index corresponds to knots, not to basis functions
     if(rLevel >1) % whar are the different cases?
         if((cBas.getIndexU(gBas.refArea(1))-cBas.p)==0 && cBas.getIndexU(gBas.refArea(2))-cBas.p==cBas.n-cBas.p)
         % Case1: refinement until both boundary
             cBas.activeIndex =[0:cFctStart cFctEnd:cBas.n-1];
         elseif((cBas.getIndexU(gBas.refArea(1))-cBas.p)==0)
             % Case2: refinement until left boundary
            cBas.activeIndex = [0:cFctStart ...
             cFctEnd:(cBas.getIndexU(gBas.refArea(2))-cBas.p -1)];
         elseif(cBas.getIndexU(gBas.refArea(2))-cBas.p==cBas.n-cBas.p)
             % Case3: refinement until right boundary
            cBas.activeIndex = [cBas.getIndexU(gBas.refArea(1)):cFctStart...%2*(gBas.getIndexU(gBas.refArea(1)))-cBas.p
             cFctEnd:cBas.n-1];
         else %cBas.getIndexU(gBas.refArea(1))
         cBas.activeIndex = [cBas.getIndexU(gBas.refArea(1)):cFctStart...
             cFctEnd:(cBas.getIndexU(gBas.refArea(2))-cBas.p -1)];
         end
     elseif(rLevel ==1)
        if(cFctStart == cBas.p-1 && cFctEnd==cBas.n - cBas.p)
         cBas.activeIndex = [];
         cBas.activeKnots = [];
     %    cBas.setCharM();
        elseif(cFctStart==cBas.p-1) % handle right border
         cBas.activeIndex = [cFctEnd:cBas.n-1];
         cBas.activeKnots = [cBas.knotVector(cFctEnd+cBas.p+1):cBas.knotspan:cBas.knotVector(end-cBas.p)];
    %     cBas.setCharM();
        elseif(cFctEnd == cBas.n-cBas.p)
         cBas.activeIndex =[0:cFctStart];
         cBas.activeKnots = [cBas.knotVector(cBas.p):cBas.knotspan:cBas.knotVector(cFctStart+2)];
        else 
            cBas.activeIndex =[0:cFctStart cFctEnd:cBas.n-1]; % active Indices w.r.t basis functions
            % active knots
             cBas.activeKnots = refArea(1):cBas.knotspan:refArea(2); %[cBas.knotVector(cBas.p):cBas.knotspan:cBas.knotVector(cFctStart+2) ...
            % cBas.knotVector(cFctEnd+cBas.p+1):cBas.knotspan:cBas.knotVector(end-cBas.p)];
        end
        % cBas.setCharM();
     end
     cBas.setCharM();
     
    fFctStart = fBas.getIndexU(refArea(1)) - 1;
    fFctEnd = fBas.getIndexU(refArea(2)) - fBas.p;
    
     if(fFctStart == cBas.p-1 && fFctEnd == fBas.n-cBas.p)
         fBas.activeIndex = [0 : fBas.n-1];
         fBas.activeKnots = fBas.getAllKnots;
     elseif(fFctStart == cBas.p-1) %Test
        fBas.activeIndex =0:(fFctEnd-1);
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:fBas.knotVector(fFctEnd+fBas.p+1)];

     elseif(fFctEnd == fBas.n-cBas.p) % Test!!
        fBas.activeIndex =(fFctStart+1):fBas.n-1;
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:fBas.knotVector(fFctEnd+fBas.p+1)];
     else
        fBas.activeIndex =(fFctStart+1):(fFctEnd-1);
        fBas.activeKnots = [fBas.knotVector(fFctStart+2):fBas.knotspan:...
        fBas.knotVector(fFctEnd+fBas.p+1)];
     end
             fBas.setCharM();
    h0 = cBas.knotspan;
    h1 = fBas.knotspan;
    
    assert(~(refArea(end) - refArea(1) < cBas.p*h1),'Error: Omega1 to small for one basis function.');


    % box aligned refinement
    X = refArea(1)+h1:h0:refArea(2) -h1;%cBas.b-h1; %% not really general, work this out later
    r= size(X,2)-1; % choose only refinement points in refArea
    [Ubar, Qw] = RefineKnotVectCurve(cBas.n,cBas.p,U,Points,X,r); % the curve points Qw are not needed
    
end
