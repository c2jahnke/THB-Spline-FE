function HbRef1D(obj,rLevel)
 %--------------------------------------------------------------
 % activates indices between rLevel and rLevel + 1
 %--------------------------------------------------------------
      
    cBas = obj.levelBas{rLevel};
    fBas = obj.levelBas{rLevel+1};
    refArea = obj.levelBas{rLevel}.refArea;
    if(rLevel > 1)
        gBas = obj.levelBas{rLevel -1}; % ground base
    end

    if(isempty(refArea))
        fBas.activeIndex = [];
        return;
    end


    cFctStart = cBas.getIndexU(refArea(1)) -1; 
    cFctEnd = cBas.getIndexU(refArea(2)) -cBas.p;
    
    fFctStart = fBas.getIndexU(refArea(1)) - 1;
    fFctEnd = fBas.getIndexU(refArea(2)) - fBas.p;

      if(rLevel >1) % middle base Works
         if((cBas.getIndexU(gBas.refArea(1))-cBas.p)==0 && cBas.getIndexU(gBas.refArea(2))-cBas.p==cBas.n-cBas.p)
         % Case1: refinement until both boundary
             cBas.activeIndex =[0:cFctStart cFctEnd:cBas.n-1];
         elseif((cBas.getIndexU(gBas.refArea(1))-cBas.p)==0)
             % Case2: refinement until left boundary
            cBas.activeIndex = [0:cFctStart ...
             cFctEnd:(cBas.getIndexU(gBas.refArea(2))-cBas.p -1)];
         elseif(cBas.getIndexU(gBas.refArea(2))-cBas.p==cBas.n-cBas.p)
             % Case3: refinement until right boundary
            cBas.activeIndex = [cBas.getIndexU(gBas.refArea(1)):cFctStart ...%2*(gBas.getIndexU(gBas.refArea(1)))-cBas.p
             cFctEnd:cBas.n-1];
         else 
         cBas.activeIndex = [cBas.getIndexU(gBas.refArea(1)):cFctStart ...
             cFctEnd:(cBas.getIndexU(gBas.refArea(2))-cBas.p -1)];
         end
     elseif(rLevel ==1) %ground base Works
        if(cFctStart == cBas.p-1 && cFctEnd==cBas.n - cBas.p)
         cBas.activeIndex = [];
        elseif(cFctStart==cBas.p-1) % handle right border
         cBas.activeIndex = [cFctEnd:cBas.n-1];

        elseif(cFctEnd == cBas.n-cBas.p)
         cBas.activeIndex =[0:cFctStart];
        else 
            cBas.activeIndex =[0:cFctStart cFctEnd:cBas.n-1]; % active Indices w.r.t basis functions
            % active knots
             %cBas.activeKnots = refArea(1):cBas.knotspan:refArea(2); 
        end
    end
      %finest base, %Works
    if(fFctStart == cBas.p-1 && fFctEnd == fBas.n-cBas.p)
         fBas.activeIndex = [0 : fBas.n-1];
         %fBas.activeKnots = fBas.getAllKnots;
        elseif(fFctStart == cBas.p-1) %Test
        fBas.activeIndex =0:(fFctEnd-1);
 
        elseif(fFctEnd == fBas.n-cBas.p) % Test!!
        fBas.activeIndex =(fFctStart+1):fBas.n-1;
        else %%
        fBas.activeIndex =(fFctStart+1):(fFctEnd-1);
    end
    cBas.setCharM(); % reset characteristic matrices
    fBas.setCharM();
    
    assert(~(refArea(end) - refArea(1) < cBas.p*fBas.knotspan),'Error: Omega1 to small for one basis function.');
    
end
