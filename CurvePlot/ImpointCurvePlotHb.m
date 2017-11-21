function ImpointCurvePlotHb(obj,n,sP,p,U,plotVector,Points)
% extension to HB-splines, TEST!
    refArea = obj.levelBas{1}.refArea;
    area = abs(obj.levelBas{1}.a - obj.levelBas{1}.refArea(1)) + abs(obj.levelBas{1}.b - obj.levelBas{1}.refArea(2));
    totArea = abs(obj.levelBas{1}.b- obj.levelBas{1}.a);
    lambda = area/totArea;
    sC = lambda*obj.levelBas{1}.sP+(1-lambda)*obj.levelBas{2}.sP;
    C = zeros(sC, 2 );
    QPoints = genQPoints(obj,Points,2);
    PPoints = genQPoints(obj,Points,1);

    % does not work so far!
    % we want to calculate the curve piecewise
    a = obj.levelBas{1}.a;
    b = obj.levelBas{1}.b;
    % how can P1 be set in a general way?
    P1 = Points(1:2,:);
    N1 = (refArea(1) -a)/obj.levelBas{1}.knotspan;
    basis1 = bSplBas(a,refArea(1),p,N1,0.1);
    ImpointCurvePlot(basis1.n,basis1.sP,basis1.p,basis1.knotVector,basis1.plotVector,P1)

    P2 = Points(2:end,:)
    N2 = (refArea(2) -refArea(1))/obj.levelBas{2}.knotspan;
    basis2 = bSplBas(refArea(1),refArea(2),p,N2,0.1);
    ImpointCurvePlot(basis2.n,basis2.sP,basis2.p,basis2.knotVector,basis2.plotVector,P2)



    for l = 1:refArea(1)/obj.levelBas{1}.resol 
        C(l,:) = CurvePoint(obj.levelBas{1}.n,obj.levelBas{1}.p,obj.levelBas{1}.knotVector,PPoints,obj.levelBas{1}.plotVector(l));

    end
    for l = refArea(1)/obj.levelBas{1}.resol+1: refArea(2)/obj.levelBas{1}.resol
        C(l,:) = CurvePoint(obj.levelBas{2}.n,obj.levelBas{2}.p,obj.levelBas{2}.knotVector,QPoints,obj.levelBas{2}.plotVector(refArea(1)/obj.levelBas{1}.resol + l));
    end
    for l = refArea(2)/obj.levelBas{1}.resol+1 : obj.levelBas{1}.b/obj.levelBas{1}.resol
        C(l,:) = CurvePoint(obj.levelBas{1}.n,obj.levelBas{1}.p,obj.levelBas{1}.knotVector,PPoints,obj.levelBas{1}.plotVector(l));
    end

    plt_curve = plot(C(:,1),C(:,2),'r.');
    hold on;

    Impoints = {};
    for i = 1:n
       Impoints{i} = impoint(gca, Points(i,1) , Points(i,2) ); 
       addNewPositionCallback( Impoints{i} , @(newpos) updateCurvePlot(i,newpos) );
    end

    function updateCurvePlot( index , newpos )
        Points(index,:) = newpos;
         for k = 1:sP
             C(k,:) = CurvePoint(n,p,U,Points,plotVector(k));
         end
        plt_curve.XData = C(:,1);
        plt_curve.YData = C(:,2);
        drawnow;
    end
end


    
