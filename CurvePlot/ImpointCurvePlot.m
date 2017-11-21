
function ImpointCurvePlot(n,sP,p,U,plotVector,Points)
    % function to plot a curve with control points
 
    C = zeros( sP , 2 );
    for l = 1:sP 
        C(l,:) = CurvePoint(n,p,U,Points,plotVector(l));
    end
    plt_curve = plot(C(:,1),C(:,2),'k','LineWidth',1.2);
    hold on;
    xlabel('x');
    ylabel('y');

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
    
    
