
    function ImpointCurvePlot(n,sP,p,U,plotVector,Points)

        C = zeros( sP , 2 );
        for l = 1:sP % work around, make sure other points are not needed
            C(l,:) = CurvePoint(n,p,U,Points,plotVector(l));
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
    
    
