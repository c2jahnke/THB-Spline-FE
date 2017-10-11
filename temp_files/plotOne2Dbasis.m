    function plotOne2Dbasis(p,m,knotVector,iU,iV,plotVector,sP,resol)
        [X,Y] = meshgrid(knotVector(1):resol:knotVector(end),knotVector(1):resol:knotVector(end));

        CU = zeros(sP,1);
        for kk = 1 : sP
        CU(kk,1) = OneBasisFun(p,m,knotVector,iU,plotVector(kk));

        end
        CV = zeros(sP,1);
        for kk = 1 : sP
        CV(kk,1) = OneBasisFun(p,m,knotVector,iV,plotVector(kk));

        end
        Z = CU(:,1)*CV(:,1)'; % tensor product structure
        surf(X,Y,Z)
        %hold off
    end