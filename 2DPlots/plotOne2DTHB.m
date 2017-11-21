function plotOne2DTHB(obj,lvl,ind)
        [X,Y] = meshgrid(obj.levelBas{lvl}.knotVector(1):obj.levelBas{lvl}.resol...
            :obj.levelBas{lvl}.knotVector(end),obj.levelBas{lvl}.knotVector(1)...
            :obj.levelBas{lvl}.resol:obj.levelBas{lvl}.knotVector(end));
        objFun1 = thbSplBasFun(ind,obj,lvl);
        CU = objFun1.generOneBasisFun();
        CV = CU;
        Z = CU(:,1)*CV(:,1)'; % tensor product structure
        %surf(X(1:end-40,(1:end-40)),Y(1:end-40,(1:end-40)),...
        %    Z(1:end-40,(1:end-40)));
       surf(X,Y,Z);
end