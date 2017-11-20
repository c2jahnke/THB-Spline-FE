    function plot2DTHB(obj,lvl,ind)
        [X,Y] = meshgrid(obj.levelBas{lvl}.knotVector(1):obj.levelBas{lvl}.resol...
            :obj.levelBas{lvl}.knotVector(end),obj.levelBas{lvl}.knotVector(1)...
            :obj.levelBas{lvl}.resol:obj.levelBas{lvl}.knotVector(end));
        objFun1 = thbSplBasFun(ind,obj,lvl);
        CU = objFun1.generOneBasisFun();
        objFun1. plotOneBasisFun(CU)
        for k =0 : obj.levelBas{lvl}.n-1
            objFun2 = bSplBasFun(k,obj.levelBas{lvl});
            CV = objFun2.generOneBasisFun();
        Z = CU(:,1)*CV(:,1)'; % tensor product structure
        %surf(X(1:end-20,(1:end-20)),Y(1:end-20,(1:end-20)),Z(1:end-20,(1:end-20)))
%         if(k == 3 && ind ~=3)
%             surf(Y,X,Z);
%         end
        surf(X,Y,Z);
        hold on;
        end
        %objFun2 = bSplBasFun(ind,obj.levelBas{1});
        %objFun2.plotOneBasisFun(objFun2.generOneBasisFun)
        %C2 = objFun2.generOneBasisFun;
        %plot(objFun2.plotVector,C2)
        
    end