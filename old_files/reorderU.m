function u_ord = reorderU(obj,u,iLvl,iBasisFctInd)

% reorder solution u according to basis function appearence
 u_ord = zeros(size(u));
 [lvl, BasisFctInd] = obj.getAllActiveFct;
 for k = 1 : obj.nOF
     for l = 1 : obj.nOF
         if(lvl(k) == iLvl(l) && BasisFctInd(k) == iBasisFctInd(l) )  
             u_ord(k+1) = u(l+1);
         end
     end 
 end


end