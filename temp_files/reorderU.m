function u_ord = reorderU(obj,u,nOF,iLvl,iBasisFctInd)

% reorder solution u according to basis function appearence
 u_ord = zeros(size(u));
 [lvl, BasisFctInd] = obj.getAllActiveFct;
 for k = 1 : nOF
     for l = 1 : nOF
         if(lvl(k) == iLvl(l) && BasisFctInd(k) == iBasisFctInd(l) )  
             u_ord(k+1) = u(l+1);
         end
     end 
 end
y =  u_ord(2:end+1-obj.levelBas{1}.p);

end