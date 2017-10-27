function [H1err,L2err,Enerr] = errCalc(obj,g,y,iBasisFctInd,iLvl)
%% Calculate H1,L2,Energy norm
L2err = 0; H1err = 0; Enerr = 0;
gDer = @(x) (g(x+obj.levelBas{1}.resol) -g(x))/obj.levelBas{1}.resol;
for el = 1:obj.nOE
    %[bVal,gradVal,lIndex,s,w] = evalEl(obj,el);
    allKnots = obj.getAllKnots;
    ngp = 10;
    [s,w]=lgwt(ngp,allKnots(el),allKnots(el+1));
    uexVal = feval(g,s);
    %calculate dervative using central difference scheme, upwind scheme

    uexGrad = feval(gDer, s);
    
    
    uappGrad = [];
    for k = 1 : length(s)
        approx = evalDersSol(obj,s(k),iBasisFctInd,iLvl,y);
        uappGrad = [uappGrad approx];
    end
    L2err = L2err + ((uappGrad(1,:)-uexVal(:)').^2*w);%L2err + (h/6)*((ue(1)-q(i))^2+4*(ue(2)-0.5*(q(i)+q(i+1)))^2+(ue(3)-q(i+1))^2);
    Enerr = Enerr +((uappGrad(2,:)- uexGrad(:)').^2*w);%Enerr + (h/6)*((udxe(1)-qdx)^2+4*(udxe(2)-qdx)^2+(udxe(3)-qdx)^2);
end
H1err = sqrt(L2err + Enerr);
L2err = sqrt(L2err);
Enerr = sqrt(Enerr);
    

end