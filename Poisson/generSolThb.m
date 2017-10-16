function uh =  generSol2(obj,y)
% fix the bug!
D = [];
for ll = 1 : obj.level
    C = zeros(obj.levelBas{ll}.sP,length(obj.levelBas{ll}.activeIndex));
    for l = 1:length(obj.levelBas{ll}.activeIndex)
        % create basis function for thb-splines
        basisF = thbSplBasFun(obj.levelBas{ll}.activeIndex(l),obj,ll);
        C(:,l) = basisF.generOneBasisFun;
    end
    D = [D C(1:2^(ll-1):size(C,1),:)];% D = [D C]
end

uh = y(1)*D(:,1);
for k = 2 : obj.nOF
    uh = uh + y(k)*D(:,k);
end

%plot(obj.levelBas{1}.plotVector(1:end-1),uh(1:end-1),'r')

end