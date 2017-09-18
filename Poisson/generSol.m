function uh =  generSol(obj,y)
D = [];
for ll = 1 : obj.level
    C = zeros(obj.levelBas{ll}.sP,length(obj.levelBas{ll}.activeIndex));
    for l = 1:length(obj.levelBas{ll}.activeIndex)
        basisF = bSplBasFun(obj.levelBas{ll}.activeIndex(l),obj.levelBas{ll});
        C(:,l) = basisF.generOneBasisFun;
    end
    D = [D C(1:ll:size(C,1),:)];
end

uh = y(1)*D(:,1);
for k = 2 : size(D,2)
    uh = uh + y(k)*D(:,k);
end

plot(obj.levelBas{1}.plotVector,uh,'r')
   
end