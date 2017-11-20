%% plot THB basis
a = 0;
b = 8;
p = 2;
N = 4;
resol = 0.5;
obj = thbSplBasML(a,b,p,N,0.1,2);
%obj.plotBas
obj.ThbRefinement1DML(1,[0,4])
%obj.plotBas
Ind = obj.levelBas{1}.getIndices;
passIndex = setdiff(Ind.basisFunctionIndex, obj.levelBas{k}.activeIndex)
for k = 1:obj.level
    for l = obj.levelBas{k}.activeIndex
        plot2DTHB(obj,k,l)
        hold on;
    end
    for l = passIndex
        plot2DTHBRef(obj,k,l)
        hold on
    end
    
    xlabel('x')
    ylabel('y')
    zlabel('value');
end