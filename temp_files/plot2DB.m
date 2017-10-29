%% plot THB basis
a = 0;
b = 10;
p = 2;
N = 5;
resol = 0.5;
obj = thbSplBasML(0,6,2,3,0.1,2)
%obj.ThbRefinement1DML(1,[2,8])

for k = 1
    for l = obj.levelBas{k}.activeIndex
        plot2DTHB(obj,k,l)
        hold on;
    end
    xlabel('x')
ylabel('y')
zlabel('value');
end
figure
plotOne2DTHB(obj,1,2)
hold on
obj2 = thbSplBasML(0,6,2,3,0.001,2)
objFun2 = bSplBasFun(2,obj2.levelBas{1});
%objFun2.plotOneBasisFun(objFun2.generOneBasisFun)
C2 = objFun2.generOneBasisFun;
scatter3(objFun2.plotVector,zeros(objFun2.sP,1),C2,'k.')
scatter3(zeros(objFun2.sP,1),objFun2.plotVector,C2,'r.')
xlabel('x')
ylabel('y')
zlabel('value');