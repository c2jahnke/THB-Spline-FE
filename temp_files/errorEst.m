function errorEst(uh,g,obj)

%
error = g(obj.levelBas{1}.plotVector)' - uh; 
figure
plot(obj.levelBas{1}.plotVector,error,'b')
ylim([-0.1 0.1])
ylabel("absolute error");

end

