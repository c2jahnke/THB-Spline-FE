function errorEst(uh,g,obj)
%
error = g(obj.levelBas{1}.plotVector)' - uh; 
error = error(1:end-1);
figure
plot(obj.levelBas{1}.plotVector(1:end-1),error,'b')
ylim([min(error) max(error)])
ylabel("absolute error");
end

