function PoissGqThbMlSin(nEL,p,fileName)
% case study for Poisson equation with sine source term
a = 0;
b = 2;
%p = 2;
N = nEL;
resol = 0.001;
lvl = 2;%parameter.a,parameter.b,parameter.p,parameter.N,parameter.resol,parameter.lvl
obj = thbSplBasML(a,b,p,N,resol,lvl);
f = @(x)  pi.^2/4*sin(pi/2*x);%source function 
g = @(x) sin(pi/2*x);

dBC = boundCond('Dirichlet','Dirichlet',0,0);% Dichlet BC
% refArea = [];
refArea = [0.5 1.5]; %refArea = [4/3 8/3];
obj.ThbRefinement1DML(1,refArea);
refArea2 = [2/3 7/3];%[3*pi/8 7*pi/8];
% obj.ThbRefinement1DML(2,refArea2);
% refArea3 = [4*pi/8 6*pi/8 ];
% obj.ThbRefinement1DML(3,refArea3);

ps = PoissSolv(obj,f);

%obj.plotBas
[Stiffn, rhs, ~,~] = ps.assembleMl();

y = ps.solveSyst(Stiffn,rhs,dBC);

% fplot(g,[obj.levelBas{1}.a obj.levelBas{1}.b],'b')
% hold on;

uh =  ps.generSolThb(y);%
plot(obj.levelBas{1}.plotVector(1:end-1),uh(1:end-1),'r')
hold on;
plot(obj.getAllKnots,0.01,'k*', 'markers',4)
xlabel("x");
ylabel("y");
% errorEst(uh,g,obj);
[H1err, L2err, Enerr] = ps.errCalc(g,y)
errExport(fileName,H1err,L2err,Enerr,obj)

end
 