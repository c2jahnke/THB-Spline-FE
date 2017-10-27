function [parameter, thbML] = param()

parameter.a = 0;
parameter.b = 2;
parameter.p = 1;
parameter.N = 5;
parameter.resol = 0.001;
parameter.lvl = 3;
thbML = thbSplBasML(parameter.a,parameter.b,parameter.p,parameter.N,parameter.resol,parameter.lvl);
parameter.f = @(x)  pi.^2/4*sin(pi/2*x);%source function 
parameter.g = @(x) sin(pi/2*x);

end