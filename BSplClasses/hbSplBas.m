classdef hbSplBas < bSplBas
     properties
        %activeKnots = 0; % for basis refinement
        activeIndex = 0;% function index
        refArea = []; % area for refinement
        charM = [];% add structure for characteristic matrix
     end
     methods (Access = public)
         function obj = hbSplBas(a,b,p,N,resol)
            %--------------------------------------------------------------
            % constructor for class hbSplBas
            % Input:
            % a, b, 
            % p: degree
            % N: number of elements at level 0
            % resol: resolution at level 0
            %--------------------------------------------------------------
          
            if nargin >0
                obj.a = a;
                obj.b = b;
                obj.N = N;
                obj.p = p;
                obj.resol = resol;
                
            end
                obj.knotspan = (obj.b-obj.a)/obj.N;
                obj.knotVector = ConstrKnotVector(obj.a,obj.b,obj.knotspan,obj.p);
                obj.plotVector = [obj.knotVector(1):obj.resol:obj.knotVector(end)];
                obj.sP = size(obj.plotVector,2);
                obj.m = size(obj.knotVector,2);
                obj.n = obj.m - obj.p - 1;
                %obj.activeKnots = [];%obj.getAllKnots;
                obj.activeIndex = [];%[0 : obj.n-1]; % function index
                obj.charM = sparse(obj.n, obj.n);
            function knotVector = ConstrKnotVector(a,b,h1,p)
                % simple constructor for knot vector
                temp = (b-a)/h1;
                assert( rem(temp,1) == 0,'Use values a,b,h1 such that b-a/h1 is an integer.')
                m = temp + 2*p+1;
                knotVector = zeros(1,m);
                knotVector(1:p) = a;
                knotVector(p+1:m-p) = a:h1:b;
                knotVector(m-p+1:end) = b;
            end
         end
         function C = generBasisRed(obj) % überladen
            %% calculate only active basis functions!
            C = zeros(obj.sP,obj.n);
            tableSpan = lookUpTableSpan(obj);
            for i = 1 : obj.sP
                startX = tableSpan(i) - tableSpan(1) +1;
                bas_temp = BasisFuns(tableSpan(i),obj.plotVector(i),obj.p,obj.knotVector);
                for j = obj.activeIndex % plot only active basis functions
                    if (j-tableSpan(i) + obj.p) < obj.p+1 && tableSpan(i) - j  < obj.p+1
                        tmp = mod(j+1 - startX, (obj.p+1))+1;  
                        C(i,j+1) = bas_temp(1,tmp);
                    end
                end
            end
            
         end
         function setCharM(obj) %characteristic matrix
             % obj.charM = zeros(obj.n, obj.n);
         for k = obj.activeIndex 
                obj.charM(k+1,k+1) = 1;
         end
         end
     end
end