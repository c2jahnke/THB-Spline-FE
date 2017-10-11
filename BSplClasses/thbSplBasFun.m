classdef thbSplBasFun < thbSplBasML
     properties
         index = [];
         truncIndex = [];
         trunc = [];
         a = 0; % left interval
         b = 10; % right interval
         N = 5;
         knotspan = 2;
         p = 2; % spline degree
         knotVector = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
         resol = 0.1; % resolution of plot vector
         % properties initialized in constructor
         plotVector = 0;
         sP = 0; % number of evaluation points in plotVector
         m = 0; % number of knots in knotvector
         n = 0; % number of basis functions, n = m - p -1
         lev = 0;
     end
     methods (Access = public)
        function obj = thbSplBasFun(index,thbML,lev) 
            % standard constructor for subclass of bSplBas 
            % index, basis
            assert(rem(index,1) == 0, 'Variable index has to be an integer between 0 and n.');
            % something strange is happening here! Should be: (index < obj.n) 
            basis = thbML.levelBas{lev}; % thbML.nOF
            assert((0 <= index) && (index < basis.n),'Variable index has to be an integer between 0 and n.');
             if nargin >0
                obj.index = index;
                obj.lev = lev;
                obj.a = basis.a;
                obj.b = basis.b;
                obj.knotspan = basis.knotspan;
                obj.p = basis.p;
                obj.resol = basis.resol;
                obj.truncIndex = basis.truncIndex;
                obj.trunc = basis.trunc;
                obj.levelBas = thbML.levelBas;
                obj.basis0 = thbML.basis0;
            end
                obj.knotVector = ConstrKnotVector(obj.a,obj.b,obj.knotspan,obj.p);
                obj.plotVector = [obj.knotVector(1):obj.resol:obj.knotVector(end)];
                obj.sP = size(obj.plotVector,2);
                obj.m = size(obj.knotVector,2);
                obj.n = obj.m - obj.p - 1;
                 
                
            function knotVector = ConstrKnotVector(a,b,h1,p)
                % simple constructor for knot vector
                temp = (b-a)/h1;
                assert( isinteger(temp) == 0,'Use values a,b,h1 such that b-a/h1 is an integer.')
                m = temp + 2*p+1;
                knotVector = zeros(1,m);
                knotVector(1:p) = a;
                knotVector(p+1:m-p) = a:h1:b;
                knotVector(m-p+1:end) = b;
            end
        end
     
        % readjust this method for THB assembly
        function C = generOneBasisFun(obj)
            
            if(ismember(obj.index,obj.truncIndex))
                tmp = obj.trunc(:,obj.index +1);
                tmp2 = find(tmp);
                C = zeros(obj.levelBas{obj.lev+1}.sP,1);
                for k = tmp2(1):tmp2(end)
                   for j = 1 : obj.levelBas{obj.lev+1}.sP
                       C(j,1) = C(j,1) + tmp(k)*OneBasisFun(obj.levelBas{obj.lev+1}.p,...
                           obj.levelBas{obj.lev+1}.m,obj.levelBas{obj.lev+1}.knotVector,...
                          k-1,obj.levelBas{obj.lev+1}.plotVector(j));
                   end
                end
                C = C(1:2:end);
            else
                C = zeros(obj.sP,1);
            for j = 1 : obj.sP
                C(j,1) = OneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j));
            end
            end
        end
%         function C = generDersOneBasisFun(obj)
%                 C = zeros(obj.sP,1);
%             for j = 1 : obj.sP % inefficient, basis functions are calculated in vain
%                 temp = DersOneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j),1);
%                 C(j,1) = temp(:,2);
%             end
%         end
%         function plotOneBasisFun(obj,C)
%             plot(obj.plotVector,C(:,1),'y-');
%         end
%         
%         function v = getSupport(obj)
%             disp('Support') % correct? Test!
%            if( obj.index > obj.p && obj.index < obj.n - obj.p)  
%             v= [obj.index: obj.index+obj.p]; % correct
%            elseif (obj.index <=obj.p)
%                    v = [obj.p: obj.index+obj.p];
%            else % all other cases to the right
%                v = [obj.index : obj.n-1];
%                
%            end
%            
%         end

     end
end
    
