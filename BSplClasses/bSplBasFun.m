classdef bSplBasFun < bSplBas
     properties
        index = 1;
     end
     methods (Access = public)
        function obj = bSplBasFun(index,basis) % there are still strange things happening
            % standard constructor for subclass of bSplinBasis 
            % index, basis
            % Nooo, overloading not possible!
            assert(rem(index,1) == 0, 'Variable index has to be an integer between 0 and n.');
            % something strange is happening here! Should be: (index < obj.n) 
            assert((0 <= index) && (index < basis.n),'Variable index has to be an integer between 0 and n.');
             if nargin >0
                obj.index = index;
                obj.a = basis.a;
                obj.b = basis.b;
                obj.knotspan = basis.knotspan;
                obj.p = basis.p;
                obj.resol = basis.resol;
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
     
        % other methods
        function C = generOneBasisFun(obj)
                C = zeros(obj.sP,1);
            for j = 1 : obj.sP
                C(j,1) = OneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j));

            end
        end
        function C = generDersOneBasisFun(obj)
                C = zeros(obj.sP,1);
            for j = 1 : obj.sP % inefficient, basis functions are calculated in vain
                temp = DersOneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j),1);
                C(j,1) = temp(:,2);
            end
        end
        function plotOneBasisFun(obj,C)
            plot(obj.plotVector,C(:,1),'y-');
        end
        
        function v = getSupport(obj)
            disp('Support') % correct? Test!
           if( obj.index > obj.p && obj.index < obj.n - obj.p)  
            v= [obj.index: obj.index+obj.p]; % correct
           elseif (obj.index <=obj.p)
                   v = [obj.p: obj.index+obj.p];
           else % all other cases to the right
               v = [obj.index : obj.n-1];
               
           end
           
        end

     end
end
    
