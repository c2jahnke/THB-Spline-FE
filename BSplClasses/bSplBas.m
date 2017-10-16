    classdef bSplBas < handle
    %% Class for BSplBas
    % so far contains: Basis generation
    % tablespan
    % derivative generation
    % Basis plotting
    % getIndices
    % TODO: write evalBasis, evalDersBasis, getIndexU
    properties %(GetAccess = public, SetAccess = private)    
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
    end
     properties (Hidden = true, SetAccess = private)

    end
     methods (Access = public)
        function obj = bSplBas(a,b,p,N,resol)
            % constructor for class
            % a, b, p, knotspan, resol
            if nargin >0
                obj.a = a;
                obj.b = b;
                obj.N = N;
                obj.p = p;
                obj.resol = resol;
                
            end
                if((obj.b==obj.a))
                    obj.knotspan = 0;
                else
                obj.knotspan = (obj.b-obj.a)/obj.N;
                obj.knotVector = ConstrKnotVector(obj.a,obj.b,obj.knotspan,obj.p);
                obj.plotVector = [obj.knotVector(1):obj.resol:obj.knotVector(end)];
                obj.sP = size(obj.plotVector,2);
                end
                obj.m = size(obj.knotVector,2);
                obj.n = obj.m - obj.p - 1;
                
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
        
        function C = generBasis(obj)
%         % 1.
%         C = zeros(obj.sP,obj.n);
%         for ll = 1:obj.n
%             for ks = 1 : obj.sP
%                 C(ks,ll) = OneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks));
%             end
%         end
        % 2. 
        %The second way is a lot faster!
            C = zeros(obj.sP,obj.n);
            tableSpan = lookUpTableSpan(obj);
            for i = 1 : obj.sP
                startX = tableSpan(i) - tableSpan(1) +1;
                bas_temp = BasisFuns(tableSpan(i),obj.plotVector(i),obj.p,obj.knotVector);
                for j = 1 : obj.n 
                    if (j-tableSpan(i) + obj.p -1 ) < obj.p+1 && tableSpan(i) +1 - j  < obj.p+1
                        tmp = mod(j - startX, (obj.p+1))+1;  % reorder the functions
                        C(i,j) = bas_temp(1,tmp);
                    end
                end
            end
            
            
        end
        
        
        function plotBasisStruct(obj,C)
        %assert(isequal(size(C),[obj.sP,obj.n]),'Warning: Basis structure may not be plottable.');
        for ll = 1:obj.n
            plot(obj.plotVector,C(:,ll));%,'-k');
            hold all
        end
        hold off
        end
        % TODO: Test!
        function x = evalBasis(obj, u)
            assert( (u >= obj.a) && (u <= obj.b),'u is out of range.') 
            % returns all active Basis functions
            x = BasisFuns(obj.getIndexU(u),u,obj.p,obj.knotVector);
        end
        
        function x = evalDersBasis(obj, u)
           assert( (u >= obj.a) && (u <= obj.b),'u is out of range.') 
           % returns all active Basis function derivatives
           x = DersBasisFuns(obj.getIndexU(u),u,obj.p,obj.knotVector,1);
        end
        % TODO:  Something is wrong here!
        function cnt = getIndexU(obj,u) % returns index as in lookUpTableSpan
            assert( (obj.a <= u) & (u <= obj.b),'u out of range of knotVector.');
            if(u == obj.a) % is this correct?
                cnt = obj.p;
               	return 
            end
            if(u == obj.b)
                cnt = obj.m-obj.p-1;%% work around, think of better one!
                return;
            end
            left = obj.a;
            cnt = 0;
            while(u >= left) % changed to >= to get a new index for first 
                cnt = cnt + 1;
                left = obj.knotVector(cnt+2);
            end
        end
        
        function allKnots = getAllKnots(obj)
            allKnots = obj.knotVector(obj.p+1: end-obj.p);
        end
 
        %%%%% Isogat algorithm
        function tableSpan = lookUpTableSpan(obj)
        %--------------------------------------------------------------
        % Input:
        % knotVector : knot vector
        % plotVector : plot vector
        % m          : size of knot vector
        % sP          : size of plot vector
        % Output
        % tableSpan : row vector of same dimension as plotVector
        %--------------------------------------------------------------
        tableSpan = zeros(1,obj.sP);
        plotIndex = 1;
        leftIndex = 1; % left boundary of knot span
        for knotIndex = 2:obj.m
            rightValue = obj.knotVector(knotIndex);
            if(obj.knotVector(leftIndex) == rightValue)
                continue
            end
            leftIndex = knotIndex -1;
            while(obj.plotVector(plotIndex) < rightValue)
                tableSpan(plotIndex) = leftIndex;
                plotIndex = plotIndex +1;
            end
        end
        tableSpan(plotIndex) = leftIndex -1;
        tableSpan = tableSpan -ones(size(tableSpan));
        end
        
        function v = getIndices(obj)
            % see also lookUpTableSpan
            v.knotVector = obj.knotVector;
            v.allKnotIndex = [ 0 : obj.m-1];
            %v.activeKnots = [obj.p : obj.n-1];
            v.basisFunctionIndex = [0:obj.n-1];
        end
        
        function C = gener1Deriv(obj)
        %--------------------------------------------------------------
        % Warning: partly own idea, test carefully!
        % Input:
        % knotVector : knot vector
        % plotVector : plot vector
        % tableSpan  : lookUpTableSpan
        % p          : polynomial degree
        % sP         : size of plot vector
        %Output
        % tableSpan : row vector of same dimension as plotVector
        %-------------------------------------------------------------
        % 1.
%         C = zeros(obj.sP,obj.n);
%         for ll = 1:obj.n
%             for ks = 1 : obj.sP
% 
%                 temp = DersOneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks),1);
%                 C(ks,ll) = temp(2); %OneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks));
%             end
%         end
        % 2.
        % For derivatives method 2 is much faster!
            C = zeros(obj.sP,obj.n);
            tableSpan = lookUpTableSpan(obj);
            % use getPlotMatrix line 1040 ff
            % use startX to switch index when entering new span index
            for i = 1 : obj.sP
                startX = tableSpan(i) - tableSpan(1) +1;
                ders = DersBasisFuns(tableSpan(i),obj.plotVector(i),obj.p,obj.knotVector,1);
                for j = 1 : obj.n % change to n
                    if (j-tableSpan(i) + obj.p -1 ) < obj.p+1 && tableSpan(i) +1 - j  < obj.p+1
                        tmp = mod(j - startX, (obj.p+1))+1;  % temporary solution
                        C(i,j) = ders(2,tmp);
                    end
                end
            end

        end


        
     end
end
    
    
        
        

    
   
