classdef hbSplBasML < handle
    properties %(GetAccess = public, SetAccess = private)    
        levelBas = [];
        basis0 = [];
        level = 0; %% add variable for the number of active basis functions
        
    end
     properties (Hidden = true, SetAccess = private)
        foo = 1
    end
     methods (Access = public)
        function obj = hbSplBasML(a,b,p,N,resol,level)
            % constructor for class
            % a, b, p, knotspan, resol
            if nargin >0
                assert(level >= 1, "The number of levels has to be a positive integer.");
            for k = 1 : level
                obj.levelBas{1,k} = hbSplBas(a,b,p,N*2^(k-1),resol/2^(k-1));
            end
                obj.basis0 = obj.levelBas{1,1};
                obj.level = level;
            end
              
        end
        

        
        function x = getAllKnots(obj) % do we want the multiple knots as marking points for level change?
            x = []; %% TEST!!!
            for k = 1: obj.level
                if(isempty(obj.levelBas{k}.refArea)) % wrong, what about trivial refinement
                    if( k > 1)
                    x = [x obj.levelBas{k-1}.refArea(1):obj.levelBas{k}.knotspan:obj.levelBas{k-1}.refArea(2)];
                    break;
                    else
                        x = obj.levelBas{k}.knotVector;
                        return;
                    end
                end
                x = [x obj.levelBas{k}.knotVector(obj.levelBas{k}.p):obj.levelBas{k}.knotspan:obj.levelBas{k}.refArea(1)];
            end
            if(obj.levelBas{k-1}.refArea(2) ~= obj.levelBas{k-1}.b)
            for l = 1:k-1
                x = [x obj.levelBas{k-l}.refArea(2):obj.levelBas{k-l}.knotspan:obj.levelBas{k-l}.b];
            end
            end    
        end
        
        function [lvl, BasisFctInd] = getActiveFctIndU(obj,u)
            lvl = [];
            BasisFctInd = [];
            for k = 1 : obj.level
                indU = obj.levelBas{k}.getIndexU(u);
                indSet = (indU-obj.levelBas{k}.p) :indU; % index of act fcts
                actSet = intersect(indSet,obj.levelBas{k}.activeIndex);
                lvl = [lvl, k*ones(1,length(actSet))];
                BasisFctInd = [ BasisFctInd , actSet];
            end
        end
        
        function [lvl,BasisFctInd,basVal] = evalBasisLvl(obj,u)
            basVal = [];
            [lvl, BasisFctInd] = getActiveFctIndU(obj,u);
            for k = 1 : length(BasisFctInd)
                index = BasisFctInd(k);%p,m,U,i,u
                basValues = OneBasisFun(obj.levelBas{lvl(k)}.p,obj.levelBas{lvl(k)}.m,...
                obj.levelBas{lvl(k)}.knotVector,index,u);
                basVal = [basVal basValues];
        
            end
        end
        
        function basVal = evalBasis(obj,u)
             [~,~,basVal] = obj.evalBasisLvl(u);
        end
        
        function [lvl,BasisFctInd,basVal] = evalDersBasisLvl(obj,u)
            basVal = [];
            [lvl, BasisFctInd] = getActiveFctIndU(obj,u);
            for k = 1 : length(BasisFctInd)
                index = BasisFctInd(k);%(p,m,U,i,u,n)
                basValues = DersOneBasisFun(obj.levelBas{lvl(k)}.p,obj.levelBas{lvl(k)}.m,...
                obj.levelBas{lvl(k)}.knotVector,index,u,1);
                basVal = [basVal basValues'];
            end
        end
        function basVal = evalDersBasis(obj,u)
           [~,~,basVal] = obj.evalDersBasisLvl(u);
        end
     end
end