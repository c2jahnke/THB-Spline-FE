classdef hbSplBasML < handle
    properties %(GetAccess = public, SetAccess = private)    
        levelBas = []; % cell array of hbSplBas
        level = 2; % number of levels
        nOF = 0;% number of functions
        nOE = 0;% number of elements
    end
     properties (Hidden = true, SetAccess = private)
        foo = 1
    end
     methods (Access = public)
        function obj = hbSplBasML(a,b,p,N,resol,level)
            %--------------------------------------------------------------
            % constructor for class hbSplBasML
            % Input:
            % a, b, 
            % p: degree
            % N: number of elements at level 0
            % resol: resolution at level 0
            % level: number of levels
            %--------------------------------------------------------------
            if nargin >0
                assert(level >= 1, 'The number of levels has to be a positive integer.');
            for k = 1 : level
                obj.levelBas{1,k} = hbSplBas(a,b,p,N*2^(k-1),resol/2^(k-1)); %resol needed for solution plot!
            end % set index of level 0 to active
                obj.levelBas{1}.activeIndex = [0 : obj.levelBas{1}.n-1];
                obj.level = level;
                obj.nOF = obj.levelBas{1}.n;
                obj.nOE = length(obj.levelBas{1}.getAllKnots) -1;
            end
              
        end
        
        function plotBas(obj)
            % plot all active functions of the basis
            figure
            for jj = 1 : obj.level-1
            obj.levelBas{jj}.plotBasisStruct(obj.levelBas{jj}.generBasisRed)
            hold on
            end
            obj.levelBas{obj.level}.plotBasisStruct(obj.levelBas{obj.level}.generBasisRed)
            xlabel('x');
            ylabel('y');
        end

        function x = getAllKnots(obj) 
            % returns active knots from all levels
            x = []; 
            if(isempty(obj.levelBas{1}.refArea))
                x = obj.levelBas{1}.knotVector;
                return;
            else
                x = [obj.levelBas{1}.a:obj.levelBas{1}.knotspan:obj.levelBas{1}.refArea(1)];
            for k = 2: obj.level
                if(isempty(obj.levelBas{k}.refArea))
                    break;
                end
                x = [x obj.levelBas{k-1}.refArea(1):obj.levelBas{k}.knotspan:obj.levelBas{k}.refArea(1)];
            end
                x = [x obj.levelBas{k-1}.refArea(1):obj.levelBas{k}.knotspan:obj.levelBas{k-1}.refArea(2)];
            for l = k-1:-1:2
                x = [x obj.levelBas{l}.refArea(2):obj.levelBas{l}.knotspan:obj.levelBas{l-1}.refArea(2)];
            end
                x = [x obj.levelBas{1}.refArea(2):obj.levelBas{1}.knotspan:obj.levelBas{1}.b];
            end
                x = unique(x);
        end
        
        function [lvl, BasisFctInd] = getAllActiveFct(obj)
            % returns all active basis functions in correct order

            BasisFctInd = [];
            lvl = [];
            disp('Generalize for more then 2 levels! Test! Debug!')
            for ll = 1 : obj.level-1
            fActInd = obj.levelBas{ll+1}.activeIndex;
            cActInd = obj.levelBas{ll}.activeIndex;
            if(isempty(obj.levelBas{ll}.refArea))
                break
            end
            if(obj.levelBas{ll}.refArea(2) == obj.levelBas{ll}.b)
                BasisFctInd = [cActInd fActInd];
                lvl = [ll*ones(1,length(cActInd)) (ll+1)*ones(1,length(fActInd))];
            elseif(obj.levelBas{ll}.refArea(1) == obj.levelBas{ll}.a)
                BasisFctInd = [fActInd cActInd];
                lvl = [ (ll+1)*ones(1,length(fActInd)) ones(1,length(cActInd))];
            else
            for k = 1 : length(cActInd)-1
                BasisFctInd = [BasisFctInd k-1];
                lvl = [lvl 1];
                if(k < cActInd(k+1))
                    BasisFctInd = [BasisFctInd fActInd];
                    lvl = [lvl 2*ones(1,length(fActInd))];
                    k = (k+1);
                    break
                end
            end
            BasisFctInd = [BasisFctInd cActInd(k:end)];
            lvl = [lvl ones(1,length(cActInd(k:end)))];
            end
            end
        end
        
        
        function [lvl, BasisFctInd] = getActiveFctIndU(obj,u)
            % returns level and function index of active functions at
            % evaluation point u
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
            % evaluate active basis function at point u
            basVal = [];
            [lvl, BasisFctInd] = getActiveFctIndU(obj,u);
            for k = 1 : length(BasisFctInd)
                index = BasisFctInd(k); % p,m,U,i,u
                basValues = OneBasisFun(obj.levelBas{lvl(k)}.p,obj.levelBas{lvl(k)}.m,...
                obj.levelBas{lvl(k)}.knotVector,index,u);
                basVal = [basVal basValues];
            end
        end
        
        function basVal = evalBasis(obj,u)
            % call for evalBasisLvl without level and index
             [~,~,basVal] = obj.evalBasisLvl(u);
        end
        
        function [lvl,BasisFctInd,basVal] = evalDersBasisLvl(obj,u)
            % evaluate active basis function and derivatives at point u
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
            % call for evalDersBasisLvl without level and index
           [~,~,basVal] = obj.evalDersBasisLvl(u);
        end
        
        function HbRefinement1DML(obj,rLevel,refArea)
            % function for  refinement for class hbSplBasML
            %--------------------------------------------------------------
            % function for  refinement for class hbSplBasML
            % Input:
            % rlevel: the coarser level which one desires to refine
            % refArea: the area \Omega^{\ell+1} as a vector of length 2
            %--------------------------------------------------------------
            
            assert(rLevel < obj.level, 'Error: rLevel >= obj.level');
            obj.levelBas{rLevel}.refArea = refArea;
            
            HbRef1D(obj,rLevel)
            
            obj.nOF = 0;
            for k = 0: rLevel
                obj.nOF = obj.nOF + length(obj.levelBas{k+1}.activeIndex);
            end
            obj.nOE = length(obj.getAllKnots)-1;
        end
     end
end