classdef thbSplBasML < hbSplBasML
    properties %(GetAccess = public, SetAccess = private)    
    end

     methods (Access = public)
        function obj = thbSplBasML(a,b,p,N,resol,level)
            %--------------------------------------------------------------
            % constructor for class thbSplBasML
            % Input:
            % a, b, 
            % p: degree
            % N: number of elements at level 0
            % resol: resolution at level 0
            % level: total number of levels
            %--------------------------------------------------------------
            if nargin >0
                assert(level >= 1, 'The number of levels has to be a positive integer.');
            for k = 1 : level
                obj.levelBas{1,k} = thbSplBas(a,b,p,N*2^(k-1),resol/2^(k-1)); %needed for solution plot!
            end% set index of level 0 to active
                obj.levelBas{1}.activeIndex = [0 : obj.levelBas{1}.n-1];
                %obj.levelBas{1}.activeKnots = obj.levelBas{1}.getAllKnots;
                %obj.basis0 = obj.levelBas{1,1};
                obj.level = level;
                obj.nOF = obj.levelBas{1}.n;
                obj.nOE = length(obj.levelBas{1}.getAllKnots) -1;
            end
              
        end
        
        function plotBas(obj)
            figure
            for jj = 1 : obj.level-1
            obj.levelBas{jj}.plotBasisStruct(obj.levelBas{jj}.generBasisRed(obj.levelBas{jj+1}))
            hold on
            end
            obj.levelBas{obj.level}.plotBasisStruct(obj.levelBas{obj.level}.generBasisRed(obj.levelBas{obj.level}))
  
        end
 
        
        function [lvl,BasisFctInd,basVal] = evalBasisLvl(obj,u)
            %% get rid of index tmp, tmp2 if possible
            basVal = [];
            [lvl, BasisFctInd] = getActiveFctIndU(obj,u);
            for k = 1 : length(BasisFctInd)
                basValues = 0;
                 % check if index is in truncInd
                if(ismember(BasisFctInd(k),obj.levelBas{lvl(k)}.truncIndex))
                    tmp = obj.levelBas{lvl(k)}.trunc(:,BasisFctInd(k)+1);
                    tmp2 = find(tmp);
                    for l = 1 : length(tmp2)
                        basValues =  basValues + tmp(tmp2(l))*OneBasisFun(obj.levelBas{lvl(k)+1}.p,...
                            obj.levelBas{lvl(k)+1}.m,...
                        obj.levelBas{lvl(k)+1}.knotVector,tmp2(l)-1,u);
                    end
                else
                    basValues = OneBasisFun(obj.levelBas{lvl(k)}.p,obj.levelBas{lvl(k)}.m,...
                    obj.levelBas{lvl(k)}.knotVector,BasisFctInd(k),u);
                end
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
                basValues = 0;
                if(ismember(BasisFctInd(k),obj.levelBas{lvl(k)}.truncIndex) )
                    tmp = obj.levelBas{lvl(k)}.trunc(:,BasisFctInd(k)+1);
                    tmp2 = find(tmp);
                    for l = 1 : length(tmp2)
                        
                        basValues =  basValues + tmp(tmp2(l))*DersOneBasisFun(obj.levelBas{lvl(k)+1}.p,...
                            obj.levelBas{lvl(k)+1}.m,...
                        obj.levelBas{lvl(k)+1}.knotVector,tmp2(l)-1,u,1);
                    end
                else
                    basValues = DersOneBasisFun(obj.levelBas{lvl(k)}.p,obj.levelBas{lvl(k)}.m,...
                    obj.levelBas{lvl(k)}.knotVector,BasisFctInd(k),u,1);
                end
                basVal = [basVal basValues'];
            end
        end
        function basVal = evalDersBasis(obj,u)
           [~,~,basVal] = obj.evalDersBasisLvl(u);
        end
        
        function ThbRefinement1DML(obj,rLevel,refArea)
            % function for THB refinement, including truncation
            %--------------------------------------------------------------
            % function for THB refinement, including truncation
            % rLevel: level-1 to be refined (starts with 1)
            % refArea: domain Omega^{\ell+1} for refinement, e.g.
            % [2 8] results in refinement between 2 and 8
            %--------------------------------------------------------------
            % Test: hb refinement!
            assert(rLevel < obj.level, 'Error: rLevel >= obj.level');

            obj.levelBas{rLevel}.refArea = refArea;
            HbRef1D(obj,rLevel); 
            obj.nOF = 0;
            for k = 0: rLevel
                obj.nOF = obj.nOF + length(obj.levelBas{k+1}.activeIndex);
            end
            obj.nOE = length(obj.getAllKnots)-1;
            ThbRef1D(obj,rLevel);            
            % check partition of unity, important for debugging
            for k = 1 : obj.levelBas{rLevel}.sP
                eVal =sum(obj.evalBasis(obj.levelBas{rLevel}.plotVector(k)));
               if(abs(eVal-1) > 10e-2 && obj.levelBas{rLevel}.plotVector(k) ~=...
                       obj.levelBas{rLevel}.b)          
                    disp('Warning, sum(obj.evalBasis(obj.levelBas{rLevel}.plotVector(k))) = ');
                    disp(eVal);
                    disp(' for k = ');
                    disp(k);
                    disp('Partition of unity is violated.');
                end
            end
        end

     end
end