classdef thbSplBasML < hbSplBasML
    properties %(GetAccess = public, SetAccess = private)    
    end

     methods (Access = public)
        function obj = thbSplBasML(a,b,p,N,resol,level)
            % constructor for class
            % a, b, p, knotspan, resol
            if nargin >0
                assert(level >= 1, 'The number of levels has to be a positive integer.');
            for k = 1 : level
                obj.levelBas{1,k} = thbSplBas(a,b,p,N*2^(k-1),resol/2^(k-1)); %needed for solution plot!
            end% set index of level 0 to active
                obj.levelBas{1}.activeIndex = [0 : obj.levelBas{1}.n-1];
                obj.levelBas{1}.activeKnots = obj.levelBas{1}.getAllKnots;
                obj.basis0 = obj.levelBas{1,1};
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
        
%        function x = getAllKnots(obj) % do we want the multiple knots as marking points for level change?
%             x = []; %% TEST!!!
%             for k = 1: obj.level
%                 if(isempty(obj.levelBas{k}.refArea)) % wrong, what about trivial refinement
%                     if( k > 1)
%                     x = [x obj.levelBas{k-1}.refArea(1):obj.levelBas{k}.knotspan:obj.levelBas{k-1}.refArea(2)];
%                     break;
%                     else
%                         x = obj.levelBas{k}.knotVector;
%                         return;
%                     end
%                 end
%                 x = [x obj.levelBas{k}.knotVector(obj.levelBas{k}.p):obj.levelBas{k}.knotspan:obj.levelBas{k}.refArea(1)];
%                 
%             end
%             if(obj.levelBas{k-1}.refArea(2) ~= obj.levelBas{k-1}.b)
%             for l = 1:k-1
%                 x = [x obj.levelBas{k-l}.refArea(2):obj.levelBas{k-l}.knotspan:obj.levelBas{k-l}.b];
%             end
%             end    
%             x = unique(x);
%        end

%        function [lvl, BasisFctInd] = getAllActiveFct(obj)
            % returns all active basis functions in correct order
            % just for 2 levels, 
            % TEST!
            % generalize for more levels
%             BasisFctInd = [];
%             lvl = [];
%             fActInd = obj.levelBas{2}.activeIndex;
%             cActInd = obj.levelBas{1}.activeIndex;
%             if(obj.levelBas{1}.refArea(2) == obj.levelBas{1}.b)
%                 BasisFctInd = [cActInd fActInd];
%                 lvl = [ones(1,length(cActInd)) 2*ones(1,length(fActInd))];
%             elseif(obj.levelBas{1}.refArea(1) == obj.levelBas{1}.a)
%                 BasisFctInd = [fActInd cActInd];
%                 lvl = [ 2*ones(1,length(fActInd)) ones(1,length(cActInd))];
%             else
%             for k = 1 : length(cActInd)-1
%                 BasisFctInd = [BasisFctInd k-1];
%                 lvl = [lvl 1];
%                 if(k < cActInd(k+1))
%                     BasisFctInd = [BasisFctInd fActInd];
%                     lvl = [lvl 2*ones(1,length(fActInd))];
%                     k = (k+1);
%                     break
%                 end
%             end
%             BasisFctInd = [BasisFctInd cActInd(k:end)];
%             lvl = [lvl ones(1,length(cActInd(k:end)))];
%             end
%        end
        
        
%        function [lvl, BasisFctInd] = getActiveFctIndU(obj,u)
%             lvl = [];
%             BasisFctInd = [];
%             for k = 1 : obj.level
%                 indU = obj.levelBas{k}.getIndexU(u);
%                 indSet = (indU-obj.levelBas{k}.p) :indU; % index of act fcts
%                 actSet = intersect(indSet,obj.levelBas{k}.activeIndex);
%                 lvl = [lvl, k*ones(1,length(actSet))];
%                 BasisFctInd = [ BasisFctInd , actSet];
%             end
%        end
        
        function [lvl,BasisFctInd,basVal] = evalBasisLvl(obj,u)
            %% get rid of index tmp, tmp2 if possible
            basVal = [];
            [lvl, BasisFctInd] = getActiveFctIndU(obj,u);
            for k = 1 : length(BasisFctInd)
                basValues = 0;
                 % check if index is in truncInd
                if(ismember(BasisFctInd(k),obj.levelBas{lvl(k)}.truncIndex))
                    %disp('Truncation necessary!')
                    tmp = obj.levelBas{lvl(k)}.trunc(:,BasisFctInd(k)+1);
                    tmp2 = find(tmp);
                    for l = 1 : length(tmp2)
                        % bug
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
                    %disp('Truncation necessary!')
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
        
        function [THB0, THB1, trunq, q, trP] = ThbRefinement1DML(obj,rLevel,refArea,f)
            % Test: hb refinement!
            % function f not needed, as Points are not needed
            assert(rLevel < obj.level, 'Error: rLevel >= obj.level');
            Points = zeros(obj.levelBas{rLevel}.n,2); % not needed
            
            obj.levelBas{rLevel}.refArea = refArea;
            [~, ~, ~, Qw] = HbRef1D(obj,rLevel, Points); 
            obj.nOF = 0;
            for k = 0: rLevel
                obj.nOF = obj.nOF + length(obj.levelBas{k+1}.activeIndex);
            end
            obj.nOE = length(obj.getAllKnots)-1;
            [THB0, THB1, trunq, q, trP] =  ThbRef1D(obj,rLevel, Qw);            % check partition of unity
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