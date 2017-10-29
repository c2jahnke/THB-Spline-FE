    classdef bSplBas2D < handle
    %% Class for BSplBas
    % so far contains: Basis generation
    % tablespan
    % derivative generation
    % Basis plotting
    % getIndices
    % TODO: write evalBasis, evalDersBasis, getIndexU
    properties %(GetAccess = public, SetAccess = private)    
        aU = 0; % left interval
        bU = 10; % right interval
        NU = 5;
        knotspanU = 2;
        pU = 2; % spline degree
        knotVectorU = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
        resolU = 0.1; % resolution of plot vector
        % properties initialized in constructor
        plotVectorU = 0;
        sPU = 0; % number of evaluation points in plotVector
        mU = 0; % number of knots in knotvector
        nU = 0; % number of basis functions, n = m - p -1
        
        aV = 0; % left interval
        bV = 10; % right interval
        NV = 5;
        knotspanV = 2;
        pV = 2; % spline degree
        knotVectorV = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
        resolV = 0.1; % resolution of plot vector
        % properties initialized in constructor
        plotVectorV = 0;
        sPV = 0; % number of evaluation points in plotVector
        mV = 0; % number of knots in knotvector
        nV = 0; % number of basis functions, n = m - p -1
    end
     properties (Hidden = true, SetAccess = private)

    end
     methods (Access = public)
        function obj = bSplBas2D(bSplBasU,bSplBasV)
            %--------------------------------------------------------------
            % constructor for class bSplBas2D
            % Input:
            % bSplBasU: in u direction
            % bSplBasV: in v direction
            %--------------------------------------------------------------
            if nargin >0
                obj.aU = bSplBasU.a;
                obj.bU = bSplBasU.b;
                obj.NU = bSplBasU.N;
                obj.pU = bSplBasU.p;
                obj.resolU = bSplBasU.resol;
                
                                obj.aV = bSplBasV.a;
                obj.bV = bSplBasV.b;
                obj.NV = bSplBasV.N;
                obj.pV = bSplBasV.p;
                obj.resolV = bSplBasV.resol;
                
            end
                
                
                if((obj.bU==obj.aU))
                    obj.knotspanU = 0;
                else
                obj.knotspanU = (obj.bU-obj.aU)/obj.NU;
                obj.knotVectorU = ConstrKnotVector(obj.aU,obj.bU,obj.knotspanU,obj.pU);
                obj.plotVectorU = [obj.knotVectorU(1):obj.resolU:obj.knotVectorU(end)];
                obj.sPU = size(obj.plotVectorU,2);
                end
                obj.mU = size(obj.knotVectorU,2);
                obj.nU = obj.mU - obj.pU - 1;

                if((obj.bV==obj.aV))
                obj.knotspanV = 0;
                else
                obj.knotspanV = (obj.bV-obj.aV)/obj.NV;
                obj.knotVectorV = ConstrKnotVector(obj.aV,obj.bV,obj.knotspanV,obj.pV);
                obj.plotVectorV = [obj.knotVectorV(1):obj.resolV:obj.knotVectorV(end)];
                obj.sPV = size(obj.plotVectorV,2);
                end
                obj.mV = size(obj.knotVectorV,2);
                obj.nV = obj.mV - obj.pV - 1;
                
            function knotVector = ConstrKnotVector(a,b,h1,p)
                % simple constructor for knot vector
                temp = (b-a)/h1; % = N
                assert( rem(temp,1) == 0,'Use values a,b,h1 such that b-a/h1 is an integer.')
                m = temp + 2*p+1;
                knotVector = zeros(1,m);
                knotVector(1:p) = a;
                knotVector(p+1:m-p) = a:h1:b;
                knotVector(m-p+1:end) = b;
            end
            

              
        end
     
         function  plotOne2Dbasis(obj,iU,iV)
            [X,Y] = meshgrid(obj.knotVectorU(1):obj.resolU:obj.knotVectorU(end),obj.knotVectorU(1):obj.resolU:obj.knotVectorU(end));

            CU = zeros(obj.sPU,1);
            for kk = 1 : obj.sPU
            CU(kk,1) = OneBasisFun(obj.pU,obj.mU,obj.knotVectorU,iU,obj.plotVectorU(kk));

            end
            CV = zeros(obj.sPV,1);
            for kk = 1 : obj.sPV
            CV(kk,1) = OneBasisFun(obj.pV,obj.mV,obj.knotVectorV,iV,obj.plotVectorV(kk));

            end
            Z = CU(:,1)*CV(:,1)'; % tensor product structure
            surf(X,Y,Z)
            %hold off
        end
     end
     
     
%         function obj = bSplBasKV(knotVector)
%             obj.a = knotVector(1);
%             obj.b = knotVector(end)
%             obj.knotVector = knotVector;
%             
%         end
%         
%         function C = generBasis(obj)
% %         % 1.
% %         C = zeros(obj.sP,obj.n);
% %         for ll = 1:obj.n
% %             for ks = 1 : obj.sP
% %                 C(ks,ll) = OneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks));
% %             end
% %         end
%         % 2. 
%         %The second way is a lot faster!
%             C = zeros(obj.sP,obj.n);
%             tableSpan = lookUpTableSpan(obj);
%             for i = 1 : obj.sP
%                 startX = tableSpan(i) - tableSpan(1) +1;
%                 bas_temp = BasisFuns(tableSpan(i),obj.plotVector(i),obj.p,obj.knotVector);
%                 for j = 1 : obj.n 
%                     if (j-tableSpan(i) + obj.p -1 ) < obj.p+1 && tableSpan(i) +1 - j  < obj.p+1
%                         tmp = mod(j - startX, (obj.p+1))+1;  % reorder the functions
%                         C(i,j) = bas_temp(1,tmp);
%                     end
%                 end
%             end
%             
%             
%         end
%         
%         
%         function plotBasisStruct(obj,C)
%         %assert(isequal(size(C),[obj.sP,obj.n]),'Warning: Basis structure may not be plottable.');
%         for ll = 1:obj.n
%             plot(obj.plotVector,C(:,ll),'LineWidth',1.2);%'color',[0.7,0.7,0.7]);%'color',[0,0,0]);%
%             hold all
%         end
%         xlabel('x');
%         ylabel('y');
%         set(gca,'ytick',[])
%         hold off
%         end
%         % TODO: Test!
%         function x = evalBasis(obj, u)
%             assert( (u >= obj.a) && (u <= obj.b),'u is out of range.') 
%             % returns all active Basis functions
%             x = BasisFuns(obj.getIndexU(u),u,obj.p,obj.knotVector);
%         end
%         
%         function x = evalDersBasis(obj, u)
%            assert( (u >= obj.a) && (u <= obj.b),'u is out of range.') 
%            % returns all active Basis function derivatives
%            x = DersBasisFuns(obj.getIndexU(u),u,obj.p,obj.knotVector,1);
%         end
%         % TODO:  Something is wrong here!
%         function cnt = getIndexU(obj,u) % returns index as in lookUpTableSpan
%             assert( (obj.a <= u) & (u <= obj.b),'u out of range of knotVector.');
%             if(u == obj.a) % is this correct?
%                 cnt = obj.p;
%                	return 
%             end
%             if(u == obj.b)
%                 cnt = obj.m-obj.p-1;%% work around, think of better one!
%                 return;
%             end
%             left = obj.a;
%             cnt = 0;
%             while(u >= left) % changed to >= to get a new index for first 
%                 cnt = cnt + 1; % numerical rounding error at u >= left!
%                 left = obj.knotVector(cnt+2) - 2*eps; % to remain in loop subtract half the knotspan!
%             end
%         end
%         
%         function allKnots = getAllKnots(obj)
%             allKnots = obj.knotVector(obj.p+1: end-obj.p);
%         end
%  
%         %%%%% Isogat algorithm
%         function tableSpan = lookUpTableSpan(obj)
%         %--------------------------------------------------------------
%         % Input:
%         % knotVector : knot vector
%         % plotVector : plot vector
%         % m          : size of knot vector
%         % sP          : size of plot vector
%         % Output
%         % tableSpan : row vector of same dimension as plotVector
%         %--------------------------------------------------------------
%         tableSpan = zeros(1,obj.sP);
%         plotIndex = 1;
%         leftIndex = 1; % left boundary of knot span
%         for knotIndex = 2:obj.m
%             rightValue = obj.knotVector(knotIndex);
%             if(obj.knotVector(leftIndex) == rightValue)
%                 continue
%             end
%             leftIndex = knotIndex -1;
%             while(obj.plotVector(plotIndex) < rightValue)
%                 tableSpan(plotIndex) = leftIndex;
%                 plotIndex = plotIndex +1;
%             end
%         end
%         tableSpan(plotIndex) = leftIndex -1;
%         tableSpan = tableSpan -ones(size(tableSpan));
%         end
%         
%         function v = getIndices(obj)
%             % see also lookUpTableSpan
%             v.knotVector = obj.knotVector;
%             v.allKnotIndex = [ 0 : obj.m-1];
%             %v.activeKnots = [obj.p : obj.n-1];
%             v.basisFunctionIndex = [0:obj.n-1];
%         end
%         
%         function C = gener1Deriv(obj)
%         %--------------------------------------------------------------
%         % Warning: partly own idea, test carefully!
%         % Input:
%         % knotVector : knot vector
%         % plotVector : plot vector
%         % tableSpan  : lookUpTableSpan
%         % p          : polynomial degree
%         % sP         : size of plot vector
%         %Output
%         % tableSpan : row vector of same dimension as plotVector
%         %-------------------------------------------------------------
%         % 1.
% %         C = zeros(obj.sP,obj.n);
% %         for ll = 1:obj.n
% %             for ks = 1 : obj.sP
% % 
% %                 temp = DersOneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks),1);
% %                 C(ks,ll) = temp(2); %OneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks));
% %             end
% %         end
%         % 2.
%         % For derivatives method 2 is much faster!
%             C = zeros(obj.sP,obj.n);
%             tableSpan = lookUpTableSpan(obj);
%             % use getPlotMatrix line 1040 ff
%             % use startX to switch index when entering new span index
%             for i = 1 : obj.sP
%                 startX = tableSpan(i) - tableSpan(1) +1;
%                 ders = DersBasisFuns(tableSpan(i),obj.plotVector(i),obj.p,obj.knotVector,1);
%                 for j = 1 : obj.n % change to n
%                     if (j-tableSpan(i) + obj.p -1 ) < obj.p+1 && tableSpan(i) +1 - j  < obj.p+1
%                         tmp = mod(j - startX, (obj.p+1))+1;  % temporary solution
%                         C(i,j) = ders(2,tmp);
%                     end
%                 end
%             end
% 
%         end
% 
% 
%         
%      end
end
    
    
        
        

    
   
