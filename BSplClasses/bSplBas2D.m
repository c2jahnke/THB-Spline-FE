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
            [X,Y] = meshgrid(obj.knotVectorV(1):obj.resolV:obj.knotVectorV(end),obj.knotVectorU(1):obj.resolU:obj.knotVectorU(end));
            CU = zeros(obj.sPU,1);
            for kk = 1 : obj.sPU
            CU(kk,1) = OneBasisFun(obj.pU,obj.mU,obj.knotVectorU,iU,obj.plotVectorU(kk));

            end
            CV = zeros(obj.sPV,1);
            for kk = 1 : obj.sPV
            CV(kk,1) = OneBasisFun(obj.pV,obj.mV,obj.knotVectorV,iV,obj.plotVectorV(kk));

            end
            Z = CU(:,1)*CV(:,1)'; % tensor product structure
            s = surf(X,Y,Z);
            s.EdgeColor = 'none';
            
            %hold off
         end
        
         function plot2Dbasis(obj) % TEST for bugs!
            for k = 0 : obj.nU
                for l = 0 : obj.nV
                    obj.plotOne2Dbasis(k,l)
                    hold on;
                end
            end
            xlabel('x');
            ylabel('y');
         end
         function x = evalBasis(obj, UV)
             % evalute basis at point UV = (u,v)
            assert( (UV(1) >= obj.aU) && (UV(1) <= obj.bU) && (UV(2) >= obj.aV) && (UV(2) <= obj.bV) ,'u is out of range.') 
            % returns all active Basis functions
            x1 = BasisFuns(obj.getIndexU(UV(1)),UV(1),obj.pU,obj.knotVectorU)
            x2 = BasisFuns(obj.getIndexV(UV(2)),UV(2),obj.pV,obj.knotVectorV);
            x = x1.*x2
         end
         function cnt = getIndexU(obj,u) % returns index as in lookUpTableSpan
            assert( (obj.aU <= u) & (u <= obj.bU),'u out of range of knotVector.');
            if(u == obj.aU) % is this correct?
                cnt = obj.pU;
               	return 
            end
            if(u == obj.bU)
                cnt = obj.mU-obj.pU-1;%% work around, think of better one!
                return;
            end
            left = obj.aU;
            cnt = 0;
            while(u >= left) % changed to >= to get a new index for first 
                cnt = cnt + 1; % numerical rounding error at u >= left!
                left = obj.knotVectorU(cnt+2) - 2*eps; % to remain in loop subtract half the knotspan!
            end
         end
         function cnt = getIndexV(obj,v) % returns index as in lookUpTableSpan
            assert( (obj.aV <= v) & (v <= obj.bV),'u out of range of knotVector.');
            if(v == obj.aV) % is this correct?
                cnt = obj.pV;
               	return 
            end
            if(v == obj.bV)
                cnt = obj.mV-obj.pV-1;%% work around, think of better one!
                return;
            end
            left = obj.aV;
            cnt = 0;
            while(v >= left) % changed to >= to get a new index for first 
                cnt = cnt + 1; % numerical rounding error at u >= left!
                left = obj.knotVectorV(cnt+2) - 2*eps; % to remain in loop subtract half the knotspan!
            end
        end
     end
     
     
    
end
    
    
        
        

    
   
