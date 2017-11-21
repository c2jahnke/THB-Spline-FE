classdef PoissSolv < handle
    properties %(GetAccess = public, SetAccess = private)  
        MlBas = []; % hierarchical multi level basis (hbSplBasML, thbSplBasML)
        f = []; % source term, function handle
        iLvl = []; % level index
        iBasisFctInd = []; % basis function index
        
    end

     methods (Access = public)
        function obj = PoissSolv(MlBas,f)
            %--------------------------------------------------------------
            % constructor for class PoissSolv
            % MlBas: multi level base such as thbSplBasML
            % f: right hand side
            %--------------------------------------------------------------
      
            if nargin >0
                obj.MlBas = MlBas;
                obj.f = f;
            end
        end
        function [Stiffn, rhs, iLvl,iBasisFctInd] = assembleMl(obj)
            % assemble stiffness matrix and rhs
            allKnots = obj.MlBas.getAllKnots;

            nOF = obj.MlBas.nOF; % number of functions
            nOE =length(allKnots) -1; % number of elements
            Stiffn = zeros(nOF); %basis.n
            ngp = max(0.5*(obj.MlBas.levelBas{1}.p-1)^2-1,5);%max(obj.MlBas.levelBas{1}.p+1,sqrt(obj.MlBas.levelBas{1}.p^2 -2*obj.MlBas.levelBas{1}.p+1));
            % first part to integrate right hand side sufficiently accurate, second
            % part to generate stiffness matrix exactly
            rhs = zeros(nOF,1); % number of basis functions!


            for el = 1 : nOE % loop over elements
                % export evaluation to seperate function
                [bVal,gradVal,lIndex,s,w] = evalEl(obj.MlBas,el);
                % assemble element stiffness matrix
                elRhs = zeros(size(bVal,2),1); 
                elStiff = zeros(size(bVal,2));
                elSInd = cell(size(bVal,2));
                elRInd =cell(size(bVal,2),1);
                for ii0 = 1 : size(bVal,2) 
                    elRhs(ii0) = sum(w.*obj.f(s).*bVal(:,ii0)); % changed to arrayfun for non vectorized functions
                    elRInd{ii0} = lIndex(:,ii0);
                    elStiff(ii0,ii0) = sum(w.*gradVal(:,ii0).^2);
                    elSInd{ii0,ii0} = [lIndex(:,ii0) lIndex(:,ii0)];
                    for jj = ii0+1 : size(bVal,2) 
                        elStiff(ii0,jj) = sum(w.*gradVal(:,ii0).*gradVal(:,jj));
                        elSInd{ii0,jj} = [lIndex(:,ii0) lIndex(:,jj)];
                        elStiff(jj,ii0) = elStiff(ii0,jj);
                        elSInd{jj,ii0} = elSInd{ii0,jj};
                    end      
                end
                % generation of element stiffness matrix and of index matrix done!

                for l = 1 : obj.MlBas.level
                    for k = 1:length(obj.MlBas.levelBas{l}.activeIndex)
                        if(l >1) % raise index according to lower levels
                            ind_1 = 0;
                            for tInd = 1 : l-1
                                ind_1 = ind_1 + length(obj.MlBas.levelBas{tInd}.activeIndex);
                            end
                            ind_1 = ind_1 + k;
                        else
                            ind_1 = k; % index for basis functions
                        end
                        for ii1 = 1 : size(bVal,2)
                            if(elRInd{ii1} == [l;obj.MlBas.levelBas{l}.activeIndex(k)])
                               rhs(ind_1) = rhs(ind_1) + elRhs(ii1);
                               iLvl(ind_1) = l; % keep track of all indices
                               iBasisFctInd(ind_1) = obj.MlBas.levelBas{l}.activeIndex(k);
                            end
                         end

                        for ll = l : obj.MlBas.level
                            for kk = 1:length(obj.MlBas.levelBas{ll}.activeIndex) 
                                if(ll >1) % raise index according to lower levels
                                    ind_2 = 0;
                                for tInd = 1 : ll-1
                                    ind_2 = ind_2 + length(obj.MlBas.levelBas{tInd}.activeIndex);
                                end
                                ind_2 = ind_2 + kk;
                                else
                                ind_2 = kk;
                                end
                                for  ii = 1 : size(bVal,2)
                                    for jj = ii : size(bVal,2)
                                    if(elSInd{ii,jj} == [l ll; obj.MlBas.levelBas{l}.activeIndex(k)...
                                            obj.MlBas.levelBas{ll}.activeIndex(kk)])
                                        Stiffn(ind_1,ind_2) = Stiffn(ind_1,ind_2) + elStiff(ii,jj);
                                        Stiffn(ind_2,ind_1) = Stiffn(ind_1,ind_2);
                                    end
                                    end
                                end
                            end 
                        end
                    end
                end
            end 
            obj.iLvl = iLvl;
            obj.iBasisFctInd = iBasisFctInd;
        end
        function y = solveSyst(obj,Stiffn,rhs,dBC)
            % solve system with stiffness matrix and rhs
            [sdLvl, sdInd] = obj.MlBas.getActiveFctIndU(obj.MlBas.levelBas{1}.a);
            sEval = obj.MlBas.evalBasis(obj.MlBas.levelBas{1}.a);
            sIndEval = find(sEval);
            sFLvl = find( sdLvl(sIndEval(1)) == obj.iLvl);
            sFInd = find(sdInd(sIndEval(1)) == obj.iBasisFctInd);
            sInd = intersect(sFLvl,sFInd); % start index

            [edLvl, edInd] = obj.MlBas.getActiveFctIndU(obj.MlBas.levelBas{1}.b-obj.MlBas.levelBas{1}.resol);
            eEval = obj.MlBas.evalBasis(obj.MlBas.levelBas{1}.b-obj.MlBas.levelBas{1}.resol);
            eIndEval = find(eEval);
            eFLvl = find( edLvl(eIndEval(end)) == obj.iLvl); % 
            eFInd = find(edInd(eIndEval(end)) == obj.iBasisFctInd);
            eInd = intersect(eFLvl,eFInd); % end index

            if(dBC.mixed == true)
                A = zeros(obj.MlBas.nOF +1);
                b = zeros(obj.MlBas.nOF +1,1);
                if(strcmp(dBC.east,'Dirichlet') && strcmp(dBC.west,'Neumann'))
                    % does not work this way!
                    A(1,eInd+1) = 1;
                    A(eInd+1,1) = 1;
                    b(1,1) = dBC.eVal; % Dichlet BC
                   A(2:end,2:end) = Stiffn;
                    b(2:end) = rhs;
                    u = A\b; 
                    y = u(2:end);   
                elseif(strcmp(dBC.west,'Dirichlet') && strcmp(dBC.east,'Neumann'))
                    A(1,sInd+1) = 1;
                    A(sInd+1,1) = 1;
                    b(1,1) = dBC.wVal; % Dichlet BC
                    A(2:end,2:end) = Stiffn;
                    b(2:end) = rhs;
                    u = A\b; 
                    y = u(2:end);
                end
            elseif(strcmp(dBC.east,'Dirichlet'))
                A = zeros(obj.MlBas.nOF +2);
                b = zeros(obj.MlBas.nOF +2,1);
                A(1,sInd+2) = 1;
                A(sInd+2,1) = 1;
                b(1,1) = dBC.wVal;
                A(2,eInd+2) = 1;
                A(eInd+2,2) = 1;
                b(2,1) = dBC.eVal;

                A(3:end,3:end) = Stiffn;
                b(3:end) = rhs;
                u = A\b; 
                y = u(3:end);
            elseif(strcmp(dBC.west,'Neumann'))
                y = Stiffn\b;

            end
        end
        
        function uh =  generSolThb(obj,y)
            % generate the solution using THB-splines
            D = [];
            for ll = 1 : obj.MlBas.level
                C = zeros(obj.MlBas.levelBas{ll}.sP,length(obj.MlBas.levelBas{ll}.activeIndex));
                for l = 1:length(obj.MlBas.levelBas{ll}.activeIndex)
                    % create basis function for thb-splines
                    basisF = thbSplBasFun(obj.MlBas.levelBas{ll}.activeIndex(l),obj.MlBas,ll);
                    C(:,l) = basisF.generOneBasisFun;
                end
                D = [D C(1:2^(ll-1):size(C,1),:)];% D = [D C]
            end
            uh = y(1)*D(:,1);
            for k = 2 : obj.MlBas.nOF
                uh = uh + y(k)*D(:,k);
            end
            %plot(obj.MlBas.levelBas{1}.plotVector(1:end-1),uh(1:end-1),'r')    
        end
        function uh =  generSol(obj,y)
            D = [];
            for ll = 1 : obj.MlBas.level
                C = zeros(obj.MlBas.levelBas{ll}.sP,length(obj.MlBas.levelBas{ll}.activeIndex));
                for l = 1:length(obj.MlBas.levelBas{ll}.activeIndex)
                    basisF = bSplBasFun(obj.MlBas.levelBas{ll}.activeIndex(l),obj.MlBas.levelBas{ll});
                    C(:,l) = basisF.generOneBasisFun;
                end
                D = [D C(1:2^(ll-1):size(C,1),:)];% changed 2^(ll-1) from ll
            end
            uh = y(1)*D(:,1);
            for k = 2 : size(D,2)
                uh = uh + y(k)*D(:,k);
            end
            plot(obj.MlBas.levelBas{1}.plotVector,uh,'r')
        end
        
        function [H1err,L2err,Enerr] = errCalc(obj,g,y)
        % Calculate H1, L2, Energy norm
            L2err = 0; H1err = 0; Enerr = 0;
            h = 10^-10;
            gDer = @(x) (g(x+h) - g(x-h))/(2*h);
            for el = 1:obj.MlBas.nOE
                %[bVal,gradVal,lIndex,s,w] = evalEl(obj,el);
                allKnots = obj.MlBas.getAllKnots;
                ngp = 10;
                [s,w]=lgwt(ngp,allKnots(el),allKnots(el+1));
                uexVal = feval(g,s);
                %calculate dervative using central difference scheme, upwind scheme
                uexGrad = feval(gDer, s);
                uappGrad = [];
                for k = 1 : length(s)
                    approx = evalDersSol(obj.MlBas,s(k),obj.iBasisFctInd,obj.iLvl,y);
                    uappGrad = [uappGrad approx];
                end
                L2err = L2err + ((uappGrad(1,:)-uexVal(:)').^2*w);%L2err + (h/6)*((ue(1)-q(i))^2+4*(ue(2)-0.5*(q(i)+q(i+1)))^2+(ue(3)-q(i+1))^2);
                Enerr = Enerr +((uappGrad(2,:)- uexGrad(:)').^2*w);%Enerr + (h/6)*((udxe(1)-qdx)^2+4*(udxe(2)-qdx)^2+(udxe(3)-qdx)^2);
            end
            H1err = sqrt(L2err + Enerr);
            L2err = sqrt(L2err);
            Enerr = sqrt(Enerr);
        end
        
     end
end
     