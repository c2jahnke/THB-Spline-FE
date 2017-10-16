function m20170801main()
    %parameters
    %addpatgit 
    %h('/home/laptop/Documents/_mat_files/')
    N = 5;
    resol = 0.1;
    basis = bSplBas(0,10,5,N,resol);
    %par = parameters();
     n = basis.n;
     p = basis.p;
     plotVector = basis.plotVector;
     sP = basis.sP;
     m = basis.m;
     knotVector = basis.knotVector;
%     N = bSplBasFun(2,basis);
%     nDer = 1;
% %     n = m - p - 1; % NURBS book: n = m - p -1, number of basis functions
% %    tableSpan = lookuptablespan(basis.knotVector,basis.m,basis.plotVector,basis.sP);
%      V0 = basis.generBasis();
%      ders = basis.gener1Deriv();

% %         %% test Curve evaluation
% %     % The minimum number of control points = p+1 + #(inner points)
% %     
% %     %hold off
    %% use 2D surface plots
     iU = 0;
     iV = 0;
     plotOne2Dbasis(p,m,knotVector,iU,iV,plotVector,sP,resol);
     figure
      plot2DBasis(p,m,knotVector,plotVector,sP,resol)

    Omega0 = [basis.a basis.b];
    Omega1 = [basis.a+floor((basis.a+basis.b)/4)+1 basis.b-2];
    h0 = 1;% generalize! So far works for pairs (1,0.5),(2,0.5)
    h1 = 0.5;
    
    
   
    TOL = 10e-5;
    [Char0, Char1, HB0, HB1, V0, V1, U, Ubar, Points, Qw] = OneDimHbRefinement(basis,Omega0,Omega1,h0,h1);
    hold off;
    
    plotPartBasis(basis.plotVector,size(HB0,2),HB0,'r-')
    hold on;
    plotPartBasis(basis.plotVector,size(HB1,2),HB1,'b-')
    
    
    [THB0, trunq, q] = OneDimThb(Omega0,Omega1,h0,h1,V0,V1,U,Ubar,Char1,TOL,basis);
    % to check partition of unity
    % sum(THB0(k,:))+sum(HB1(k,:))

    %% how to choose functions for truncation?
    % min basis function 
    % depending on degree -> for loop
    % calculate index for knotVector for which Omega1 begins/ends
    % calculate lowB = OneBasisFun highB for functions to truncate
    % solve as qlow = V1\lowB
    % calculate index for plotVector for which Omega1 begins/ends
    % a1 = (floor(Omega1(1)/resol)) a2 = (floor(Omega1(end)/resol))
    % set qlow to zero after a1 qhigh before a2
    movegui('east')
    figure
    plotPartBasis(basis.plotVector,size(THB0,2),THB0,'r-')
    hold on;
    plotPartBasis(basis.plotVector,size(HB1,2),HB1,'b-')
    %plotPartBasis(plotVector,size(V1,2),V1)
    
    %% check CurveKnotIns
    Points = zeros(basis.n,2);
    Points(:,1) = linspace(0,1, size(Points,1) );
    Points(:,2) = (Points(:,1) -0.1).^3 .*(Points(:,1) -0.9 ).^3;
    % does not work!
    ImpointCurvePlot(basis.n,basis.sP,basis.p,basis.knotVector,basis.plotVector,Points)
    %% knot insertion
%     np = n;
%     r = 1;
%     [k,s] = FindSpanMult(n,p,u,U);
%     [nq,UQ,QPoints] = CurveKnotIns(n,p,knotVector,Points,u,k,s,r); % works nicely!
%      figure
%      ImpointCurvePlot(nq,sP,p,UQ,plotVector,QPoints)
    %% use knot refinement
    X = [0.5 1.5 2.5 3.5 4.5];
    r= size(X,2)-1;
    [Ubar, Qw] = RefineKnotVectCurve(basis.n,basis.p,basis.knotVector,Points,X,r);
    %     figure
    % does not work so far
    ImpointCurvePlot(n+r,sP,p,Ubar,plotVector,Qw)
     
     
     
     
     
function erg = gQuadBasis(basis)
if (basis.p<5)
ngp = basis.p;
else
    ngp = 5;
end

switch(ngp)
    case 1
        s = 0;
        w = 2;
    case 2
        s = [-sqrt(1/3), sqrt(1/3)]';
        w = [1 ,1 ]';
    case 3
        s = [-0.774596669241, 0, 0.774596669241]';
        w = [5/9, 8/9, 5/9]';
    case 4
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856,0.861136311594053]';
        w = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]';
    case 5
        s = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664]';
        w = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]';
end

    erg = zeros(basis.n - basis.p, basis.p+1);
    for k = basis.p+1:(basis.n) % for loop over knots
        hS = (basis.knotVector(k+1)-basis.knotVector(k))/2;
        mid = (basis.knotVector(k+1)+basis.knotVector(k))/2;
        gPoints = s*hS+mid; % scale/transform gaussian points
        qGauss=zeros(size(gPoints,1),basis.p+1); % evaluate all p+1 functions on one element
        for i=1:size(gPoints,1)
           %  DersOneBasisFun(p,basis.m,basis.knotVector,k,gPoints(i),n)
          % OneBasisFun(p,basis.m,basis.knotVector,k,gPoints(i))
          % plotOneBasisFun(basis.sP,p,basis.plotVector,basis.knotVector,basis.m,k-p);
        qGauss(i,:) = BasisFun(k-1,basis.p,gPoints(i),basis.knotVector); % s(i) = v
    end
    erg(k-basis.p,:) = hS*(w'*qGauss)
    end
end
    function [THB0, trunq, q] = OneDimThb(Omega0,Omega1,h0,h1,V0,V1,U,Ubar,Char1,TOL,parameters)
    %V0 = generBasis(parameters.m,parameters.p,parameters.knotVector,parameters.knotVector,parameters.n,parameters.m);
    %V1 = generBasis(parameters.m,parameters.p,Ubar,Ubar,parameters.n,parameters.m);
    q = zeros(size(V1,2),size(V0,2));
    for k = 1 : size(V0,2) %% use Greville points or something similar so 
        %                   that this system is still uniquely solvable 
        q(:,k) = V1\V0(:,k);
    end
    % stabilize q
    % Konditionszahlabhängig
    q(q(:) < TOL) = 0;
    % charactaristic matrix Char0 Char1
    % coefficient matrix Points Qw
    % refinement matrix Q
    % For THB Splines
    p= parameters.p;
    trunq = q;
    % no good solution, so far works only for h1 = 0.5, h0 = 1
     if(h0 == 1 && h1 == 0.5 ) 
        trunq((Omega1(1)+p)/(h1)-(p-1):(Omega1(end))/(h1),:) = 0; % causes errors if h1 is very small;
      % trunq(2*(Omega1(1)+p)-(p-1):(Omega1(end))*2,:) = 0; % causes errors if h1 is very small;
     elseif(h0 == 2 && h1 == 1)
         trunq((Omega1(1)+p)/(h1)-(p-1)+p:(Omega1(end))/(h1),:)= 0;
%     generalize! write function that returns the support in knot indices    
     elseif(h0 == 4 && h1 == 2 ) % just works for p = 2
         trunq(floor((Omega1(1)+p+1)/(h1) + p):floor((Omega1(end)+1)/h1),:)= 0;
     else
         error('Truncation not yet implemented for h1,h0');
     end
    
    if(Omega1(end) == Omega0(end))
        trunq(Omega1(end)/h1:end,:) = 0;
    end
    
    if(Omega1(1) == Omega0(1))
        trunq(1:Omega1(end)/h1);
    end

    THB0 = V1*(trunq);
    THB0(THB0(:) < TOL) = 0;
    end
    
    
    function [Char0, Char1, HB0, HB1, V0, V1, U, Ubar, Points, Qw] = OneDimHbRefinement(parameters,Omega0,Omega1,h0,h1)
    % still not finished
    % treat case when Omega1(end) = Omega0(end) differently
    %switch nargin
    %case 4
    %    resol = h1/10;
    %case 3
    %    if (h0 ~= 0)
    %        h1 = h0/2;
    %        resol = h1/10;
    %    else
    %        fprintf('Error h0 ==0.')
    %    end
    %case 2
    %    fprintf('Not enough input arguments.')
    %end
    assert(~(Omega1(end) - Omega1(1) < parameters.p*h1),'Error: Omega1 to small for one basis function.');
        
    
    a = Omega0(1);
    b = Omega0(end);
    
    U = ConstrKnotVector(a,b,h0,parameters.p);
    m = size(U,2);
    n = m - parameters.p - 1;
    
    plotVector = a:parameters.resol:b;
    sP = size(plotVector,2);
    Points = zeros(n,2); % not needed
    Points(:,1) = linspace(0,1, size(Points,1) );
    Points(:,2) = sin(2*pi*Points(:,1));
    
    X = Omega0(1)+h1:h0:Omega0(end)-h1; %% not really general, work this out later
    r= size(X,2)-1;
    [Ubar, Qw] = RefineKnotVectCurve(n,parameters.p,U,Points,X,r); % the curve points Qw are not needed
    plotVectorUbar = a:parameters.resol:b; % work this out! Scale by factor h1/h0 -> Problem with truncation!
    V0 = generBasis(sP,parameters.p,plotVector,U,n,m); 
    V1 = generBasis(size(plotVectorUbar,2),parameters.p,plotVectorUbar,Ubar,size(Ubar,2)-parameters.p-1,size(Ubar,2));
    

    % has to depend on degree
    % so far works well for all degrees 
    temp1 = U(:) < Omega1(1)+(parameters.p-1)*h0;% not general!!! Change to h1 or h0, test!
    temp2 = U(:) >= Omega1(end)-(parameters.p-1)*h0;
    if (Omega1(end) == Omega0(end))
        temp2 = zeros(size(U))';
    end
    % calculate characteristic matrices
    Char0 = zeros(size(U,2)-parameters.p-1);
    Char0 =Char0 + diag(temp1(parameters.p:end -2)) + diag(temp2(2:end-parameters.p));

    temp3 = (Ubar(:) < Omega1(1) + (parameters.p-1)*h1);
    temp4 = (Ubar(:) >= Omega1(end) -(parameters.p-1)*h1 );
    if (Omega1(end) == Omega0(end))
        temp4 = zeros(size(Ubar))';
    end
    Char1 = eye(size(Ubar,2)-parameters.p-1);
    Char1 = Char1 -diag(temp3(parameters.p:end-2)) - diag(temp4(2:end-parameters.p));
    
    % calculate hierarchical basis
    HB1 = V1*Char1;
    HB0 = V0*Char0;
    
    %plotPartBasis(plotVector,size(HB0,2),HB0)
    %hold on;
    %plotPartBasis(plotVector,size(HB1,2),HB1)
    % truncation
    end
    
    function knotVector = ConstrKnotVector(a,b,h1,p)
    % simple constructor for knot vector
    temp = (b-a)/h1;
    m = temp + 2*p+1;
    knotVector = zeros(1,m);
    knotVector(1:p) = a;
    knotVector(p+1:m-p) = a:h1:b;
    knotVector(m-p+1:end) = b;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ploting ((T)H)BSpline-Basis
    function plotPartBasis(plotVector,n,C0,string)
        for ll = 1:n
            plot(plotVector,C0(:,ll),string);
            hold all
        end
        hold off
    end    
    
   
    %%%%%%%%%%%%%%%%%%%%%%
    % Num PDE 1, p 45 generate stiffness matrix and rhs vector
    function [A,b] = assemlin(x,f)
    N = length(x);
    n = N-1;
    Ae = [1, -1; -1, 1];
    A = zeros(N); b = zeros(N,1);
    for i = 1:n
        k = i-1; ke = k+1:k+2; % ??
        h = x(i+1) - x(i);
        A(ke,ke) = A(ke,ke) + 1./h * Ae;
        fe = feval(f,[x(i), x(i) + h/2, x(i+1)]);
        be = h/6*[fe(1) + 2*fe(2); 2*fe(2) + fe(3)];
        b(ke) = b(ke) + be;
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot2DBasis(p,m,knotVector,plotVector,sP,resol) %% arrange indices
        temp = 0;
        for iU = 0:knotVector(end)
            for iV = 0:knotVector(end)
                temp = temp +1;
                subplot(knotVector(end)+1,knotVector(end)+1,temp)
                plotOne2Dbasis(p,m,knotVector,iU,iV,plotVector,sP,resol);
                shading interp
            end
        end
        hold off
    end
    
    
    function plotOne2Dbasis(p,m,knotVector,iU,iV,plotVector,sP,resol)
        [X,Y] = meshgrid(knotVector(1):resol:knotVector(end),knotVector(1):resol:knotVector(end));

        CU = zeros(sP,1);
        for kk = 1 : sP
        CU(kk,1) = OneBasisFun(p,m,knotVector,iU,plotVector(kk));

        end
        CV = zeros(sP,1);
        for kk = 1 : sP
        CV(kk,1) = OneBasisFun(p,m,knotVector,iV,plotVector(kk));

        end
        Z = CU(:,1)*CV(:,1)'; % tensor product structure
        surf(X,Y,Z)
        %hold off
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot curve with the help of Impoints
    % thanks to Steffen
    function ImpointCurvePlot(n,sP,p,U,plotVector,Points)

        C = zeros( sP , 2 );
        for l = 1:sP % work around, make sure other points are not needed
            C(l,:) = CurvePoint(n,p,U,Points,plotVector(l));
        end
        plt_curve = plot(C(:,1),C(:,2),'r.');
        hold on;

        Impoints = {};
        for i = 1:n
           Impoints{i} = impoint(gca, Points(i,1) , Points(i,2) ); 
           addNewPositionCallback( Impoints{i} , @(newpos) updateCurvePlot(i,newpos) );
        end

        function updateCurvePlot( index , newpos )
            Points(index,:) = newpos;
             for k = 1:sP
                 C(k,:) = CurvePoint(n,p,U,Points,plotVector(k));
             end
            plt_curve.XData = C(:,1);
            plt_curve.YData = C(:,2);
            drawnow;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function plotDersOneBasisFun(sP,p,plotVector,knotVector,m,i,n)

    C = zeros(sP,n+1);
    % C(:) = DersOneBasisFun(p,m,knotVector,i,plotVector(:),n);

    for j = 1 : sP
        C(j,:) = DersOneBasisFun(p,m,knotVector,i,plotVector(j),n);
        
    end
    figure
    plot(plotVector,C(:,2),'y-');
    end
    % generate one basis function w r to plotVector
    function C = generOneBasisFun(sP,p,plotVector,knotVector,m,i)

    C = zeros(sP,1);
    for j = 1 : sP
        C(j,1) = OneBasisFun(p,m,knotVector,i,plotVector(j));

    end
    end    
    % plot one basis function
    function C = plotOneBasisFun(sP,p,plotVector,knotVector,m,i)

    C = zeros(sP,1);
    for j = 1 : sP
        C(j,1) = OneBasisFun(p,m,knotVector,i,plotVector(j));

    end
    plot(plotVector,C(:,1),'y-');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function innerKnotVector = inKnotVector(knotVector,p)

    innerKnotVector = knotVector(p+2:size(knotVector,2) -p-1);

    end
    
    % generate Basis array with respect to plotVector
    % highly inefficient! Instead use NURBS Book function A2.2 BasisFun and switch
    % indeces as with derivatives
    function C = generBasis(sP,p,plotVector,knotVector,n,m)
        C = zeros(sP,n);
        %C = BasisFun(i,p,u,U)
        for ll = 1:n
            for ks = 1 : sP
                C(ks,ll) = OneBasisFun(p,m,knotVector,ll-1,plotVector(ks));
            end
        end
    end

    % plot whole basis
    function C = plotBasis(sP,p,plotVector,knotVector,n,m)
        C = zeros(sP,n);
        for ll = 1:n
            for ks = 1 : sP
                C(ks,ll) = OneBasisFun(p,m,knotVector,ll-1,plotVector(ks));
            end
            plot(plotVector,C(:,ll));
            hold all
        end
        hold off
    end

    
    % plot all 1st derivatives
    function C = plot1Deriv(sP,p,tableSpan,plotVector,U,nDer)
    n = size(U,2) - p -1;
    C = zeros(sP,n);
    % use getPlotMatrix line 1040 ff (isogat)
    % use startX to switch index when entering new span index
    for i = 1 : sP
        startX = tableSpan(i) - tableSpan(1) +1;
        ders = DersBasisFuns(tableSpan(i),plotVector(i),p,nDer,U);
        % how can the functions be splitted such that several functions are
        % obtained? Observation: split C(:,j) at the j-th inner point!
        % Why
        for j = 1 : n % change to n
            if (j-tableSpan(i) + p -1 ) < p+1 && tableSpan(i) +1 - j  < p+1
                tmp = mod(j - startX, (p+1))+1;  % temporary solution
                C(i,j) = ders(2,tmp);
            end
        end
    end
    for j = 1: n
        plot(plotVector,C(:,j));
        hold all
    end

     hold off;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function tableSpan = lookuptablespan(knotVector,k,plotVector,s)
    %--------------------------------------------------------------
    % knotVector : knot vector
    % plotVector : plot vector
    % k          : size of knot vector
    % s          : size of plot vector
    %Output
    % tableSpan : row vector of same dimension as plotVector
    %--------------------------------------------------------------

    tableSpan = zeros(1,s);
    plotIndex = 1;
    leftIndex = 1; % left boundary of knot span
    for knotIndex = 2:k
        rightValue = knotVector(knotIndex);
        if(knotVector(leftIndex) == rightValue)
            continue
        end
        leftIndex = knotIndex -1;
        while(plotVector(plotIndex) < rightValue)
            tableSpan(plotIndex) = leftIndex;
            plotIndex = plotIndex +1;
        end
    end

    tableSpan(plotIndex) = leftIndex -1;
    tableSpan = tableSpan -ones(size(tableSpan));
    end
    % how to use lookuptablespan? Lookuptablespan returns a vector of the same
    % size as plotVector. The entries correspond to the span indices in
    % knotVector.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.1
    function knotSpanIndex = FindSpan(n,p,u,U)
    %--------------------------------------------------------------
    %function knotSpanIndex = FindSpan(n,p,u,U)
    % NURBS-Book (algorithm A2.1)
    % find the knot span index for one variable u
    %INPUT:
    % n          : number of basis function -1
    %        NURBS-Book: np+1 # basis, np max index (startindex 0)
    %        here        np   # basis and max index (startindex 1)
    % p          : degree of the basis functions
    % u          : evaluation point
    % U          : knot vector (row vector)
    %OUTPUT:
    % knotSpanIndex : index of knot span
    %--------------------------------------------------------------
    if (u == U(n+2))
        knotSpanIndex= n;
        return
    end
    low = p;
    high = n+1;
    mid = floor((low + high)/2);
    while (u <U(mid+1) || u >= U(mid+2) )
        if( u < U(mid+1))
            high = mid;
        else
            low = mid;
        end
        mid = floor((low+high)/2);
    end
    knotSpanIndex = mid;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.2
    % evaluate basis function
    function N = BasisFun(i,u,p,U) 
    % change order at every invocation
    %--------------------------------------------------------------
    %function N = BasisFun(i,p,u,U)
    % NURBS-Book (algorithm A2.2)
    % evalute nonzero basis functions
    %INPUT:
    % i          : current knotspan
    % u          : evaluation point
    % p          : degree of the basis functions
    % U          : knot vector (row vector)
    %OUTPUT:
    % N          : row vector (dim p+1)
    %              values of the basis function N_(i-p) ... N_(i)
    %              at the evaluation point u
    %--------------------------------------------------------------
    N=zeros(1,p+1);
    N(1)=1;
    left=zeros(1,p+1);
    right=zeros(1,p+1);
    for j=1:p
        left(j+1) = u-U(i+1-j+1); %% why j+1 and not j-1?
        right(j+1) = U(i+j+1)-u;
        saved = 0;
        for r=0:j-1
            temp = N(r+1)/(right(r+2)+left(j-r+1));
            N(r+1) = saved + right(r+2)*temp;
            saved = left(j-r+1)*temp;
        end
        N(j+1) = saved;
    end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Nurbs Book algorithm A2.3
    % create array of first k derivatives of all basis functions
    % there might be some mistakes in this algorithm
    function ders = DersBasisFuns(i,u,p,n,U)
    % Input:    i current knotspan
    %           u evaluation point
    %           p degree
    %           n degree of derivative (n<=p)
    %           U knotVector
    % Output    ders array, ders(k+1,j+1) stores the k-th derivative of N_i-p+j,p
    %           where 1 <=k <=n+1 and 1 <= j <= p+1
    ders = zeros(n+1,p+1);
    ndu = zeros(p+1,p+1); 
    ndu(1,1) = 1;
    left = zeros(1,p+1);
    right = zeros(1,p+1);
    for j=1:p
        left(j+1) = u- U(i+1-j+1);
        right(j+1) = U(i+j+1)-u;
        saved = 0;
        for r = 0 : (j-1) 
            % index shift wrt NURBS Book
            ndu(j+1,r+1) = right(r+2) + left(j-r+1);
            temp = ndu(r+1,j)/ndu(j+1,r+1);

            ndu(r+1,j+1) = saved + right(r+2)*temp;
            saved =left(j-r+1)*temp;
        end
        ndu(j+1,j+1) = saved;
    end
    for j= 0 : p % load basis functions
        ders(1,j+1) = ndu(j+1,p+1);
    end 
    % compute derivatives
    for r = 0 : p
        s1 = 0;
        s2 = 1;
        a = zeros(p+n,p+n);
        a(1,1) = 1;
        % compute kth derivative
        for k = 1 : n
            d = 0;
            rk = r - k;
            pk = p - k;
            if (r >=k )
                a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1); % Division by 0 happens! Index mistake? The plots seem reasonable.
                d =  a(s2+1,1)*ndu(rk+1,pk+1);
            end
            if rk >= -1
                j1 = 1;
            else
                j1 = -rk;
            end
            if r-1 <= pk
                j2 = k-1;
            else
                j2 = p-1;
            end
            for j = j1:j2
                a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1); % maybe rk+j+2? No. What else could be wrong?
                d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
            end
            if r <= pk
                a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                d = d+ a(s2+1,k+1)*ndu(r+1,pk+1);
            end
            ders(k+1,r+1) = d;
            j = s1;
            s1 = s2;
            s2 = j;
        end    
    end
    r= p;
    for k = 1:n
        for j = 0:p
            ders(k+1,j+1) =ders(k+1,j+1)* r;
        end
        r = r*(p-k );
    end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.4
    % evaluate one basis function
    function Nip = OneBasisFun(p,m,U,i,u)
    % compute the basis function Nip
    % Input     p,m,U,i,u
    % Output    Nip
    if ( (i == 0 && u == U(1)) || (i == m-p && u == U(m+1)))
        Nip = 1;
        return
    end
    if ( u < U(i+1) || u >= U(i+p+2))
        Nip = 0;
        return
    end
    N = zeros(1,p+1);
    % initialise basis functions
    for j= 0:p
        if ( u >= U(i+j+1) && u < U(i+j+2))
            N(j+1) = 1;
        else
            N(j+1) = 0;
        end
    end

    for k=1:p
        if( N(1) == 0 )
            saved = 0;
        else
            saved = ((u-U(i+1)) * N(1))/(U(i+k+1) - U(i+1));
        end
        for j = 0:(p-k)
            Uleft = U(i+j+2);
            Uright = U(i+j+k+2);
            if(N(j+2) == 0)
                N(j+1) = saved;
                saved =0;
            else
                temp = N(j+2)/(Uright - Uleft);
                N(j+1) = saved + (Uright -u)*temp;
                saved = (u - Uleft)*temp;
            end
        end
    end
    Nip = N(1);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.5
    % evaluate first k derivatives of one basis function
    % maybe there is still some index error
    function ders = DersOneBasisFun(p,m,U,i,u,n)
    % compute derivatives of basis function Nip
    % Input     p,m,U,i,u,n
    % Output    ders
    ders = zeros(1,n+1);
    if ( u< U(i+1) || u >= U(i+p+2) )
            return;
    end
    N = zeros(p+2,n+1); % changed to p+2, maybe p+1 suffices
    for j=0:p % initialize 0-degree functions
        if (u >= U(i+j+1) && u < U(i+j+2))
            N(j+1,1) = 1;
        else
            N(j+1,1) = 0;
        end
    end
    for k = 1:p % compute full triangular table?
        if (N(1,k) == 0)
            saved = 0;
        else
            saved = ((u-U(i+1))*N(1,k))/(U(i+k+1)-U(i+1));
        end
        for j = 0:(p-k)
            Uleft = U(i+j+2);
            Uright = U(i+j+k+2);
            if( N(j+2,k) == 0)
                N(j+1,k+1) = saved;
                saved = 0;
            else
                temp = N(j+2,k)/(Uright - Uleft);
                N(j+1,k+1) = saved + (Uright - u) * temp;
                saved = (u-Uleft)*temp;
            end
        end
    end
    ders(1) = N(1,p+1);
    ND = zeros(1,n+1);
    for k=1:n
        for j=0:k
            ND(j+1) = N(j+1,p-k+1);
        end
        for jj = 1:k
            if (ND(1) == 0)
                saved = 0;
            else
                saved = ND(1)/(U(i+p-k+jj+1) - U(i+1));
            end
            for j=0:(k-jj) % index as in Nurbs Book
                Uleft = U(i+j+2);
                Uright = U(i+j+p+jj+1); % changed from +2 to +1
                if( ND(j+2) == 0)
                    ND(j+1) = (p-k+jj)*saved;
                    saved = 0;
                else
                    temp = ND(j+2)/(Uright -Uleft);
                    ND(j+1) = (p-k+jj)*(saved-temp);
                    saved = temp;
                end
            end
        end
        ders(k+1) = ND(1);
    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A3.1 adjusted
    % evaluate a point on a B-Spline curve

    function C = CurvePoint(n,p,U,P,u)
    % Compute a curve point
    % Input:    n number of basis functions -1, (n = m - p -1)
    %           p spline degree
    %           U knotVector
    %           P vector of control points
    % Output:   C curve point

    span = FindSpan(n,p,u,U);
    N = BasisFun(span,u,p,U); %BasisFun(i,u,p,U) 
    if( span+1 > length(P) )
        C = P(end,:);
    else
        C = N * P( span+1-p : span+1 , : );
    end
    %C = 0;
    %for i=0:p
    %    C = C + N(i+1)*P(span - p+i+1);
    %end
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A5.1 (partly vectorized)
    % evaluate a point on a B-Spline curve
    function [nq,UQ,Qw] = CurveKnotIns(np,p,UP,Pw,u,k,s,r,nq)
    
    %% compute new curve from knot insertion
    % Input:    np = mq-p-1 (number of basis functions before insertion)
    %           p  degree
    %           UP knot vector before insertion
    %           Pw 2d control points
    %           u  knot to insert
    %           k  knot span such that u \in [u_k,u_k+1] % check this
    %           s  multiplicity of knot to be inserted ? r+s <= p
    %           r  r times insertion
    % Output:   nq = mp - p -1 (number of basis functions after insertion)
    %           UQ knot vector after insertion
    %           Qw control points after insertion
    
    mp = np+p+1;
    nq = np+r; %??
    Qw = zeros(size(Pw,1)+r,2);
    UQ = zeros(1,size(UP,2)+r);
    % Load new knot vector
    UQ(1:k+1) = UP(1:k+1);
    UQ(k+2:k+r+1) = u;
    UQ(k+r+2:mp+1) = UP(k+2:mp);
    % save unaltered control points
    Qw(1:k-p+1,:) = Pw(1:k-p+1,:);
    Qw(k-s+r+1:np+r,:) = Pw(k-s+1:np,:);
    Rw = zeros(p-s+1,2);
    Rw(1:p-s+1,:) = Pw(k-p+1: k-s+1,:);
    
    for jj = 1:r
        L = k-p+jj;
        for i =0:(p-jj-s)
            alpha = (u-UP(L+i+1))/(UP(i+k+2)-UP(L+i+1));
            Rw(i+1,:) = alpha* Rw(i+2,:) + (1-alpha)*Rw(i+1,:);
        end
        Qw(L+1,:) = Rw(1,:);
        Qw(k+r-jj-s+1,:) = Rw(p-jj-s+1,:);
    end
    
        Qw(L+2:k-s,:) = Rw(2: k-s-L,:);
    end
    
    
    %%%%%%%%%%%%%%%
    % Algorithm to find span k and multiplicity s of u within knot Vector U
    % Input:    n = m-p-1
    %           p degree
    %           u current knot
    %           U knot vector
    % Output:   k span index
    %           s multiplicity
    function [k,s] = FindSpanMult(n,p,u,U)
        k = FindSpan(n,p,u,U);
        s = 0;
        if u == U(1)
            for i = 1 : p+1
                if u == U(i)
                    s = s+1;
                end
            end
        else
            for i = 1 : p+1
                if u == U(k+i)
                s = s+1;
                end
            end
        end
    end
   

   
end