 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Algorithm A5.4

    function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
    % refine curve knot vector
    % Input:    n = m - p -1
    %           p degree
    %           U knot vector
    %           Pw control points
    %           X set(vector) of knots to be inserted into U
    %           r size of X
    % Output    Ubar new knot vector
    %           Qw knot control points
    n = n-1; %ok, but why???
    m = n+p+1; %ok
    a = FindSpan(n,p,X(1),U); %ok
    b = FindSpan(n,p,X(r+1),U)+1; %ok   
    % initializations
    Ubar = zeros(1,m+r+2); %ok
    Qw = zeros(size(Pw,1)+r+1,2); % ?
    
    Qw(1:a-p+1,:) = Pw(1:a-p+1,:); %correct
    Qw(b+r+1:n+r+2,:) = Pw(b:n+1,:); % not sure
    Ubar(1:a+1) = U(1:a+1); % correct
    Ubar(b+p+r+1:m+r+2) = U(b+p:m+1); % not sure, check indices
    
    i = b+p;%ok
    k = b+p+r+1;%ok
    for j = r:-1:0 % index as in NURBS book
        while(X(j+1) <= U(i+1) && i > a)
            Qw(k-p,:) = Pw(i-p,:); % 
            Ubar(k+1) = U(i+1);
            k = k-1; i = i-1;
        end
        Qw(k-p,:) = Qw(k-p+1,:); % 
        for l = 1:p
            ind = k-p+l; % ell
            alfa = Ubar(k+l+1) - X(j+1);
            if (abs(alfa) == 0)
                Qw(ind,:) = Qw(ind+1,:);
            else
                alfa = alfa/(Ubar(k+l+1) - U(i-p+l+1)); %
                Qw(ind,:) = alfa*Qw(ind,:) + (1-alfa)*Qw(ind+1,:);
            end
        end
        Ubar(k+1) = X(j+1);
        k = k-1;
    end
    end
    