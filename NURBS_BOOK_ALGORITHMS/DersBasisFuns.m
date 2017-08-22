    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Nurbs Book algorithm A2.3
    % create array of first k derivatives of all basis functions
    % there might be some mistakes in this algorithm
    function ders = DersBasisFuns(i,u,p,U,n) 
    % careful: changed order of U and n to be consistent with BasisFun
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