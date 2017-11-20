
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
    