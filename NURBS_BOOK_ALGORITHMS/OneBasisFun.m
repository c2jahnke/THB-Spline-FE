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
