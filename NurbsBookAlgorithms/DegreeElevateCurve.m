
 % Test
    % NURBS Book, Algorithm A5.9 

    function [nh,Uh,Qw] =DegreeElevateCurve(n,p,U,Pw,t) % TODO: Test and Debug!!!
    % degree refine curve
    % Input:    n = m - p -1
    %           p degree
    %           U knot vector
    %           Pw control points
    % Output    t refinement to degree p+t
    %           nh number of functions
    %           Uh 
    %           Qw new control points

    m = n+p+1;
    ph = p+t;
    ph2 = ph/2;
    % Compute Bezier degree elevation coefficients
    bezalfs = zeros(ph+1,p+1);
    bezalfs(1,1) = 1; bezalfs(end,end) = 1;

    for i = 1:ph2
        inv = 1/nchoosek(ph,i);
        mpi = min(p,i);
        for j =Max(0,i-t):mpi
            bezalfs(i+1,j+1) = inv*nchoosek(p,j)*nchoosek(t,i-j);
        end
    end
    for i = ph2+1:ph-1
        mpi = Min(p,i);
        for j = Max(0,i-t):mpi
            bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);
        end
    end
    mh = ph; kind = ph+1;
    r = -1; a = p;
    b = p+1; cind = 1;
    ua = U(1);
    Qw(1,:) = Pw(1,:);
    for i = 0:ph
        Uh(i+1) = ua;
    end
    bpts = []; % fix size later!
    for i = 0:p
        bpts(i+1,:) = Pw(i+1,:);
        while (b < m)
            i = b
            while (b<m && U(b) == U(b+1))
                b = b+1;
            end
            mu1 = b-i+1;
            mh = mh + mu1 + t;
            ub = U(b);
            oldr = r;
            r = p-mu1;
            % insert knot u(b) r times
            if(oldr > 0)
                lbz = (oldr+2)/2;
            else
                lbz = 1;
            end
        if(r>0)
            rbz = ph-(r+1)/2;
        else
            rbz = ph;
        end
        if(r>0)
            % insert knot to get Bezier segment
            numer = ub - ua;
            for (k = mu1:-1:p
                alfs(k-mu1) = number/(U(a+k+1) -ua);
            end
            for j=1:r
                save = r-j; s = mu1+j;
                for k= p:-1:s
                    bpts(k+1) ) alfs(k-s+1)*bpts(k+1) + (1-alfs(k-s+1)*bpts(k);
                end
                Nextbpts(save) = bpts(p);
            end
        end
        for (i = lbz:ph) % Degree elevate Bezier
            ebpts(i+1) = 0;
            mpi = min(p,i);
            for j = max(0,i-t):mpi
                ebpts(i+1) = ebpts(i+1) + bezalfs(i+1,j+1)*bpts(j+1);
            end
        end
        if(oldr > 1) % remove knot u = U(a) oldr times
            first = kind-2; last = kind;
            den = ub -ua;
            bet = (ub-Uh(kind-1))/den;
            for tr=1:oldr-1% knot removal
                i = first; j = last; kj = j-kind+1;
                while(j-i > tr) %compute new control points
                    if(i<cind)
                        alf = (ub-Uh(i+1)/(ua-Uh(i+1));
                        Qw(i+1) = alf*Qw(i+1) + (1-alf)*Qw(i);
                    end
                    if(j>= lbz)
                        if (j-tr <= kind-ph+oldr)
                            gam = (ub - Uh(j-tr+1))/den;
                            ebpts(kj) = gam*ebpts(kj+1) + (1-gam)*ebpts(kj+2);
                        else
                            ebpts(kj+1) = bet*ebpts(kj+1)+(1-bet)*ebpts(kj+2);
                        end
                    end
                    i = i+1; j = j-1; kj = kj-1;
                end
                first = first-1; last = last +1;
            end
        end
        if( a ~= p ) % load knot ua
            for (i = 0: ph-oldr-1)
                Uh(kind+1) = ua; kind = kind+1;
            end
        end
        for j = lbz:rbz
            Qw(cind+1) = ebpts(j+1); cind = cind+1;
        end
        if(b<m)
            for j = 0:r-1
                bpts(j+1) = Nextbpts(j+1);
            end
            for j = r:p
                bpts(j+1) = Pw(b-p+j+1);
            end
            a = b; b = b+1; ua = ub;
        else % end knot
            for i = 0:ph
                Uh(kind+i) = ub;
            end
        end
        end
    nh = mh - ph -1;

    end
    
end