function [val] = bisection(x,part,tol,iter)
    n1 = 0;
    n2 = 1;
    for i = 1:iter
        
        % use decasteljau algorithm to find xvalue for a given param value
        u = (n1+n2)/2;
        for i = 1:4
                B(1,i) = part(i);
        end

        for i = 2:4
            for j = 1:4-i+1
                B(i,j) = ((1-u)*B(i-1,j))+(u*B(i-1,j+1));
            end
        end

        xmid = B(4,1,1);        
        
        % continue with rest of bisection algorithm
        mid = (n1+n2)/2;
        if(x < xmid)
            n2 = mid;
        elseif(x > xmid)
            n1 = mid;
        else
            val = mid;
            break;
        end
        if(abs(n2-n1) <= tol) % need to use relative tolerance values here
            val = (n1+n2)/2;
            break;
        end
        val = (n1+n2)/2;
    end
end