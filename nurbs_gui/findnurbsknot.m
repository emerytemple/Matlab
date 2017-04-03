function [u] = findnurbsknot(choice,clamped,n,k,P)

    u(k) = 0.0;
    u(n+1) = 1.0;
    if(choice == 1) % uniform
        for i = k+1:n
            u(i) = (i-k)/(n+1-k);
        end
    else
        for i = 1:n-1
            el(i) = sqrt(((P(i+1,1)-P(i,1))^2)+((P(i+1,2)-P(i,2))^2));
            elt(i) = sqrt(el(i));
        end

        L = 0.0;
        Lt = 0.0;
        for i = 1:n-1
            L = L + el(i);
            Lt = Lt + elt(i);
        end
        for i = k+1:n+1
            m = i-k;
            top = 0.0;
            topt = 0.0;
            for j = 1:k-1
                top = top + el(m-1+j);
                topt = topt + elt(m-1+j);
            end
            if(choice == 2) % chordlength
                u(i) = u(i-1) + (top/L);
            elseif(choice == 3) % centripetal
                u(i) = u(i-1) + (top/Lt);
            else
                error = 1
            end
        end
    end
	for i = n+2:n+k
        if(clamped == 1) % clamped == true
            u(i) = u(n+1);
        elseif(clamped == 0)
            u(i) = u(i-1)+u(n+1)-u(n);
        else
            error = 1
        end
    end
	for i = k-1:-1:1
        if(clamped == 1) % clamped = true
            u(i) = u(k);
        elseif(clamped == 0)
            u(i) = u(i+1)-u(k+1)-u(k);
        else
            error = 1
        end
	end    
    if(choice ~= 1)
        val = u(n+1);
        for i = 1:n+k
            u(i) = u(i)/val;
        end
    end    
end