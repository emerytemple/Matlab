function [u] = findknot(choice,n,P)
    if(choice == 1) % uniform
        for i = 1:n
            u(i) = (i-1)/(n-1);
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
        u(1) = 0.0;
        for i = 2:n
            if(choice == 2) % chordlength
                u(i) = u(i-1)+(el(i-1)/L);
            elseif(choice == 3) % centripetal
                u(i) = u(i-1)+(elt(i-1)/Lt);
            else
                error = 1
            end
        end
    end
end
