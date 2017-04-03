function [C] = decasteljau(u,n,np,P)
    
    C = zeros(np+1,2);
    B = zeros(n,n,np+1,2);

    for i = 1:n
        for j = 1:2
            B(1,i,:,j) = P(i,j);
        end
    end

    for i = 2:n
        for j = 1:n-i+1
            for k = 1:np+1
                for m = 1:2
                    B(i,j,k,m) = ((1-u(k))*B(i-1,j,k,m))+(u(k)*B(i-1,j+1,k,m));
                end
            end
        end
    end
    
    for i = 1:np+1
        for j = 1:2
            C(i,j) = B(n,1,i,j);
        end
    end
end