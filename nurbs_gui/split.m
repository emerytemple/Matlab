function [C1,C2] = split(t,n,P) % split a bezier curve
    
    B = zeros(n,n,2);
    C1 = zeros(n,2);
    C2 = zeros(n,2);
    
    for i = 1:n
        for j = 1:2
            B(1,i,j) = P(i,j);
        end
    end
    
    for i = 2:n
        for j = 1:n-i+1
            for m = 1:2
                B(i,j,m) = ((1-t)*B(i-1,j,m))+(t*B(i-1,j+1,m));
            end
        end
    end
    
    for i = 1:n
        for j = 1:2
            C1(i,j) = B(i,1,j);
            C2(i,j) = B(n-i+1,i,j);
        end
    end

end