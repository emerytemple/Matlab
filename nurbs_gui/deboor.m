function [c] = deboor(n,k,np,t,u,w,cp)
     t(1) = 0.00001;
     t(np+1) = 0.9999;

     % homogeneous coordinates
     for i = 1:n
         cp(i,1) = w(i)*cp(i,1);
         cp(i,2) = w(i)*cp(i,2);         
         cp(i,3) = w(i);
     end

     for i = 1:n
          for j = 1:np+1
               for m = 1:3
                    d(1,i,j,m) = cp(i,m);
               end
          end
     end
     for m = 1:np+1
          for s = 1:n+k
               if((u(s) <= t(m)) && (t(m) < u(s+1)))
                    I(m) = s;
               end
          end
          for j = 1:3
               for r = 2:k
                    for i = (I(m)-(k-1)):(I(m)-(r-1))
                         d(r,i,m,j) = (((u(i+k)-t(m))/(u(i+k)-u(i+r-1)))*d(r-1,i,m,j))+(((t(m)-u(i+r-1))/(u(i+k)-u(i+r-1)))*d(r-1,i+1,m,j)); 
                    end
               end
          end
     end
     for j = 1:np+1
          for m = 1:3
               ct(j,m) = d(k,I(j)-(k-1),j,m);
          end
     end
     for j = 1:np+1
          c(j,1) = ct(j,1)/ct(j,3);
          c(j,2) = ct(j,2)/ct(j,3);
     end
end