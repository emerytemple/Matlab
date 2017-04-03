function [b,d] = c2cubic(select,u,n,P)
     r = zeros(n,2);
     A = zeros(n,n);
     b = zeros(n-1,2);     
         
     for i = 1:n-1
         del(i) = u(i+1)-u(i);
     end  
     del(n) = 0.0;
     
     for i = 2:n-1
         if(i == 2)
             del_im2 = 0.0;
         else         
             del_im2 = del(i-2);
         end
         alpha(i) = ((del(i))^2)/(del_im2+del(i-1)+del(i));
         beta(i) = ((del(i)*(del_im2+del(i-1)))/(del_im2+del(i-1)+del(i)))+((del(i-1)*(del(i)+del(i+1)))/(del(i-1)+del(i)+del(i+1)));
         gamma(i) = ((del(i-1))^2)/(del(i-1)+del(i)+del(i+1));
     end

     for i = 2:n-1
     	A(i,i-1) = alpha(i);
     	A(i,i) = beta(i);
     	A(i,i+1) = gamma(i); 
         for j = 1:2
             r(i,j) = (del(i-1)+del(i))*P(i,j);
         end             
     end
     
     if(select == 1) % "natural" BC
         A(1,1) = 1.0;
         A(n,n) = 1.0;

         b(2,1) = P(1,1)+((1/3)*(P(2,1)-P(1,1)));
         b(2,2) = P(1,2)+((1/3)*(P(2,2)-P(1,2)));
         r(1,1) = b(2,1);
         r(1,2) = b(2,2);    

         b((3*(n-1)),1) = P(n,1)-((1/3)*(P(n,1)-P(n-1,1)));
         b((3*(n-1)),2) = P(n,2)-((1/3)*(P(n,2)-P(n-1,2))); 
         r(n,1) = b((3*(n-1)),1);
         r(n,2) = b((3*(n-1)),2);   
     elseif(select == 2) % "not-a-knot" BC
         A(1,1) = -1.0;
         A(1,2) = 1.0;
         A(n,n-1) = -1.0;
         A(n,n) = 1.0;

         r(1,1) = P(2,1) - P(1,1);
         r(1,2) = P(2,2) - P(1,2);
         r(n,1) = P(n,1) - P(n-1,1);
         r(n,2) = P(n,2) - P(n-1,2);
     else
         error = 1;
     end

     % calculate NURBS control points
     d(:,1) = A\r(:,1);
     d(:,2) = A\r(:,2);
 
     % calculate Bezier control points
	 b(1,1) = P(1,1);
	 b(1,2) = P(1,2);     
     b((3*(n-1)+1),1) = P(n,1);
     b((3*(n-1)+1),2) = P(n,2);     
     for i = 1:n-2         
          for j = 1:2               
               if(i == 1)
                   del_im2 = 0.0;
               else         
                   del_im2 = del(i-1);
               end             
               b((3*i)+1,j) = P(i+1,j); 
               b(3*i,j) = ((del(i+1)*d(i,j))+((del_im2+del(i))*d(i+1,j)))/(del_im2+del(i)+del(i+1));
               b((3*i)+2,j) = (((del(i+1)+del(i+2))*d(i+1,j))+(del(i)*d(i+2,j)))/(del(i)+del(i+1)+del(i+2));
          end
     end
     
     % shift array and add end NURBS control points
    for i = length(d)+1:-1:2
        for j = 1:2
            d(i,j) = d(i-1,j);
        end
    end
    d(1,1) = P(1,1);
    d(1,2) = P(1,2);
    nnn = length(d) + 1;
    d(nnn,1) = P(n,1);
    d(nnn,2) = P(n,2);     
end