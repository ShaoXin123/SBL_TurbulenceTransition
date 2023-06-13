function [s]=structure_f(a,n,maxlag)
nn=length(a);
for i=1:maxlag;   
   b1=a(1:nn-i+1);
   b2=a(i:nn);
   s(i)=mean((b1-b2).^n);
end

   
      
