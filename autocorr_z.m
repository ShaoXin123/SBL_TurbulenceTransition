function [s]=autocorr_z(a,iz,zbin,L)
nn=length(a(:,1));nz=length(a(1,:));


s=zeros(1,nz)*nan;
for i=zbin
    n=0;s(i)=0;
    for inn=1:L:nn-L
        b1=a(inn:inn+L-1,iz);
        b2=a(inn:inn+L-1,i);
        b=corrcoef(b1,b2);
        s(i)=s(i)+b(1,2);
        n=n+1;
    end    
   s(i)=s(i)/n;  
end

   
      
