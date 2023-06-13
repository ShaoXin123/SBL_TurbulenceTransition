function [s]=autocorr(a,maxlag,fs)
nn=length(a);
m=maxlag;
a1=(a-mean(a))/std(a);
n=1;
for i=1:fs:m;
   b1=(a1(1:nn-i+1)-mean(a1(1:nn-i+1)))/std(a1(1:nn-i+1));
   b2=(a1(i:nn)-mean(a1(i:nn)))/std(a1(i:nn));
   s(n)=mean(b1.*b2);n=n+1;
end

   
      
