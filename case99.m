
zcase99=[1.5 5 10 20 30 40 50 55];zcase99tidu=[5 15 25 35 45 55];
L=15*60*20;Ltidu=L/20;fs=20;
zcase99array=ones(600*30*60*20/L,8)*(zcase99'*ones(1,8).*eye(8,8));


z_L=zeros(fix(length(u(:,1))/L),8)*nan;buoyancy=zeros(fix(length(u(:,1))/L),8)*nan;
wnd=zeros(fix(length(u(:,1))/L),8)*nan;
epi=zeros(fix(length(u(:,1))/L),8)*nan;
by=zeros(fix(length(u(:,1))/L),8)*nan;
theta=zeros(fix(length(u(:,1))/L),8)*nan;

urms=zeros(fix(length(u(:,1))/L),8)*nan;vrms=zeros(fix(length(u(:,1))/L),8)*nan;
wrms=zeros(fix(length(u(:,1))/L),8)*nan;thetarms=zeros(fix(length(u(:,1))/L),8)*nan;
ustar=zeros(fix(length(u(:,1))/L),8)*nan;

wndshear=zeros(fix(length(u(:,1))/L),8)*nan;
thetashearEC=zeros(fix(length(u(:,1))/L),8)*nan;thetasheartidu=zeros(fix(length(u(:,1))/L),8)*nan;
Rww=zeros(8,fix(length(u(:,1))/L),8)*nan;Ruu=zeros(L/5/fs,fix(length(u(:,1))/L),8)*nan;
Nbrunt=zeros(fix(length(u(:,1))/L),8)*nan;
Loz=zeros(fix(length(u(:,1))/L),8)*nan;
Lkn=zeros(fix(length(u(:,1))/L),8)*nan;
Lh=Loz*nan;Lv=Loz*nan;
Fh=zeros(fix(length(u(:,1))/L),8)*nan;Fv=zeros(fix(length(u(:,1))/L),8)*nan;
Rig=zeros(fix(length(u(:,1))/L),8)*nan;


for i=1:fix(length(u(:,1))/L)
    for iz=2:length(zcase99)
        bin=(i-1)*L+1:i*L;
        if length( find(isnan([u(bin,iz) v(bin,iz) w(bin,iz) T(bin,iz)])) )/L<0.1
            data=rotate1([u(bin,iz),v(bin,iz),w(bin,iz)]);
            theta(i,iz)=nanmean(T(bin,iz))+273.15;
            [z_L(i,iz),buoyancy(i,iz)]=zL(u(bin,iz),v(bin,iz),w(bin,iz),T(bin,iz),zcase99(iz));
            wnd(i,iz)=nanmean(   data(:,1)       );
            urms(i,iz)=nanstd(data(:,1));vrms(i,iz)=nanstd(data(:,2));
            wrms(i,iz)=nanstd(data(:,3));thetarms(i,iz)=nanstd(T(bin,iz));
            ustar(i,iz)=-nanmean(data(:,1).*data(:,3));

            Rauto=autocorr(data(:,1),L/5,fs);  %1
            xulie=find(Rauto<0.05);
            if isempty(xulie)
                Lh(i,iz)=nan;
            else
                indice=xulie(1);xulie=find(Rauto(1:indice)>0);indice=xulie(end);
                Lh(i,iz)=regress(wnd(i,iz)*(0:indice-1)',-log(Rauto(1:indice)'));
            end
            Ruu(:,i,iz)=Rauto;
            zbin=1:8;
            Rauto_z=autocorr_z(w(bin,:),iz,zbin,1*60*20);    %3
            xulie=find(Rauto_z>0.05);
            Lv(i,iz)=regress(abs(zcase99(xulie)'-zcase99(iz)),-log(Rauto_z(xulie)'));
            if Lv(i,iz)==0;Lv(i,iz)=nan;end
            clear Rauto_z Rauto
            Rww(:,i,iz)=Rauto_z;


            zbin=2:8;
            b=regress(nanmean((u(bin,zbin).^2+v(bin,zbin).^2).^0.5,1)',...
                [ones(length(zbin),1),log(zcase99(zbin)'),log(zcase99(zbin)').^2,log(zcase99(zbin)').^3]);
            wndshear(i,iz)=b(2)/zcase99(iz)+2*b(3)*log(zcase99(iz))/zcase99(iz)+...
                3*b(4)*(log(zcase99(iz)))^2/(zcase99(iz));
            zbin=2:8;
            b=regress(nanmean((T(bin,zbin))+273.15,1)',...
                [ones(length(zbin),1),log(zcase99(zbin)'),zcase99(zbin)',zcase99(zbin)'.^2]);
            thetashearEC(i,iz)=b(2)/zcase99(iz)+b(3)+2*b(4)*zcase99(iz);

            zbintidu=1:6;
            y=[nanmean((Ttidu(bintidu,zbintidu))+273.15,1)'];
            x=[zcase99tidu(zbintidu)'];
            b=regress(y,[ones(length(x),1),log(x),x,x.^2]);
            thetasheartidu(i,iz)=b(2)/zcase99(iz)+b(3)+2*b(4)*zcase99(iz);
            clear x y

            fstrc=structure_f(data(:,1),2,floor(zcase99(iz)/4/nanmean(data(:,1))*20));
            epi(i,iz)=regress(0.3535*fstrc'.^1.5,[((1:floor(zcase99(iz)/4/nanmean(data(:,1))*fs))/fs*nanmean(data(:,1)))']);
            clear fstrc

            Nbrunt(i,iz)=(9.8./theta(i,iz).*thetasheartidu(i,iz)).^0.5;
            Loz(i,iz)=(epi(i,iz)./Nbrunt(i,iz).^3).^0.5;
            Lkn(i,iz)=urms(i,iz)./Nbrunt(i,iz);
            Fh(i,iz)=urms(i,iz)./Nbrunt(i,iz)./Lh(i,iz);
            Fv(i,iz)=urms(i,iz)./Nbrunt(i,iz)./Lv(i,iz);
            Rig(i,iz)=Nbrunt(i,iz)^2/wndshear(i,iz)^2;


        end
    end
end
TKE=0.5*(urms.^2+vrms.^2+wrms.^2);


clear u;clear v;clear w;clear T Ttidu wndtidu
clear uf;clear vf;clear wf;clear data;










