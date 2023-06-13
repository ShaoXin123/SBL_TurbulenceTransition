function data=rotate1(Rawdata)
datamean=zeros(2,0);datamean=nanmean(Rawdata);
U=datamean(1);V=datamean(2);W=datamean(3);
 
YawMatrix=[U/(sqrt(U^2+V^2)),-V/(sqrt(U^2+V^2)),0;V/(sqrt(U^2+V^2)),U/(sqrt(U^2+V^2)),0;0,0,1];
YawData=Rawdata(:,1:3)*YawMatrix;
datamean=nanmean(YawData);
U=datamean(1);V=datamean(2);W=datamean(3);

 
% PitchMatrix=[U/(sqrt(U^2+W^2)),0,-W/(sqrt(U^2+W^2)); 0,1,0;W/(sqrt(U^2+W^2)),0,U/(sqrt(U^2+W^2));];
% PitchData=YawData*PitchMatrix;
% datamean=zeros(2,0);datamean=nanmean(PitchData);
% datacov=zeros(2,0);datacov=cov(PitchData);
 
% Beta=0.5*atan(2*datacov(2,3)/(-datacov(2,2)+datacov(3,3)));
% RollMatrix=[1,0,0;0,cos(Beta),sin(Beta);0,-sin(Beta),cos(Beta)];
% RollData=PitchData*RollMatrix;

data=Rawdata;
data(:,1:3)=YawData;