clear
clc

%data=importdata('../../../Sumatra/SM04model/SM04earthmodel.dat');
%data=importdata('../../../Sumatra/SM05model/SM05earthmodelfmt.dat');
%data=importdata('../../../Sumatra/SM07model/SM07earthmodelfmt.dat');
%data=importdata('../../../Sumatra/SM12model/SM12earthmodelfmt.dat');
% data=importdata('SM07earthmodelfmt.dat');
filename='crust.dat';
fid = fopen(filename);
i=0;
while ~feof(fid)
 i = i+1;
 tline1=fgetl(fid);
end
frewind(fid);
lenf=i;

% header=fscanf(fid, '%s', 1);
header=fgetl(fid);
nl=lenf-2; %7;
data1=zeros(nl,4);data=zeros(nl+1,4);
for i=1:nl
    data1(i,:)=fscanf(fid, '%f', [4, 1])'; 
    fgetl(fid);
end
data2=fscanf(fid, '%f', [3, 1])'; 
data(1:nl,:)=data1;
data(nl+1,2:4)=data2;

% if abs(vs(1))<1e-9  
if abs(data(1,3))<1e-9  
 data(1,:)=[];
end
id=find(abs(data(:,1))<1e-9);
if id(end) == length(data(:,1)); id(end)=[];end
data(id,:)=[];

thi=data(:,1);
%vs=data.data(:,2);
%vp=data.data(:,3);
%rou=data.data(:,4)*1e3;
vp=data(:,2);
vs=data(:,3);
%rou=data.data(:,4);
rou=data(:,4)*1e3;

len=length(rou);
  j=1;i=1;
  dep(j)=0.;
  vp2(j)=vp(i);vs2(j)=vs(i);rou2(j)=rou(i);
  j=j+1;
  dep(j)=dep(j-1)+thi(i);
  vp2(j)=vp(i);vs2(j)=vs(i);rou2(j)=rou(i);
for i=2:len
  j=j+1;
  dep(j)=dep(j-1)+0.;
  vp2(j)=vp(i);vs2(j)=vs(i);rou2(j)=rou(i);
  j=j+1;
  dep(j)=dep(j-1)+thi(i);
  vp2(j)=vp(i);vs2(j)=vs(i);rou2(j)=rou(i);
end
len2=j-1;
dep=dep(1:len2);vp2=vp2(1:len2);vs2=vs2(1:len2);rou2=rou2(1:len2);
no=1:len2;
eta1=zeros(len2,1);
eta2=zeros(len2,1);
alpha=ones(len2,1);

out=[no',dep', vp2',vs2',rou2',eta1,eta2,alpha];
%# no  depth[km]  vp[km/s]  vs[km/s]  rho[kg/m^3] eta1[Pa*s] eta2[Pa*s] alpha
% 1     0.0      2.5000    1.2000    2100.0     0.0E+00    0.0E+00    1.000
% 2     0.2      2.5000    1.2000    2100.0     0.0E+00    0.0E+00    1.000
fid=fopen('tmp','w');
fprintf(fid,' %3d\n', len2);
fprintf(fid,' %3d   %7.4f   %7.4f   %7.4f   %8.4f   %1.1E   %1.1E   %1.1E\n', out');
fclose(fid);

%% Test
%A1 = [9.9, 9900];
%A2 = [8.8,  7.7 ; ...
%      8800, 7700];
%formatSpec = 'X is %4.2f meters or %8.3f mm\n';
%fprintf(formatSpec,A1,A2)
