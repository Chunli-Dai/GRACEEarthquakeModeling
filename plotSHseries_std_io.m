clear
clc
close all

%READ INPUT PARAMETERS
filename='plotSHseries_std_io.txt';
fid = fopen(filename);
[homedir, count]=fscanf(fid, '%s', [1, 1]);
[dir, count]=fscanf(fid, '%s', [1, 1]);
[filelist, count]=fscanf(fid, '%s', [1, 1]);
ymdse=fscanf(fid, '%f', [3, 2])'; %eg. 2004 1 1; 2013 12 31
stdflag= fscanf(fid, '%d', 1); %0, use the std in the data; 1, std=1;
numo = fscanf(fid, '%d', 1);
latlono=fscanf(fid, '%f', [2, numo])'; % 
fclose(fid);

% addpath(genpath('~/Osuwork/ChunliTools/Matlab_shorelines/'));
% setenv('IFILES','~/Osuwork/ChunliTools/Matlab_shorelines/');
%'F:\';%'~/';%'/Users/chunli/Desktop/tt/';
%addpath(genpath([macdir,'ChunliTools/Matlab_shorelines/']));
addpath(genpath([homedir,'/Osuwork/ChunliTools/Matlab_shorelines/']));
%dir=[hmdir,'../Software/Subtract_Reference_SHCs/CSRRL05mGOCOwstd/'];

compo={'C','S'}; 
ngra=2;

epochse=[datenum([ymdse(1,1) ymdse(1,2) ymdse(1,3)]),datenum([ymdse(2,1) ymdse(2,2) ymdse(2,3)])];

%fid=fopen('filelist_cs_CSRRL05Feb13','r');

% Generating eptick
%eptick=datenum({'1/1/03','1/1/04','1/1/05','1/1/06','1/1/07','1/1/08','1/1/09','1/1/10','1/1/11','1/1/12','1/1/13','1/1/14'});
nyrs=ymdse(2,1)-ymdse(1,1);
wof=8.; ibt=1.6; %Width of figure 8.0 inches; % inches between ticks >= 1.6
%nt<=8./ibt+1; % 8/(nt-1)>=ibt
maxnt=wof/ibt+1;
nt=maxnt;
dyr=int32((nyrs-1)/(nt-1.)); %nyrs-1=dyr*(nt-1)
if dyr < 1; dyr=1;end
nt=(nyrs-1)/dyr+1;
if nt > maxnt
    dyr=dyr+1;
    nt=(nyrs-1)/dyr+1;
end
fnt=floor(nt);
eptick=zeros(fnt,1);
for j=1:fnt
%    yeart=double(ymdse(1,1)+dyr*(j-1.));
    yeart=double(ymdse(2,1)+dyr*(j-fnt));
    eptick(j)=datenum([yeart 1 1]);
end 
% End of generating eptick

% Get data
fid=fopen(filelist,'r');
i=0;
while ~feof(fid)
 i = i+1;
% if i > lenf
%     break
% end
 tline1{i:i}=fgetl(fid);
%time(i) = (str2num(tline1{i}(6:8)))/365.25 + str2num(tline1{i}(3:4)) + 2000;   %yr
 epoch(i)=datenum(str2num(tline1{i}(1:4)),1, str2num(tline1{i}(6:8)));
 lmsh{i:i}=load([dir,tline1{i}]);
end
fclose(fid);
lenf=i;

% Delete some files
idd=[find(epoch < epochse(1)), find(epoch > epochse(2))];
if not(isempty(idd))
lmsh(idd)=[];epoch(idd)=[];tline1(idd)=[];
lenf=lenf-length(idd);
end
display(['number of files used: ',num2str(lenf)]);
% End of Delete some files

% epd=datenum({'2011/1/75'});
%lend=length(epd);
%idd=zeros(lend,1);
%
%for j=1:lend
%    idd(j)=find(epoch == epd(j));
%end
%lenf=lenf-lend;
%lmsh(idd)=[];epoch(idd)=[];tline1(idd)=[];

[nr,nc]=size(lmsh{1});
lat=zeros(nr,1);lon=zeros(nr,1);
stdcs=zeros(nr,ngra);
%jumpstd=zeros(nr,3);
% slope=zeros(nr,3);
% mean=  -0.0245    0.0633    0.1535 microGal/yr
% std(slope)=0.3571    1.3891    1.4359
% mean(abs(slope))=0.2853    1.0827    1.1223
% plot each grids 

%%%  Similar function as fit3.m
%%%  Fit the annual semiannual and 161 tidal S2 aliasing signals (Periods fixed)
%%%  fitting model: y=a + b(t-tm) + c1*cos(2pi(t-tm-t01)/T1)
%      + c2*cos(2pi(t-tm-t02)/T2)+ c3*cos(2pi(t-tm-t03)/T3)
%      +d11*(0 or 1);
%    9 parameters:  
%%%  Linear fitting model: y=a + b(t-tm) + A1*cos(2pi(t-tm)/T1)+B1*sin(2pi(t-tm)/T1)
%      + A2*cos(2pi(t-tm)/T2)+B2*sin(2pi(t-tm)/T2)
%      + A3*cos(2pi(t-tm)/T3)+ B3*sin(2pi(t-tm)/T3)
%      +d11*(0 or 1);
%    9 parameters:  
%    c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2pi);
%    c2=sqrt(A2^2+B2^2); t02=atan(B2/A2)*T2/(2pi);
%    c3=sqrt(A3^2+B3^2); t03=atan(B3/A3)*T3/(2pi);
%

% id11=find(epoch >= eqepoch); % Chile earthquake epoch

% id=find(epoch >= eqepoch); iddiag=1:lenf;
T6=zeros(lenf,1);T6std=zeros(lenf,1);
yr=365.25;
epochorg=epoch';
epoch=epochorg/yr;  % Change the unit of time from day to year, to avoid singular matrix
tm=mean(epoch);

% Find the observation point for plotting
lat=lmsh{1}(:,1);lon=lmsh{1}(:,2);

lono=latlono(:,2);lato=latlono(:,1);
dlat=1; dlon=1;

obsid=[];
for i=1:numo
    ido=find(abs(lmsh{1}(:,1)-lato(i))<=dlat/2.);
    idoo=find(abs(lmsh{1}(ido,2)-lono(i))<=dlon/2.);
    id=ido(idoo);
    obsid=[obsid,id];
end
numo=length(obsid);

for i=1:nr  %44937; %1:nr; 
%     i
%   lat(i:i)=lmsh{1}(i,1);lon(i:i)=lmsh{1}(i,2); %degree, order
    latt=lat(i);lont=lon(i);
    for k=1:ngra		% cnm snm
        for j=1:lenf
            T6(j)=lmsh{j}(i,3+k-1);
            if stdflag == 0
%     		T6std(j)=lmsh{j}(i,5+k-1);
                T6std(j)=lmsh{j}(i,3+ngra+k-1);
            elseif stdflag ==1
                T6std(j)=1.;
            end
        end
        P=diag(T6std(:).^-2); 

% Liner model
 %%%  Linear fitting model: y=a + b(t-tm) + A1*cos(2pi(t-tm)/T1)+B1*sin(2pi(t-tm)/T1)
%      + A2*cos(2pi(t-tm)/T2)+B2*sin(2pi(t-tm)/T2)
%      + A3*cos(2pi(t-tm)/T3)+ B3*sin(2pi(t-tm)/T3)
%      +d11*(0 or 1);
%    8 parameters:  
%    c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2pi);
%    c2=sqrt(A2^2+B2^2); t02=atan(B2/A2)*T2/(2pi);
%    c3=sqrt(A3^2+B3^2); t03=atan(B3/A3)*T3/(2pi);
        T1=365.25/yr;
        T2=T1/2.;
        T3=161./yr;
        epoch=epochorg/yr;
        clear AM
        AM=zeros(lenf,8);  %
        AM(:,1)=1.;AM(:,2)=epoch-tm;
        AM(:,3)=cos(2.*pi*(epoch-tm)/T1);
        AM(:,4)=sin(2.*pi*(epoch-tm)/T1);
        AM(:,5)=cos(2.*pi*(epoch-tm)/T2);
        AM(:,6)=sin(2.*pi*(epoch-tm)/T2);
        AM(:,7)=cos(2.*pi*(epoch-tm)/T3);
        AM(:,8)=sin(2.*pi*(epoch-tm)/T3);
%       AM(:,9)=0.; AM(id11,9)=1.;
        y=T6;

        var=inv(AM'*P*AM);
        
        for m=1:8 
            eststd(m:m)=sqrt(var(m,m)); 
        end
        est=var*AM'*P*y;
        
        A1=est(3);B1=est(4);A2=est(5);B2=est(6);A3=est(7);B3=est(8);
        c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2.*pi);
        c2=sqrt(A2^2+B2^2); t02=atan(B2/A2)*T2/(2.*pi);
        c3=sqrt(A3^2+B3^2); t03=atan(B3/A3)*T3/(2.*pi);
%         display([est(1:2)' c1 t01 c2 t02 c3 t03 est(9)])

%       jump(i,k)=est(end);jumpstd(i,k)=sqrt(var(end,end));	

        fit=est(1) + est(2)*(epoch-tm) +A1*cos(2.*pi*(epoch-tm)/T1)+B1*sin(2.*pi*(epoch-tm)/T1) + A2*cos(2.*pi*(epoch-tm)/T2)+B2*sin(2.*pi*(epoch-tm)/T2) ...
     + A3*cos(2.*pi*(epoch-tm)/T3)+ B3*sin(2.*pi*(epoch-tm)/T3);%+AM(:,9:end)*est(9:end);
%         
        
        etilde=y-AM*est;
        sigma02hat=etilde'*P*etilde/(lenf-m);
        stdcs(i,k)=sqrt(sigma02hat); %Special case for equal weight and identity dispersion matrix for y
        T6std(:)= sqrt(sigma02hat)*T6std(:);  %stdcs(i,k);
        var=sigma02hat*var;
       
%        jumpt=jump(i,k);jumpstdt=jumpstd(i,k);
%        slope(i,k)=est(2);
        csign=0;
        for j=1:numo
%           if i == obsid(j)  %Use the next, Nov 2014
            if i == obsid(j) && (lont ~= 0 || k~=2)
                csign=1;
            end
        end
        if csign==0; continue;end
 
        fitstdall=AM*var*AM';
        fitstd=zeros(lenf,1);
        for j=1:lenf
            fitstd(j)=sqrt(fitstdall(j,j));
        end

        epoch=epochorg;
 	formatSpec = '%d';
        figure
        set(gcf,'Color','white')
        set(gca,'FontSize', 12);
        set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
        hold all
        plot(epoch,T6(:),'b.-','linewidth',2,'markersize',18)
        plot(epoch,fit,'r.-','linewidth',2,'markersize',18)
        hold on
 	minmax=(abs(min(T6))+abs(max(T6)))/2.;
        shadedErrorBar(epoch,T6(:),T6std(:),'b',1)
        shadedErrorBar(epoch,fit,fitstd,'r',1)
        set(gca,'XTickLabel',[])
        set(gca,'XTick',eptick)
        set(gca,'XTickLabel',eptick,'FontSize',18,'XTickLabelMode','manual')
        datetick('x','mm/yy','keepticks')
        box on
        grid on
        title([compo{k}(:)',num2str(int32(latt)),',',num2str(int32(lont))],'FontSize',12)
%       title(['Grid, lat: ',num2str(latt),', lon:',num2str(lont),'; Component:',compo{k}(:)'])
	axis([min(epoch)-30 max(epoch)+30*4 min(T6)-minmax*0.2 max(T6)+minmax*0.3])
        ofile=['Degree',num2str(int32(latt)),'Order',num2str(int32(lont)),compo{k}(:)'];
        saveas(gcf,ofile,'fig')
%       print('-dpdf','-r500',ofile)
  	print('-dtiff','-r500',ofile)
   end
end

%outdir=[hmdir,'../Software/Subtract_Reference_SHCs/CSRRL05mGOCOwstdsigma0/'];
outdir=dir;
for j=1:lenf
   out=[ lmsh{j}(:,1:4)  stdcs(:,1:2)];
   ofile=[outdir,tline1{j}];
   fid = fopen(ofile,'w');
   fprintf(fid,'%5d%5d %23.15E %23.15E %23.15E %23.15E\n',out');  %2I5,4(1x,E23.15)
   fclose(fid);
% save -ascii jumpCSRN96ex1p3_FIT_wperiods_graNMAX96ncut70_grid_SM12.dat out
end
