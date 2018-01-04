clear
clc
close all

%READ INPUT PARAMETERS
filename='plotGRAseriesFIT_wperiods_io.txt';
fid = fopen(filename);
[datatype, count]=fscanf(fid, '%s', [1, 1]);
[gridtype, count]=fscanf(fid, '%s', [1, 1]);
[gra, count]=fscanf(fid, '%s', [1, 1]);
[homedir, count]=fscanf(fid, '%s', [1, 1]);
[dir, count]=fscanf(fid, '%s', [1, 1]);
[filelist, count]=fscanf(fid, '%s', [1, 1]);
ymdse=fscanf(fid, '%f', [3, 2])'; %eg. 2004 1 1; 2013 12 31
neq = fscanf(fid, '%d', 1);fgetl(fid);
% fscanf(fid, '[^\n]')
for j=1:neq
    texteq{j}=fgetl(fid);
end
ymdeq=fscanf(fid, '%f', [3, neq])'; 
latloneq=fscanf(fid, '%f', [2, neq])'; % 
lend = fscanf(fid, '%d', 1);
if lend == 0
 ymdd=[];
 fgetl(fid);
else
ymdd=fscanf(fid, '%f', [3, lend])'; % 2011 1 75
end
stdflag= fscanf(fid, '%d', 1); %0, use the std in the data; 1, std=1;
[linear, count]=fscanf(fid, '%s', [1, 1]);
numo = fscanf(fid, '%d', 1);
latlono=fscanf(fid, '%f', [2, numo])' % 
fclose(fid);

if strcmp(gra ,'gra') 
ngra=3;
gratext='Gravity ';
graunit='\muGal';
compo={'g_N','g_E','g_D'}; 
formatSpec = '%6.1f';
elseif strcmp(gra ,'grad')
ngra=6;
gratext='Gravity gradient ';
graunit='mE';
compo={'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
formatSpec = '%6.2f';
elseif strcmp(gra ,'SH')
ngra=2;
gratext='SH';
graunit='';
compo={'C','S'};
formatSpec = '%6.2e';
end

%addpath(genpath('F:\Osuwork/ChunliTools/Matlab_shorelines/'));
% addpath(genpath([homedir,'/Euler/Osuwork/ChunliTools/Matlab_shorelines/']));
addpath(genpath([homedir,'/Osuwork/ChunliTools/Matlab_shorelines/']));
%dir='~/share/CHL/progdcl/RL05CSRmEGM08GRAJP/'; 
%datatype='CSRRL05ex1';
%gridtype='grid_bw61';
% outdir='RL05GFZmGOCOGRAbw91FIT_wperiods_pic/';
%compo={'North','East','Down'};  %{'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
%lenf=106;
%neq=3;
epochse=[datenum([ymdse(1,1) ymdse(1,2) ymdse(1,3)]),datenum([ymdse(2,1) ymdse(2,2) ymdse(2,3)])];
eqepoch=zeros(neq,1);
for i=1:neq
eqepoch(i)=datenum([ymdeq(i,1) ymdeq(i,2) ymdeq(i,3)]);
end

% Generating eptick
% eptick=datenum({'1/1/04','1/1/05','1/1/06','1/1/07','1/1/08','1/1/09','1/1/10','1/1/11','1/1/12','1/1/13','1/1/14'});
% eptick=datenum({'1/1/03','1/1/04','1/1/05','1/1/06','1/1/07','1/1/08','1/1/09','1/1/10','1/1/11','1/1/12','1/1/13'});
% eptick=datenum({'1/1/05','1/1/07','1/1/09','1/1/11','1/1/13'});
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
%  if i > lenf
%      break
%  end
 tline1{i:i}=fgetl(fid);
% time(i) = (str2num(tline1{i}(6:8)))/365.25 + str2num(tline1{i}(3:4)) + 2000;   %yr
 epoch(i)=datenum(str2num(tline1{i}(1:4)),1, str2num(tline1{i}(6:8)));
 lmsh{i:i}=load([dir,tline1{i}]);
end
fclose(fid);
lenf=i;

% Delete some files
%epd=datenum({'2011/1/75'});
epd=zeros(lend,1);
for i=1:lend
epd(i)=datenum([ymdd(i,1) ymdd(i,2) ymdd(i,3)]);
end
%lend=length(epd);
idd=zeros(lend,1);
cnt=0;
for j=1:lend
%   idd(j)=find(epoch == epd(j))
%     id=find(abs (epoch - epd(j)) <= 20 );
    id=find(abs (epoch - epd(j)) <= 16 );
    if not(isempty(id))
        cnt=cnt+1;
        idd(cnt)=id;
    end
end
idd=idd(1:cnt);
if cnt >=1
lend=cnt;lenf=lenf-lend;
lmsh(idd)=[];epoch(idd)=[];tline1(idd)=[];
end
idd=[find(epoch < epochse(1)), find(epoch > epochse(2))];
if not(isempty(idd))
lmsh(idd)=[];epoch(idd)=[];tline1(idd)=[];
lenf=lenf-length(idd);
end
display(['number of files used: ',num2str(lenf)]);
% End of Delete some files


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

% id11=find(epoch >= eqepoch); % Japan earthquake epoch
% ideq=zeros(lenf,neq);
cnt=0;
ideqd=[];
for i=1:neq
    id=find(epoch >= eqepoch(i));
    if not(isempty(id)) && length(id)<lenf
        cnt=cnt+1;
        ideq{cnt}=id;
    else
        ideqd=[ideqd, i];
    end    
end
neq=cnt;
if not(isempty(ideqd))
    texteq(ideqd)=[]; eqepoch(ideqd)=[]; latloneq(ideqd,:)=[];
end
% texteq(1:neq)={'2004/12/26 Mw 9.1 2005/03/28 Mw 8.6', ...
%     '2007/09/12 Mw 8.5', ...
%     '2012/04/11 Mw 8.6+Mw8.2'};
[nr,nc]=size(lmsh{1});
jump=zeros(nr,ngra,neq);jumpstd=zeros(nr,ngra,neq);
slope=zeros(nr,ngra,neq);

phis(1:neq)=latloneq(:,2);
lats(1:neq)=latloneq(:,1);
thetas=90.-lats;

T6=zeros(lenf,1);T6std=zeros(lenf,1);
yr=365.25;
epochorg=epoch';
epoch=epochorg/yr;  % Change the unit of time from day to year, to avoid singular matrix
tm=mean(epoch);
%lon=141.2;lat=36;i=3589;
%lon=139.6;lat=36.4;i=3484
% lono=299.8; lato=-25.75; %i=; 
lat=lmsh{1}(:,1);lon=lmsh{1}(:,2);
id=find(lat==lat(1));
nlon=length(id); nlat=nr./nlon;

lono=latlono(:,2);lato=latlono(:,1);
dlat=(max(lat)-min(lat))./(nlat-1.); dlon=(max(lon)-min(lon))./(nlon-1.);

if strcmp(gra ,'SH'); dlat=1;dlon=1; end %hi
obsid=[];
for i=1:numo
%     ido=find(abs(lmsh{1}(:,1)-lato(i))<=dlat/2.); % a -dlat/2. = 1e-16
%     idoo=find(abs(lmsh{1}(ido,2)-lono(i))<=dlon/2.);
    ido=find(abs(lmsh{1}(:,1)-lato(i))- dlat/2. < 1e-14 );
    idoo=find(abs(lmsh{1}(ido,2)-lono(i)) - dlon/2. < 1e-14);
    id=ido(idoo);
%     obsid=[obsid,id]; %CAT arguments dimensions are not consistent.
    obsid=[obsid;id]; 
end
numo=length(obsid) 

mp=8+neq;

for i=1:nr %3484 %3461:3461  %2571:2571 %1:nr %2767:2767  % EPcenter 3086:3086
%   i
%     lat(i:i)=lmsh{1}(i,1);lon(i:i)=lmsh{1}(i,2);
    latt=lat(i);lont=lon(i);
    for k=1:ngra		% North East Down
        for j=1:lenf
            T6(j)=lmsh{j}(i,3+k-1);
            if stdflag == 0 
                T6std(j)=lmsh{j}(i,3+ngra+k-1);
            elseif stdflag ==1
                T6std(j)=1.;
            end
        end
        
        t1=max(abs(T6));t2=max(abs(T6std));
        if t1 < 1e-16 && t2 < 1e-16 %If e.g. S20=0 std=0
           continue  %skip remaining
        end
        
        P=diag(T6std(:).^-2); 

	% Non-linear fitting
	if not(linear)
        estpre=0.;
        t01=0.;T1=365.25/yr;c1=0.2;
        t02=0.;T2=T1/2.;c2=0.2;
        t03=0.;T3=161./yr;c3=0.2;
        count=0;
        epoch=epochorg/yr;
        while true
        clear AM
        AM=zeros(lenf,mp);  %
        AM(:,1)=1.;AM(:,2)=epoch-tm;
        AM(:,3)=cos(2.*pi*(epoch-tm-t01)/T1);
        AM(:,4)=c1*sin(2.*pi*(epoch-tm-t01)/T1)*2.*pi/T1;
        AM(:,5)=cos(2.*pi*(epoch-tm-t02)/T2);
        AM(:,6)=c2*sin(2.*pi*(epoch-tm-t02)/T2)*2.*pi/T2;
        AM(:,7)=cos(2.*pi*(epoch-tm-t03)/T3);
        AM(:,8)=c3*sin(2.*pi*(epoch-tm-t03)/T3)*2.*pi/T3;
        for j=1:neq
        AM(:,8+j)=0.; AM(ideq{j},8+j)=1.;
        end
        y=T6-c1*cos(2*pi*(epoch-tm-t01)/T1)-c2*cos(2*pi*(epoch-tm-t02)/T2)-c3*cos(2*pi*(epoch-tm-t03)/T3);
        var=inv(AM'*P*AM);
        
        for m=1:mp 
            eststd(m:m)=sqrt(var(m,m)); 
        end
        est=var*AM'*P*y;
%       est=AM\y; %No weight

        count=count+1;
        if max(abs(est-estpre)) < 1e-6
            break;
        end
        c1=c1+est(3);
        t01=t01+est(4);
        c2=c2+est(5);
        t02=t02+est(6);
        c3=c3+est(7);
        t03=t03+est(8);
        estpre=est;
        end
%         jump(i,k)=est(end);jumpstd(i,k)=sqrt(var(end,end));	
        jump(i,k,1:neq)=est(9:end);
        for k1=1:neq
            jumpstd(i,k,k1)=sqrt(var(k1+8,k1+8));	
        end

        fit=est(1) + est(2)*(epoch-tm) +c1*cos(2*pi*(epoch-tm-t01)/T1)+c2*cos(2*pi*(epoch-tm-t02)/T2)+c3*cos(2*pi*(epoch-tm-t03)/T3)+AM(:,9:end)*est(9:end);

% Liner model
 %%%  Linear fitting model: y=a + b(t-tm) + A1*cos(2pi(t-tm)/T1)+B1*sin(2pi(t-tm)/T1)
%      + A2*cos(2pi(t-tm)/T2)+B2*sin(2pi(t-tm)/T2)
%      + A3*cos(2pi(t-tm)/T3)+ B3*sin(2pi(t-tm)/T3)
%      +d11*(0 or 1);
%    9 parameters:  
%    c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2pi);
%    c2=sqrt(A2^2+B2^2); t02=atan(B2/A2)*T2/(2pi);
%    c3=sqrt(A3^2+B3^2); t03=atan(B3/A3)*T3/(2pi);
	elseif linear
        if strcmp(datatype , 'OSUyr')
 %%%  Linear fitting model: y=a + b(t-tm) 
%      +d11*(0 or 1);
%    3 parameters:  
%    c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2pi);
        mpnoj=2;
        mp=mpnoj+neq;
        epoch=epochorg/yr;
        clear AM
        AM=zeros(lenf,mp);  %
        AM(:,1)=1.;AM(:,2)=epoch-tm;
%         AM(:,9)=0.; AM(id11,9)=1.;
        for j=1:neq
        AM(:,mpnoj+j)=0.; AM(ideq{j},mpnoj+j)=1.;
        end

        y=T6;

        cdA=cond(AM);
        if(cdA>1e5); display(['condtion number of AM is large: ',num2str(cdA)]);end
        
        var=inv(AM'*P*AM);
        
        for m=1:mp 
            eststd(m:m)=sqrt(var(m,m)); 
        end
        est=var*AM'*P*y;

%         jump(i,k)=est(end);jumpstd(i,k)=sqrt(var(end,end));	
        jump(i,k,1:neq)=est(mpnoj+1:end);
        for k1=1:neq
            jumpstd(i,k,k1)=sqrt(var(k1+mpnoj,k1+mpnoj));	
        end

        fit=est(1) + est(2)*(epoch-tm) ...
     +AM(:,mpnoj+1:end)*est(mpnoj+1:end);
 
        else  % Default choice 
        T1=365.25/yr;
        T2=T1/2.;
        T3=161./yr;
        epoch=epochorg/yr;
        clear AM
        AM=zeros(lenf,mp);  %
        AM(:,1)=1.;AM(:,2)=epoch-tm;
        AM(:,3)=cos(2.*pi*(epoch-tm)/T1);
        AM(:,4)=sin(2.*pi*(epoch-tm)/T1);
        AM(:,5)=cos(2.*pi*(epoch-tm)/T2);
        AM(:,6)=sin(2.*pi*(epoch-tm)/T2);
        AM(:,7)=cos(2.*pi*(epoch-tm)/T3);
        AM(:,8)=sin(2.*pi*(epoch-tm)/T3);
%         AM(:,9)=0.; AM(id11,9)=1.;
        for j=1:neq
        AM(:,8+j)=0.; AM(ideq{j},8+j)=1.;
        end

        y=T6;

        cdA=cond(AM);
        if(cdA>1e10); display(['condtion number of AM is large: ',num2str(cdA)]);end
        
        var=inv(AM'*P*AM);
        
        for m=1:mp 
            eststd(m:m)=sqrt(var(m,m)); 
        end
        est=var*AM'*P*y;
        A1=est(3);B1=est(4);A2=est(5);B2=est(6);A3=est(7);B3=est(8);
        c1=sqrt(A1^2+B1^2); t01=atan(B1/A1)*T1/(2.*pi);
        c2=sqrt(A2^2+B2^2); t02=atan(B2/A2)*T2/(2.*pi);
        c3=sqrt(A3^2+B3^2); t03=atan(B3/A3)*T3/(2.*pi);

        if stdflag == 1 %estimate the variance from postfit residuals. %referes to plotSHseries_std_io.m
        etilde=y-AM*est;
        sigma02hat=etilde'*P*etilde/(lenf-m);
%         stdcs(i,k)=sqrt(sigma02hat); %Special case for equal weight and identity dispersion matrix for y
        T6std(:)= sqrt(sigma02hat)*T6std(:);  %stdcs(i,k);
        var=sigma02hat*var;
        end
        
%         jump(i,k)=est(end);jumpstd(i,k)=sqrt(var(end,end));	
        jump(i,k,1:neq)=est(9:end);
        for k1=1:neq
            jumpstd(i,k,k1)=sqrt(var(k1+8,k1+8));	
        end

        fit=est(1) + est(2)*(epoch-tm) +A1*cos(2.*pi*(epoch-tm)/T1)+B1*sin(2.*pi*(epoch-tm)/T1) + A2*cos(2.*pi*(epoch-tm)/T2)+B2*sin(2.*pi*(epoch-tm)/T2) ...
     + A3*cos(2.*pi*(epoch-tm)/T3)+ B3*sin(2.*pi*(epoch-tm)/T3)+AM(:,9:end)*est(9:end);
        end
	end
 %         
 
        csign=0;
        for j=1:numo
            if i == obsid(j)
                csign=1;
            end
        end
        if csign==0; continue;end
        
% % For plotting figures       
%         jumpt=jump(i,k);jumpstdt=jumpstd(i,k);
%         slope(i,k)=est(2);
        epoch=epochorg;
        epochfit=epoch;
        lenfit=200;
        if ( lenf < 100 && strcmp(linear,'true') ) %establish a denser points.
            clear AM
            epochfit=linspace(min(epoch),max(epoch),lenfit);
            AM=zeros(lenfit,mp);  %
            AM(:,1)=1.;AM(:,2)=epochfit/yr-tm;
            AM(:,3)=cos(2.*pi*(epochfit/yr-tm)/T1);
            AM(:,4)=sin(2.*pi*(epochfit/yr-tm)/T1);
            AM(:,5)=cos(2.*pi*(epochfit/yr-tm)/T2);
            AM(:,6)=sin(2.*pi*(epochfit/yr-tm)/T2);
            AM(:,7)=cos(2.*pi*(epochfit/yr-tm)/T3);
            AM(:,8)=sin(2.*pi*(epochfit/yr-tm)/T3);
    %         AM(:,9)=0.; AM(id11,9)=1.;
            for ieq=1:neq
                id=find(epochfit >= eqepoch(ieq));
                if not(isempty(id)) && length(id)<lenfit
                    ideqfit{ieq}=id;
                end    
            end
            for j=1:neq
            AM(:,8+j)=0.; AM(ideqfit{j},8+j)=1.;
            end
            fit=AM*est;
        end
        fitstdall=AM*var*AM';
        fitstd=zeros(length(fitstdall(:,1)),1);
        for j=1:length(fitstdall(:,1))
            fitstd(j)=sqrt(fitstdall(j,j));
        end
        
%         formatSpec = '%6.1f';
        figure  
        set(gcf,'Color','white')
        set(gca,'FontSize', 12);
        set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
        hold all
        plot(epoch,T6(:),'b.-','linewidth',2,'markersize',18)
        plot(epochfit,fit,'r.-','linewidth',2,'markersize',18)
        hold on
%         if stdflag ~= 1
        shadedErrorBar(epoch,T6(:),T6std(:),'b',1)
        shadedErrorBar(epochfit,fit,fitstd,'r',1)
%         end

%        minmax=(abs(min(T6))+abs(max(T6)))/2.;
        minmax=max(abs(T6));
        meanT=mean(T6);dev=max(T6)-min(T6);
        uplmt=ceil(max(T6)+dev*0.3);lowerlmt=floor(min(T6)-dev*0.2);
        for k1=1:neq
            jumpt=jump(i,k,k1);jumpstdt=jumpstd(i,k,k1);
%             plot([eqepoch(k1),eqepoch(k1)],[min(T6)*0.8, max(T6)*0.8],'k')
            plot([eqepoch(k1),eqepoch(k1)],[lowerlmt uplmt],'k','Linewidth',2)
            text(eqepoch(k1)+15,max(T6),texteq(k1),'FontSize',10)
            text(eqepoch(k1)+15,max(T6)+dev*0.2 ,['Jump=',num2str(jumpt,formatSpec),'\pm',num2str(jumpstdt,formatSpec),graunit],'FontSize',10)
        end
%         ylabel([gratext,compo{k}(:)',' (',graunit,')'],'FontSize',12)
        ylabel([compo{k}(:)',' (',graunit,')'],'FontSize',12)
        set(gca,'XTickLabel',[])
        set(gca,'XTick',[])
        set(gca,'XTick',eptick)
        set(gca,'XTickLabel',eptick,'FontSize',12,'XTickLabelMode','manual')
        datetick('x','mm/yy','keepticks')
        box on
        grid on
        if (uplmt> lowerlmt)
        axis([min(epoch)-30 max(epoch)+30*4 lowerlmt uplmt])
        end

	    ofile=['Lat',num2str(int32(latt*10)),'Lon',num2str(int32(lont*10)),compo{k}(:)'];
        saveas(gcf,ofile,'fig')
       	print('-dpdf','-r500',ofile)
        print('-dtiff','-r500',ofile) %better for older matlab before 2016, but larger size
    end
end


for k1=1:neq
% out=[lat lon jump jumpstd];
out=[lat lon jump(:,:,k1) jumpstd(:,:,k1)];
% save -ascii jumpCSRRL05ex1_FIT_wperiods_graNMAX60_grid_JPFeb13.dat out
% save -ascii jumpCSRRL05ex1_FIT_wperiods_graNMAX60_grid_JPApr12.dat out
% save -ascii jumpGFZex1_FIT_wperiods_graNcut60_grid_JP.dat out
% save -ascii jumpGFZex1_FIT_wperiods_graNcut60_grid_bw61.dat out
% save -ascii jumpCSRRL05ex1Jun_FIT_wperiods_graNMAX60_grid_JPFeb13.dat out
% save -ascii jumpCSRRL05ex1_FIT_wperiods_grastdNMAX60_grid_JPFeb13.dat out
% save -ascii jumpCSRRL05ex1_FIT_wperiods_grastdNMAX60_grid_bw61Feb13.dat out
% save -ascii jumpCSRRL05mEGM08ex1_FIT_wperiods_graNMAX60_grid_JPFeb13.dat out
ofile=['jump',datatype,'_FIT_wperiods_',gra,'_',gridtype,texteq{k1}(1:4),'.dat'];
if strcmp(gra ,'SH')
    fid = fopen(ofile, 'w');
    fprintf(fid, '%5d%5d %23.15e %23.15e %23.15e %23.15e\n', out');
    fclose(fid);
else
    save(ofile,'out','-ascii')
end
% save('ofilename','ofile','-ascii')
end


 % nlon=122;nlat=122;
%  nlon=101;nlat=64;
id=find(lat==lat(1));
nlon=length(id); nlat=nr./nlon;
 % nlon=122;nlat=122;
 % nlon=242;nlat=242;
 % nlon=182;nlat=182;
 
 LON=reshape(lon,nlon,nlat)';
 LAT=reshape(lat,nlon,nlat)';
 TH = 20;
%  homedir='/Users/chunli/';
 addpath(genpath([homedir,'/share/Osuworks/UltraSH/Slepian/SimonsSoftware']));
 setenv('IFILES',[homedir,'/share/Osuworks/UltraSH/Slepian/IFILES']);
 
for k1=1:neq
figure  
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
nrsub=ngra/3.;
 for k=1:ngra
 JUMP=reshape(jump(:,k,k1),nlon,nlat)';
 cm=int32(max(max(abs(JUMP)))*0.8);
 subplot(nrsub,3,k)
 pcolor(LON,LAT,JUMP)
 ylabel('Latitude');xlabel('Longitude')
 colormap(jet);
 shading interp
 colorbar
%  title([gratext,compo{k}(:)',' ',texteq{k1}])
 title([texteq{k1},', ',compo{k}(:)'])
 hold on
 [ax1,ch,XY1]=plotcont([],[],1);
 [ph,XY2]=plotplates([],[],1);
 hold on 
 circ(TH,2*pi,[phis(k1) 90-thetas(k1)]); 
 caxis([-cm cm]);
 axis([min(lon) max(lon) min(lat) max(lat)]);
 end
ofile=[gra,'jump_',datatype,'_',gridtype,texteq{k1}(1:4)];
saveas(gcf,ofile,'fig')
% print('-dtiff','-r50',ofile)
% print('-dpdf','-r500',ofile)        %Segmentation fault 
 
figure  
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
 for k=1:ngra
 JUMPSTD=reshape(jumpstd(:,k,k1),nlon,nlat)';
 cm=int32(max(max(abs(JUMPSTD)))*0.8);
subplot(nrsub,3,k)
 pcolor(LON,LAT,JUMPSTD)
 colormap(jet);shading interp
 colorbar
 caxis([-cm cm]);
 hold on 
 circ(TH,2*pi,[phis(k1) 90-thetas(k1)]); 
 ylabel('Latitude');xlabel('Longitude')
%  title([gratext,'std ',compo{k}(:)',' ',texteq{k1}])
 title([texteq{k1},', ',compo{k}(:)',' std'])
 axis([min(lon) max(lon) min(lat) max(lat)]);
 end
ofile=[gra,'jumpstd_',datatype,'_',gridtype,texteq{k1}(1:4)];
saveas(gcf,ofile,'fig')
% print('-dpdf','-r500',ofile)
end
