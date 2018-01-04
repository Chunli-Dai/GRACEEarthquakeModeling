clear all;
clc
close('all');

%READ INPUT PARAMETERS
filename='Grid2PSD_io.txt';
fid = fopen(filename);
[datatype, count]=fscanf(fid, '%s', [1, 1]);
[homedir, count]=fscanf(fid, '%s', [1, 1]);  fgetl(fid);
for j=1:2
    ifile{j}=fgetl(fid);
end
[texteq, count]=fscanf(fid, '%s', [1, 1]);
latloneq=fscanf(fid, '%f', [2, 1])';
bw = fscanf(fid, '%d', 1)
fclose(fid);

phi = double(int32(latloneq(2))); latc = double(int32(latloneq(1)));  %Japan

% dir='/Users/chunli/Euler/';
% dir2='/Users/chunli/';%'/0/home/dai.56/';
addpath(genpath([homedir,'/share/Osuworks/UltraSH/Slepian/SimonsSoftware']));
setenv('IFILES',[homedir,'/share/Osuworks/UltraSH/Slepian/IFILES']);


TH = 20; L = 15;  omega = 0; 
% phi = 200; latc = 0;
% phi = 15; latc = -10;
% phi = 270; latc = 75;
%phi =140; latc =40;
%phi = 143; latc = 37;  %Japan
% phi = 287; latc = -36; %Chile
% phi = 94; latc = 2;  %Sumatra 2004 2005 2012
% phi = 101; latc = -4;  %Sumatra  2007

omega =110;  %20; %110; % Lei Wang's dissertation % Use integers to be save to IFILES/GLMALPHAPTO/
windcenter=['_Lon',num2str(int32(phi),'%3d'),'Lat',num2str(int32(latc),'%2d')];
theta=90-latc;
% bw=96+1;
%  bw=70+1;
%  bw=60+1;
%bw=90+1;
%bw=120+1;
% bw=900;
dlat=180./2./bw;dlon=360./2./bw;
lat_start=90-0.5*dlat;
lon_start=0;
latwi=lat_start:-dlat:-lat_start;
lonwi=lon_start:dlon:360-dlon;
ilat=int32((90-theta-lat_start)/(-dlat))+1;
jlon=int32((phi-lon_start)/dlon)+1;

outdir='./';

Comtopog{1}=load([ifile{1}]);
datatypeg{1}=[datatype,'_Gra_GRACE'];
compog{1}={'North','East','Down'};  %{'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
ntypeg{1}=3;

Comtopog{2}=load([ifile{2}]);
datatypeg{2}=[datatype,'_Grad_GRACE'];
compog{2}={'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
ntypeg{2}=6;
% 
Comtopog{3}=load(['SA_graD_sigma_flat_bw900_potential_gra.txt']);
datatypeg{3}=[datatype,'_Gra_model'];
compog{3}={'North','East','Down'};
ntypeg{3}=3;

Comtopog{4}=load(['SA_graD_sigma_flat_bw900_potential_grad.txt']);
datatypeg{4}=[datatype,'_Grad_model'];
compog{4}={'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
ntypeg{4}=6;


% Comtopog{3}=load([dir,'Osuwork/GRACE/Sumatra/SM05model/SA_graD_sigma_flat_bw900_potential_graN90bw91_SM0405.dat']);
% datatypeg{3}='N90_SM04_Gra_';
% Comtopog{3}=load([dir,'Osuwork/GRACE/Sumatra/SM12model/SA_graD_sigma_flat_bw900_potential_graN90bw91.dat']);
% datatypeg{3}='N90_SM12_Gra_';
% Comtopog{3}=load([dir,'Osuwork/GRACE/Sumatra/SM07model/SA_graD_sigma_flat_bw900_potential_graN90bw91.dat']);
% datatypeg{3}='N90_SM07_Gra_';
% outdir='MN96GRAChile/';

% Comtopog{4}=load([dir,'Osuwork/GRACE/Seismic2gravityByLei/uf/coseism_slipWei_po_crustWei_bw900/SA_graD_sigma_flat_bw900_potential_grad_gdbw900_orb.txt']);
% datatypeg{4}='coseism_slipWei_po_crustWei_grad_bw900_orb';
% compog{4}={'Txx','Tyy','Tzz','Txy','Txz','Tyz'};
% ntypeg{4}=6;

for idt=1:4 %[1,3] %idt=1;3 %1:4 %id of data type
ntype=ntypeg{1,idt};

% For case of regional input data, padding zeros to be global
% Refers to 2006Nov2007JanKurilIslandEQ/paddingzeros2gauss.m
%  lat=(90-0.2):-0.4:(-90.+0.2);
%  lon=0:0.4:(360.-0.4);
 [nr,nc]=size(Comtopog{idt});
 [LON,LAT]=meshgrid(lonwi,latwi);
 LATT=LAT';LONT=LON';[m,n]=size(LAT);zbg=zeros(m*n,ntype);
 gridbg=[LATT(:), LONT(:), zbg];
 len=nr;
 for i=1:len
%      id=find(abs(gridbg(:,1)-Comtopo(i,1))<1e-9); %Slow
%      id2=find(abs(gridbg(id,2)-Comtopo(i,2))<1e-9);
%      idbg=id(1)-1+id2; 
     lati=Comtopog{idt}(i,1);loni=Comtopog{idt}(i,2);
     ilat=int32((lati-lat_start)/(-dlat))+1;
     jlon=int32((loni-lon_start)/dlon)+1;
     idbg=(ilat-1)*2*bw+jlon;
     gridbg(idbg,2+(1:ntype))=Comtopog{idt}(i,2+(1:ntype));
 end
 Comtopog{idt}=gridbg;
 % End of padding zeros


Comtopo=Comtopog{1,idt}(:,:);
datatypet=datatypeg{1,idt};
compo=compog{1,idt};

lat=Comtopo(:,1)';lon=Comtopo(:,2)';

for itype=1:ntype  %[1,4,5] %1:ntype
% for itype=1:6

data = reshape(Comtopo(:,2+itype),bw*2,bw*2);
data = data';
% clear Comtopo

% figure, imagef([0 lat_start],[360-dlon -lat_start],data);
figure, imagef([lon(1) lat(1)],[lon(end) lat(end)],data);
 minmax=mean([abs(min(min((data)))),max(max(data))]);
 caxis([-minmax minmax]);
colorbar;grid on;hold on;axis equal;axis([0 360 -90 90]);
% % change the circ function, using caploc.m, dcl, Sep 5, 2012
circ(TH,2*pi,[phi 90-theta]);
% load coast
% id=find(long<0);
% long(id)=long(id)+360.;
% plot(long,lat,'b.')
load('topo');
contour(0:359,-89:90,topo,[0 0],'k')
hold off;
title([datatypet,compo{itype}(:)'])
saveas(gcf,[outdir,datatypet,'window',windcenter],'tiff');
% [m,n]=size(data);
% data=ones(m,n);
% 

[G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalphapto(TH,L,phi,theta,omega);
display(['Shannon number:', num2str(N)])
[dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=addmon(L);
% %G1 = G(:,1);
[Y,ORD] = sort(V,'descend');
% save ([outdir,datatypet,compo{itype}(:)','Eigenvalues',windcenter], 'Y',  '-ascii') 

sdl_iwe=zeros(bw,1); weight=0.;
 
 for rank = 1:floor(N)%3 %7
% %G1(:,k) = G(:,rank);
 %figure,plotslep(G,ORD(rank),1);
 coe = lmcosi(:,3:4);
 coe(ronm) = G(:,ORD(rank));
 lmcosi1 = [lmcosi(:,1:2) coe];
% [window1,lonw,latw,Plm]=plm2xyz(lmcosi1,lat',lon',L);  % in this case, lack of memory.
% [window1,lonw,latw,Plm]=plm2xyz(lmcosi1,lat',lon',L); % This is slow (20 minutes).
% [window2,lonw,latw,Plm]=plm2xyz(lmcosi1,dlat,[0 lat_start 360-dlat -lat_start],L); % this is 2 minutes
[window2,lonw,latw,Plm]=plm2xyz(lmcosi1,dlat,[lon(1) lat(1) lon(end) lat(end)],L); % this is 2 minutes
window2r=window2(:,1:2:end);
clear Comtopo_w
 Comtopo_w(:,:) = data .* window2r;
 clear window2r window2
% outname=[outdir,datatypet,compo{itype}(:)','w',num2str(rank,'%1d'),'_LatLonD_L15',windcenter];
 
% figure, imagef([0 lat_start],[360-dlon -lat_start],Comtopo_w(:,:));%caxis([-2 0.5]);
% colorbar;grid on;hold on;axis equal;axis([phi-20 phi+20 latc-20 latc+20]);
% minmax=mean([abs(min(min((Comtopo_w)))),max(max(Comtopo_w))]);
% caxis([-minmax minmax]);
% circ(TH,2*pi,[phi 90-theta]);
% load('topo');
% contour(0:359,-89:90,topo,[0 0],'k');hold off;
% title(outname)
%  saveas(gcf,[outname,'_global'],'jpg');
 
% [lmcosi_inv]=xyz2plm(Comtopo_w',bw-1,'im',lat,lon); %Slow, accurate; %2m40s for bw=61
clock
Comtopo_w=[Comtopo_w,Comtopo_w(:,1)];
[lmcosi_inv]=xyz2plm(Comtopo_w,bw-1,'im'); %Fast, but the Gauss grid is taken as [[0 360 -90 90]]. The difference is 3%
        % 36minutes for bw900
clock
[sdl_i,l]=plm2spec(lmcosi_inv,3); % 2 for PSD; 3 for Degree variances
sdl_iwe=sdl_iwe+sdl_i*V(rank);
weight=weight+V(rank);

% Comtopo_w1 = Comtopo_w(:,:)';
% clear Comtopo_w
% Comtopo_w1_LatLonD = [lat' lon' Comtopo_w1(:)]; 
%  save (outname, 'Comtopo_w1_LatLonD',  '-ascii')
%save 'Comtopo_w1_LatLonD_L15' Comtopo_w1_LatLonD  -ascii

  
 end %rank
 sdl_iwe=sdl_iwe/weight;
%  grapsd{idt}(:,itype)=sdl_iwe/sum(sdl_iwe)*100.;
 grapsd{idt}(:,itype)=sdl_iwe;
 grapsdpercpl{idt}(:,itype)=sdl_iwe/sum(sdl_iwe)*100.;
 nl=length(l);
 for ip=1:nl
 grapsdperc{idt}(ip,itype)=sum(sdl_iwe(1:ip))/sum(sdl_iwe)*100.;
 end
end % itype
end % idt=1:4 %id of data type

Darkbrown=[0.40,0.26,0.13];

% Compare with PSD using Dr. Guo's Grid2SHCsLC.exe: coseism_slipWei_po_crustWei_gra_bw61CSRRL05N60_Gra_degree_variance_Lon143Lat37_pub.pdf 
% The difference doesn't change the conclusion.
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
%set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
set(gcf, 'PaperPosition', [0 0 6 4.5],'Papersize',[6 4.5]);
hold all
plot(l,grapsdpercpl{3}(1:end,1),'b--','linewidth',2)
plot(l,grapsdpercpl{1}(1:end,1),'b-','linewidth',2)
plot(l,grapsdpercpl{4}(1:end,1),'--','Color',Darkbrown,'linewidth',2)
plot(l,grapsdpercpl{2}(1:end,1),'-','Color',Darkbrown,'linewidth',2)
plot(l,grapsdpercpl{4}(1:end,4),'g--','linewidth',2)
plot(l,grapsdpercpl{2}(1:end,4),'g-','linewidth',2)
plot(l,grapsdpercpl{4}(1:end,5),'r--','linewidth',2)
plot(l,grapsdpercpl{2}(1:end,5),'r-','linewidth',2)
legend('Predicted g_N','GRACE observed g_N','Predicted Txx','GRACE observed Txx ','Predicted Txy','GRACE observed Txy','Predicted Txz','GRACE observed Txz','location','NorthWest')
xlabel('Spherical Harmonic Degree');
ylabel('Degree Variance Percentage');
box on
grid on
ofile=['PSD_',datatype,'GraGrad','_',texteq];
saveas(gcf,ofile,'fig');
print('-dpdf','-r500',ofile)

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
%set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
set(gcf, 'PaperPosition', [0 0 6 4.5],'Papersize',[6 4.5]);
hold all
plot(l,grapsd{3}(1:end,1),'b--','linewidth',2)
plot(l,grapsd{1}(1:end,1),'b-','linewidth',2)
plot(l,grapsd{3}(1:end,2),'g--','linewidth',2)
plot(l,grapsd{1}(1:end,2),'g-','linewidth',2)
plot(l,grapsd{3}(1:end,3),'r--','linewidth',2)
plot(l,grapsd{1}(1:end,3),'r-','linewidth',2)
% legend('Predicted g_N','GRACE observed g_N','Predicted g_E','GRACE observed g_E ','Predicted g_D','GRACE observed g_D','Location','NorthWest')% legend('location','NorthEast')
legend('Predicted g_N','GRACE g_N','Predicted g_E','GRACE g_E ','Predicted g_D','GRACE g_D','Location','NorthWest')% legend('location','NorthEast')
xlabel('Spherical Harmonic Degree');
ylabel('Degree Variance \muGal^2');
legend('location','NorthWest')
box on
grid on
% axis([0 60 0 1])
% axis([0 100 0 12])
% axis([0 100 0 0.4])
% axis([min(l) max(l)*1.2 0 max()])
% ofile='PSD_JPcopo_Gra_bw900_R0';
ofile=['PSD_',datatype,'gra','_',texteq];
saveas(gcf,ofile,'fig');
print('-dpdf','-r500',ofile)

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
%set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
set(gcf, 'PaperPosition', [0 0 6 4.5],'Papersize',[6 4.5]);
hold all
plot(l,grapsd{4}(1:end,1),'b--','linewidth',2)
plot(l,grapsd{2}(1:end,1),'b-','linewidth',2)
plot(l,grapsd{4}(1:end,4),'g--','linewidth',2)
plot(l,grapsd{2}(1:end,4),'g-','linewidth',2)
plot(l,grapsd{4}(1:end,5),'r--','linewidth',2)
plot(l,grapsd{2}(1:end,5),'r-','linewidth',2)
% legend('Predicted g_N','GRACE observed g_N','Predicted g_E','GRACE observed g_E ','Predicted g_D','GRACE observed g_D','Location','NorthWest')% legend('location','NorthEast')
legend('Predicted Txx','GRACE Txx','Predicted Txy','GRACE Txy ','Predicted Txz','GRACE Txz','Location','NorthWest')% legend('location','NorthEast')
xlabel('Spherical Harmonic Degree');
ylabel('Degree Variance mE^2');
legend('location','NorthWest')
box on
grid on
% axis([0 100 0 0.2])
% ofile='PSD_JPcopo_Gra_bw900_R0';
ofile=['PSD_',datatype,'grad','_',texteq];
saveas(gcf,ofile,'fig');
print('-dpdf','-r500',ofile)


% idt=3;
% figure
% set(gcf,'Color','white')
% set(gca,'FontSize', 12);
% set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
% hold all
% plot(l,grapsd{idt}(1:end,1),'b--','linewidth',2)
% plot(l,grapsd{idt}(1:end,2),'g--','linewidth',2)
% plot(l,grapsd{idt}(1:end,3),'r--','linewidth',2)
% legend('Predicted g_N','Predicted g_E','Predicted g_D','location','NorthWest')
% legend('Predicted g_N','Predicted g_E','Predicted g_D','location','NorthEast')
% legend('location','NorthEast')
% xlabel('Spherical Harmonic Degree');
% ylabel('Degree Variance \muGal^2');
% box on
% grid on
% % ofile='PSD_JPcopo_Gra_bw900_R0';
% ofile=['PSD_',datatypeg{1,idt}];
% saveas(gcf,ofile,'fig');
% print('-dpdf','-r500',ofile)
% 
% figure
% set(gcf,'Color','white')
% set(gca,'FontSize', 12);
% hold all
% plot(l,grapsdperc{idt}(1:end,1),'b--','linewidth',2)
% plot(l,grapsdperc{idt}(1:end,2),'g--','linewidth',2)
% plot(l,grapsdperc{idt}(1:end,3),'r--','linewidth',2)
% legend('Predicted g_N','Predicted g_E','Predicted g_D','location','NorthWest')
% legend('Predicted g_N','Predicted g_E','Predicted g_D','location','SouthEast')
% xlabel('Spherical Harmonic Degree');
% ylabel('Degree Variance accumulated percentage');
% plot([1, nl],[90, 90],'y-')
% box on
% grid on
% % ofile='PSD_perc_JPcopo_Gra_bw900_R0';
% ofile=['PSD_perc_',datatypeg{1,idt}];
% saveas(gcf,ofile,'fig');
% print('-dpdf','-r500',ofile)

% idt=4;
% figure
% set(gcf,'Color','white')
% set(gca,'FontSize', 12);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
% hold all
% semilogy(l,grapsd{idt}(1:end,1),'b-','linewidth',2)
% semilogy(l,grapsd{idt}(1:end,4),'g-','linewidth',2)
% semilogy(l,grapsd{idt}(1:end,5),'r-','linewidth',2)
% legend('Predicted Txx','Predicted Txy','Predicted Txz','location','NorthEast')
% xlabel('Spherical Harmonic Degree');
% ylabel('Degree Variance mE^2');
% box on
% grid off
% % ofile='PSD_JPcopo_Gra_bw900_R0';
% ofile=[outdir,'PSD_',datatypeg{1,idt}];
% saveas(gcf,ofile,'fig');
% print('-dpdf','-r500',ofile)

% 
% figure
% set(gcf,'Color','white')
% set(gca,'FontSize', 12);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
% hold all
% plot(l,grapsdperc{idt}(1:end,1),'b-','linewidth',2)
% plot(l,grapsdperc{idt}(1:end,4),'g-','linewidth',2)
% plot(l,grapsdperc{idt}(1:end,5),'r-','linewidth',2)
% legend('Predicted Txx','Predicted Txy','Predicted Txz','location','SouthEast')
% xlabel('Spherical Harmonic Degree');
% ylabel('Degree Variance accumulated percentage');
% plot([1, nl],[90, 90],'y-')
% box on
% grid on
% % ofile='PSD_perc_JPcopo_Gra_bw900_R0';
% ofile=[outdir,'PSD_perc_',datatypeg{1,idt}];
% saveas(gcf,ofile,'fig');
% print('-dpdf','-r500',ofile)

% Refers to lunar_grav_topo/Slep_expansion_topo_vs_tapering.m
% [test,dw]=xyz2plm(Comtopo_w,180,'im',[],[]);
% 
    

