%================ plots for the first draft ===============================

%--- fig1 model configurations 
load manuscriptvort_ideal.mat
load figureplot2.mat
x1=x*1e-3;
y1=y*1e-3;
[X,Y]=meshgrid(x1(2:end-1),y1(2:end-1));

close all
x0=10;
y0=10;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
% ax(1)=subplot(1,2,1);
% ax(1).Position=[0.1 0.12 0.2 0.8];% [left bottom width height]
ax(1)=subplot(1,3,1);
ax(1).Position=[0.1 0.12 0.18 0.8];% [left bottom width height]
rho0=999.8;
plot(rho0*windstress,y1,'LineWidth',2)
hold on
plot(rho0*windstressdw,y1,'LineWidth',2)
grid on
xlabel('Wind stress [N m^{-2}]')
ylabel('Meridional distance [km]')
legend({'reference-wind','double-wind'},'fontsize',10,'Location','best')
title('(a)','position',[0.012 y1(end-1)+20])
set(gca,'fontsize',14,'ytick',200:200:1800)


ax(2)=subplot(1,3,2);
ax(2).Position=[0.31 0.12 0.15 0.8];
plot(sflx(96,:),y1,'LineWidth',2)
grid on
xlabel('Heat flux [W m^{-2}]')
title('(b)','position',[-9 y1(end-1)+20])
set(gca,'fontsize',14,'ytick',[])

ax(3)=subplot(1,3,3);
ax(3).Position=[0.5 0.12 0.42 0.8];
imagesc(x1(2:end-1),y1(2:end-1),d(2:end-1,2:end-1)'),axis xy
colorbar
caxis([2000 3000])
colorbar('Ticks',2000:200:3000)
ylabel(colorbar,'Bathymetry [m]','fontsize',12)
% colormap(ax(1),'hot')
set(ax(3),'xtick',[],'ytick',[]);
h2=axes('position',get(ax(3),'position'),'color','none','fontsize',14); 
hold on
contour(X,Y,etabar(2:end-1,2:end-1)','LineWidth',2,'ShowText','on','parent',h2)
colormap(h2,[0,0,0])
% ylabel('Meridional distance [km]','fontsize',14)
title('(c)','position',[x1(5) y1(end-1)+1])
set(gca,'fontsize',14,'xtick',200:400:1800,'ytick',[])
xlabel('Zonal distance [km]')
print -dpng modeldomain.png

%--- fig2 sst snapshot 
load figureplot2.mat
x=x*1e-3;
y=y*1e-3;
close all
imagesc(x,y(2:end),sst(:,2:end)'),axis xy
colorbar
xlabel('Zonal distance [km]')
ylabel('Meridional distance [km]')
ylabel(colorbar,'[^{o}C]','fontsize',14)
set(gca,'fontsize',14)
caxis([2 8])
print -dpng sstsnap.png



%------ Fig3 surface Ekman layer establishment of stage0-------------------
% surface Ekman transport across mid-latitude where wind stress=0.1 N/m^2
load medtrans.mat
t=1:120;
upd=8;
surfref=sum(medtransref(:,1:upd),2)/1e6;
surfhom=sum(medtransuni(:,1:upd),2)/1e6;
close all
x0=10;
y0=10;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])

plot(t,surfref,'LineWidth',2)
hold on
plot(t,surfhom,'LineWidth',2)
grid on
plot(t,-(mean(tauLx)/1e6)*ones(size(t)),'k--','LineWidth',2)
xlabel('Time [day]','fontsize',14)
legend({'stratified','homogenous','theoretical'},'fontsize',14)
set(gca,'fontsize',14)
ylabel('[Sv]')
print -dpng stage0.png

%-------- Fig4 SSH hovomollers -----------------
%-- zonally averaged SSH
load figureplot2.mat
y1=y*1e-3;
t=1:120;
sth=39;
nth=160;
close all
x0=10;
y0=10;
width=800;
height=700;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(3,1,1);
imagesc(t,y1(2:end),etameanhom(t,2:end)'),axis xy
colorbar
caxis([-0.04 0.04])
hold on
plot(t,y1(sth)*ones(size(t)),'k--')
plot(t,y1(nth)*ones(size(t)),'k--')
set(ax(1),'ytick',0:500:2000,'fontsize',14)
ylabel('Meridional distance [km]')
ylabel(colorbar,'[m]')
title('(a)','position',[2 2001])

ax(2)=subplot(3,1,2);
imagesc(t,y1(2:end),etameanref(t,2:end)'),axis xy
colorbar
caxis([-0.04 0.04])
hold on
plot(t,y1(sth)*ones(size(t)),'k--')
plot(t,y1(nth)*ones(size(t)),'k--')
set(ax(2),'ytick',0:500:2000,'fontsize',14)
ylabel('Meridional distance [km]')
ylabel(colorbar,'[m]')
title('(b)','position',[2 2001])

ax(3)=subplot(3,1,3);
plot(t,etameanref(:,nth)-etameanref(:,sth),'LineWidth',2)
hold on
plot(t,etameanhom(:,nth)-etameanhom(:,sth),'LineWidth',2)
grid on
legend({'stratified','homogenous'},'fontsize',14)
ylabel('[m]')
xlabel('Time [day]')
set(ax(3),'fontsize',14)
title('(c)','position',[2 0.101])
print -dpng sshstage1.png


%---- fig 5 zonal velocity profiles (zonally averaged at mid latitude)
load zonalu.mat
t=1:120;
[T,Z]=meshgrid(t,zc);
close all
x0=10;
y0=10;
width=800;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
homprof=nan(size(uprofmidlat_hom));
homprof(1:end-1,1:end-1)=uprofmidlat_hom(2:end,2:end);
pcolor(T,Z,homprof')
shading flat
colorbar
caxis([0 1e-2])
set(ax(1),'fontsize',14,'layer','top')
ylim([zc(end) zc(1)])
ylabel('Depth [m]')
ylabel(colorbar,'[m s^{-1}]','fontsize',14)
title('(a)','position',[2 0.01])

ax(2)=subplot(2,1,2);
refprof=nan(size(uprofmidlat_ref));
refprof(1:end-1,1:end-1)=uprofmidlat_ref(2:end,2:end);
pcolor(T,Z,refprof')
shading flat
colorbar
caxis([0 1e-2])
set(ax(1),'fontsize',14,'layer','top')
ylim([zc(end) zc(1)])
ylabel('Depth [m]')
ylabel(colorbar,'[m s^{-1}]','fontsize',14)
title('(b)','position',[2 0.01])
set(ax(2),'fontsize',14,'layer','top')
xlabel('Time [day]')
print -dpng fig3_uprofile




%------ Fig6 adjustment of bottom Ekman layer---------------
load medtrans.mat
rgb = get(gca,'colororder');
wid=1;
t=1:wid:120;
upd=8;
lowd=36;
surfref=sum(medtransref(:,1:upd),2)/1e6; % meridionally averaged
surfhom=sum(medtransuni(:,1:upd),2)/1e6; % meridionally averaged
botref=sum(medtransref(:,lowd:end),2)/1e6;% meridionally averaged
bothom=sum(medtransuni(:,lowd:end),2)/1e6;% meridionally averaged


refmeansurf=mean(reshape(surfref,wid,[]),1);
refmeanbot=mean(reshape(botref,wid,[]),1);%*
homeansurf=mean(reshape(surfhom,wid,[]),1);
homeanbot=mean(reshape(bothom,wid,[]),1);%*


reftot=surfref+botref;
homtot=surfhom+bothom;
refmean=mean(reshape(reftot,wid,[]),1);% cross mid-latitude
homean=mean(reshape(homtot,wid,[]),1);% cross mid-latitude


close all
x0=10;
y0=10;
width=800;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
plot(t,-refmeanbot,'LineWidth',2)
hold on
plot(t,-homeanbot,'LineWidth',2)
plot(t,surfhom(end)*ones(size(t)),'k--','LineWidth',2)

grid on
legend({'stratified','homogenous','Ekman transport'},'fontsize',14)

set(ax(1),'fontsize',14)
ylabel('[Sv]')
ylim([0.4 2.2])
title('(a)','position',[2 2.21])


ax(2)=subplot(2,1,2);
plot(t,refmean,'LineWidth',2)
hold on
plot(t,homean,'LineWidth',2)
grid on

set(ax(2),'fontsize',14)
ylabel('[Sv]')
xlabel('Time [day]')
legend({'stratified','homogenous'},'fontsize',14,'location','best')
ylim([-0.6 1])
title('(b)','position',[2 1.01])

print -dpng fig4_volume.png



%------ Fig.7 time series of TFS ----------------------
% area-averaged stresses 
raw=10420*10420*192*192;
rho0=999.8;
load mmtdaily.mat
load figureplot2.mat
close all
t=1:120;
x0=10;
y0=10;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
plot(t,-rho0*tfsref/raw,'LineWidth',2)
hold on
plot(t,-rho0*tfshom/raw,'LineWidth',2)
plot(t,rho0*ext/raw,'k--','LineWidth',2)
grid on 
xlabel('Time [day]')
ylabel('[N m^{-2}]')
% ylabel('[m^{4}s^{-2}]')
legend({'-TFS stratified','-TFS homogenous','wind stress'},'fontsize',14)
set(gca,'fontsize',14,'Layer','top')
print -dpng TFScomp.png

%-- Fig8 isotherm slopes (stage0-2a)
%- computed from rhoanmination.m
load Tcontours.mat
yled=1000:100:1300;
zled1=-2200;
zled2=-2400;
zled3=-2600;
y=x;
y1=y*1e-3;
[Y,Z]=meshgrid(y1(2:end),zc(9:end));

close all
figure
ax(1)=gca;
tc=0:1:6;
contour(Y,Z,T1day(2:end,9:end)',tc,'--','ShowText','off','LineWidth',2)
colormap(ax(1),[0,0,0])
set(ax(1),'XTick',[],'YTick',[])
h1=axes('position',get(ax(1),'position'),'color','none','fontsize',14);
hold on
contour(Y,Z,T2month(2:end,9:end)',tc,'ShowText','on','LineWidth',2)
colormap(h1,[1,0,0])
h2=axes('position',get(h1,'position'),'color','none','fontsize',14);
hold on
contour(Y,Z,T6year(2:end,9:end)',tc,'ShowText','off','LineWidth',2)
colormap(h2,[0,0,0])
hold on
plot(yled,zled1*ones(size(yled)),'k--','LineWidth',2)
plot(yled,zled2*ones(size(yled)),'color',[1,0,0],'LineWidth',2)
plot(yled,zled3*ones(size(yled)),'color',[0,0,0],'LineWidth',2)
text(1400,zled1,'Initial','fontsize',12)
text(1400,zled2,'Month 2','fontsize',12)
text(1400,zled3,'Year 6','fontsize',12)

xlabel('Meridional distance [km]')
ylabel('Depth [m]')
print -dpng fig6_isotherms.png




%-- Fig 9 transport separation druing first 4 months and long time ---------
load timescale.mat
load figureplot2.mat
close all
tshort=1:120;
tlong=1:100;
x0=10;
y0=10;
width=800;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
plot(tshort,transref,'LineWidth',2)
hold on
plot(tshort,BCref,'LineWidth',2)
plot(tshort,BTref,'LineWidth',2)
grid on
% legend({'total','baroclinic','barotropic'},'fontsize',14,'location','best')
set(gca,'fontsize',14)
xlabel('Time [day]')
ylabel('[Sv]')
title('(a)','position',[2 15.1])
% print -dpng transportspinup.png

ax(2)=subplot(2,1,2);

cbar=get(gca,'colororder');
plot(1:100,trans100,'color',cbar(1,:),'LineWidth',2)
hold on
plot(1:100,BC100,'color',cbar(2,:),'LineWidth',2)
plot(1:100,BT100,'color',cbar(3,:),'LineWidth',2)
grid on
plot(101:130,transdw,'-.','color',cbar(1,:),'LineWidth',2)
plot(101:130,BCdw,'-.','color',cbar(2,:),'LineWidth',2)
plot(101:130,BTdw,'-.','color',cbar(3,:),'LineWidth',2)
legend({'total','baroclinic','barotropic','total double-wind','baroclinic double-wind','barotropic double-wind','Year 6'},'fontsize',10,'location','best')
plot(6*ones(1,100),1:100,'k--','LineWidth',2)
% set(ax(1),'fontsize',14)
xlabel('Time [year]')
ylabel('[Sv]')
title('(b)','position',[2 100.1])
set(gca,'fontsize',14)
print -dpng transportspinup.png



%-- Fig 10 time series of heat transport after doubling wind (separate time-mean and time varying field)
load figureplot2.mat
cbar=get(gca,'colororder');
% y1=y*1e-3;
close all
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
la=96;

% ax(1)=subplot(2,1,1);
plot(1:100,HTEts,'LineWidth',2)
hold on
plot(1:100,Hekmants+HSEts,'LineWidth',2)
plot(1:100,Htotalts,'LineWidth',2)
plot(1:130,-Hsuf(la)*ones(size(1:130)),'LineWidth',2)
plot(101:130,HTEtsdw,'-.','color',cbar(1,:),'LineWidth',2)
plot(101:130,Hekmantsdw+HSEtsdw,'-.','color',cbar(2,:),'LineWidth',2)
plot(101:130,Htotaltsdw,'-.','color',cbar(3,:),'LineWidth',2)
plot(6*ones(size(-30:60)),-30:60,'k--','LineWidth',2)
% ylim([-50 100])
grid on
legend({'eddy','mean','total','surface',...
    'eddy double-wind','mean double-wind','total double-wind','Year 6'},'fontsize',12,'location','best')
% legend({'H_{eddy}','H_{mean}','H_{total}','H_{surface}',...
%     'H_{eddy} double wind','H_{mean} double wind','H_{total} double wind'},'fontsize',12,'location','bestoutside')
set(gca,'fontsize',14)
ylabel('[TW]')
xlabel('Time [year]','fontsize',14)
% print -dpng heatdw30yr.png
% print -dpng heatransportdweddy.png
print -dpng fig8_heatflux.png

%===== fig 11 long term meridional transport and SSH hovemolloer============
% remove long time meridional transport
% load medtrans.mat
load figureplot2.mat
t=1:100;
y1=y*1e-3;
close all
x0=10;
y0=10;
width=900;
height=400;
set(gcf,'position',[x0,y0,width,height])

imagesc(t,y1(2:end),etameanref100(t,2:end)'),axis xy
colorbar
caxis([-0.5 0.5])
xlabel('Time [year]')
ylabel('Meridional distance [km]')
set(gca,'fontsize',14)
ylabel(colorbar,'[m]')

print -dpng volumeandssh.png

%-- alternatively fig11 (add one more panel to have the line plot)
%----------------------
load figureplot2.mat
sth=39;
nth=160;
t=1:100;
y1=y*1e-3;
close all
x0=10;
y0=10;
width=900;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
imagesc(t,y1(2:end),etameanref100(t,2:end)'),axis xy
colorbar
caxis([-0.5 0.5])
hold on
plot(t,y1(sth)*ones(size(t)),'k--')
plot(t,y1(nth)*ones(size(t)),'k--')
xlabel('Time [year]')
ylabel('Meridional distance [km]')
set(gca,'fontsize',14)
ylabel(colorbar,'[m]')
title('(a)','position',[2 2001])


ax(2)=subplot(2,1,2);
plot(t,etameanref100(:,nth)-etameanref100(:,sth),'LineWidth',2)
grid on
ylabel('[m]')
xlabel('Time [year]')
set(ax(2),'fontsize',14)
title('(b)','position',[2 0.801])
print -dpng fig9_sshslope.png


%-- Fig12 isotherm slopes (doubling wind)
load Tcontours.mat
yled=1000:100:1300;
zled1=-2200;
zled2=-2400;
zled3=-2600;
y=x;
y1=y*1e-3;
[Y,Z]=meshgrid(y1(2:end),zc(9:end));
tc=0:1:6;
% h3=subplot(1,2,2);
figure
h3=gca;
contour(Y,Z,Tfinal(2:end,9:end)',tc,'ShowText','on','LineWidth',2)
set(h3,'XTick',[],'Ytick',[])
colormap(h3,[0,0,1])
h4=axes('position',get(h3,'position'),'color','none','fontsize',14);
hold on
contour(Y,Z,Tdwind(2:end,9:end)',tc,'ShowText','off','LineWidth',2)
colormap(h4,[0.5,0.5,0.5])
xlabel('Meridional distance [km]')
ylabel('Depth [m]')
% title('(b)','position',[y1(5) zc(9)+0.1])
hold on
% plot(yled,zled1*ones(size(yled)),'color',[0,0,0],'LineWidth',2)
plot(yled,zled2*ones(size(yled)),'color',[0,0,1],'LineWidth',2)
plot(yled,zled3*ones(size(yled)),'color',[0.5,0.5,0.5],'LineWidth',2)
% text(1400,zled1,'Year 6','fontsize',12)
text(1400,zled2,'reference-wind','fontsize',12)
text(1400,zled3,'double-wind','fontsize',12)
print -dpng fig11_isotherms.png

%Fig 13 hovomller of SSH and meridional SSH difference after doubling wind
load figureplot2.mat
y1=y*1e-3;
% original is sth=39; north=154
% sth=30;
% nth=154;
sth=39;
nth=160;
% nth=165;
close all
x0=10;
y0=10;
width=800;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
t=1:120;
imagesc(t,y1(2:end),etameanrefdwhf(t,2:end)'),axis xy
colorbar
caxis([-0.5 0.5])
hold on
% plot(t,y1(5)*ones(size(t)),'k--','LineWidth',2)
% plot(t,y1(end-5)*ones(size(t)),'k--','LineWidth',2)
plot(t,y1(sth)*ones(size(t)),'k--','LineWidth',2)
plot(t,y1(nth)*ones(size(t)),'k--','LineWidth',2)
set(ax(1),'fontsize',14)
xlabel('Time [day]')
ylabel('Meridional distance [km]')
ylabel(colorbar,'[m]')
title('(a)','position',[2 2000.1])
% print -dpng sshmeanovdbhf.png

ax(2)=subplot(2,1,2);
% plot(t,etameanrefdwhf(:,end-5)-etameanrefdwhf(:,5),'LineWidth',2)
plot(t,etameanrefdwhf(:,nth)-etameanrefdwhf(:,sth),'LineWidth',2)
grid on
ylabel('[m]')
xlabel('Time [day]')
% ylim([0.8 1])
% ylim([0.5 0.7])
ylim([0.65 0.75])
% title('(b)','position',[2 1.001])
title('(b)','position',[2 0.751])
set(ax(2),'fontsize',14)
print -dpng sshdiffhov.png
% print -dpng sshmeanovdbhfdif.png

% %--Fig 14 time seires of TFS in both simulations after doubling wind
%-+meridional transport for stratified & homogenous after doubling wind
load medtrans.mat
upd=8;
lowd=36;
surfrefdw=sum(medtransdbhf(:,1:upd),2)/1e6; % meridionally averaged 
botrefdw=sum(medtransdbhf(:,lowd:end),2)/1e6;% meridionally averaged 
surfhomdw=sum(medtranshomdwhf(:,1:upd),2)/1e6;% meridionally averaged 
bothomdw=sum(medtranshomdwhf(:,lowd:end),2)/1e6;% meridionally averaged 

close all
x0=10;
y0=10;
width=800;
% height=700;
height=500;
set(gcf,'position',[x0,y0,width,height])

t=1:120;
ax(1)=subplot(2,1,1);
plot(t,-botrefdw,'LineWidth',2)
hold on
plot(t,-bothomdw,'LineWidth',2)
plot(t,surfhomdw(end)*ones(size(t)),'k--','LineWidth',2)
grid on
ylim([1.5 3.6])
legend({'stratified','homogenous','Ekman transport'},'fontsize',14,'location','best')
xlabel('Time [day]')
ylabel('[Sv]')
set(ax(1),'fontsize',14)
title('(a)','position',[2 3.61])

% print -dpng medtransdw.png
raw=10420*10420*192*192;
rho0=999.8;
load mmtdaily.mat
load figureplot2.mat
ax(2)=subplot(2,1,2);
plot(1:120,-rho0*tfsrefdbhf/raw,'LineWidth',2)
hold on
plot(1:120,-rho0*tfshomdbhf/raw,'LineWidth',2)
plot(1:120,rho0*extdb/raw,'k--','LineWidth',2)
legend({'-TFS stratified','-TFS homogenous','wind stress'},'fontsize',14,'Location','best')
ylim([0.05 0.18])
grid on
xlabel('Time [day]')
ylabel('[N m^{-2}]')
title('(b)','Position',[2,0.181])
set(gca,'fontsize',14)
print -dpng windtfsdw.png

%======= barotropic structure adjustment ==================================
% %--Fig15 cross-stream heat fluxes (divergent component)
load crosflxbc.mat
load heatmaps.mat
load figureplot2.mat
phi2dnewextref=zeros(192,188);
phi2dnewextref(2:end,:)=phinewref;
phi2dnewextref(1,:)=phinewref(end,:);

Fdivxref=(phi2dnewextref(2:end,:)-phi2dnewextref(1:end-1,:))/10420;%191*188 u grid 
Fdivxref1=Fdivxref(:,2:end);% 191*187
Fdivyref=(phinewref(:,2:end)-phinewref(:,1:end-1))/10420;% 191*187 v grid 
Fdivyref1=Fdivyref;

% Tdznewref=Tdzref(:,3:end-2);% 192*188
Tdznewref=etabaref(:,3:end-2);% 192*188
DXref=(Tdznewref(2:end,:)-Tdznewref(1:end-1,:))/10420;% 191*188 u grid 
DYref=(Tdznewref(:,2:end)-Tdznewref(:,1:end-1))/10420;% 192*187 v grid 
DXref1=DXref(:,2:end);% 191*187
DYref1=DYref(2:end,:);% 191*187
nmref=sqrt(DXref1.^2+DYref1.^2);% 191*187
dotabref=Fdivxref1.*DXref1+Fdivyref1.*DYref1;
crosshflxref=dotabref./nmref;

phi2dnewextdw=zeros(192,188);
phi2dnewextdw(2:end,:)=phinewdw;
phi2dnewextdw(1,:)=phinewdw(end,:);

Fdivxdw=(phi2dnewextdw(2:end,:)-phi2dnewextdw(1:end-1,:))/10420;%191*188 u grid 
Fdivxdw1=Fdivxdw(:,2:end);% 191*187
Fdivydw=(phinewdw(:,2:end)-phinewdw(:,1:end-1))/10420;% 191*187 v grid 
Fdivydw1=Fdivydw;

% Tdznewdw=Tdzdw(:,3:end-2);% 192*188
Tdznewdw=etabardw(:,3:end-2);% 192*188
DXdw=(Tdznewdw(2:end,:)-Tdznewdw(1:end-1,:))/10420;% 191*188 u grid 
DYdw=(Tdznewdw(:,2:end)-Tdznewdw(:,1:end-1))/10420;% 192*187 v grid 
DXdw1=DXdw(:,2:end);% 191*187
DYdw1=DYref(2:end,:);% 191*187
nmdw=sqrt(DXdw1.^2+DYdw1.^2);% 191*187
% dotab=Fdivxref.*DX+Fdivyref.*DY;
dotabdw=Fdivxdw1.*DXdw1+Fdivydw1.*DYdw1;
crosshflxdw=dotabdw./nmdw;


y1=y*1e-3;
x1=x*1e-3;
% tc=1:0.2:2.4;
tc=-0.3:0.1:0.5;
[X,Y]=meshgrid(x1(1:end-1),y1(4:end-2));


close all
figure
x0=10;
y0=10;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(1,2,1);
ax(1).Position=[0.1 0.12 0.38 0.8];% [left bottom width height]
imagesc(x1(2:end),y1(4:end-2),crosshflxref'),axis xy
% colorbar
% caxis([-20 20])
caxis([-10 10])
set(ax(1),'XTick',[],'YTick',[])
% caxis([-25 25])
% ylabel(colorbar,'[^{o}C m^{2} s^{-1}]','fontsize',14)
h2=axes('position',get(ax(1),'position'),'color','none','fontsize',10);
hold on
[ca,ha]=contour(X,Y,etabaref(1:end-1,4:end-2)',tc,'showtext','off','LineWidth',0.5);
ha.LevelList=round(ha.LevelList,1);
colormap(h2,[1,1,1])
clabel(ca,ha,'manual','fontsize',12)
% title('cross stream F^{div}_{eddy} reference wind')
xlabel('Zonal distance [km]')
ylabel('Meridional distance [km]')
set(gca,'fontsize',14)
title('(a)','position',[x1(5) y1(end-1)+5])

ax(2)=subplot(1,2,2);
ax(2).Position=[0.5 0.12 0.42 0.8];% [left bottom width height]
imagesc(x1(2:end),y1(4:end-2),crosshflxdw'),axis xy
colorbar
% caxis([-20 20])
caxis([-10 10])
set(ax(2),'XTick',[],'YTick',[])
% caxis([-25 25])
ylabel(colorbar,'[^{o}C m^{2} s^{-1}]','fontsize',14)
h3=axes('position',get(ax(2),'position'),'color','none','fontsize',10);
hold on
[cb,hb]=contour(X,Y,etabardw(1:end-1,4:end-2)',tc,'showtext','off','LineWidth',0.5);
hb.LevelList=round(hb.LevelList,1);
colormap(h3,[1,1,1])
clabel(cb,hb,'manual','fontsize',12)
set(h3,'YTick',[])
% title('cross stream F^{div}_{eddy} reference wind')
xlabel('Zonal distance [km]')
% ylabel('Y [km]')
set(gca,'fontsize',14)
title('(b)','position',[x1(5) y1(end-1)+5])
print -dpng Fdivtheta.png





%-- Fig16 time series of SSH contour length ---------
load sshsenseideal.mat
close all
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
plot(1:100,mean(sshtsref(:,5:8),2,'omitnan')/1e3,'LineWidth',2)% average over 0.1 to 0.4
% plot(1:100,mean(sshtsref(:,1:end-1),2,'omitnan'),'LineWidth',2)% average over -0.3 to 0.4
hold on
plot(101:130,mean(sshtsdw(:,5:8),2,'omitnan')/1e3,'LineWidth',2)% average over 0.1 to 0.4
% plot(101:130,mean(sshtsdw(:,1:end-1),2,'omitnan'),'LineWidth',2)% average over -0.3 to 0.4
grid on
legend({'reference-wind','double-wind'},'fontsize',14,'Location','best')
xlabel('Time [year]')
ylabel('Contour length [km]')
set(gca,'fontsize',14)
print -dpng sshtsideal.png



%-fig15 long time series of SSH and meridional transport of stratified simulation
load figureplot2.mat
load medtrans.mat
close all
% x0=10;
% y0=10;
% width=800;
% height=500;
ax(1)=subplot(2,1,1);
plot(1:100,etameanref100(:,end-5)-etameanref100(:,5),'LineWidth',2)
hold on
plot(101:130,etameanrefdw(:,end-5)-etameanrefdw(:,5),'LineWidth',2)
grid on
% legend({'reference','double wind'},'fontsize',14)
% xlabel('Time [year]')
ylabel('[m]')
set(gca,'fontsize',14)
title('(a)','position',[2 1.01])
% print -dpng sshmeandiff.png

ax(2)=subplot(2,1,2);
plot(1:100,sum(medtransmid100,2)/1e6,'LineWidth',2)
hold on
plot(101:130,sum(medtransmidw,2)/1e6,'LineWidth',2)
grid on
ylim([-2e-3 5e-3])
legend({'stratified reference wind','stratified double wind'},'fontsize',14)
xlabel('Time [year]')
ylabel('[Sv]')
set(gca,'fontsize',14)
title('(b)','position',[20 5.1e-3])
print -dpng BTdiag.png
% print -dpng totmedmidtransdw.png

%-- fig 16 time series of heat transport after doubling wind
% computed from heatts.m
load figureplot2.mat
cbar=get(gca,'colororder');
% y1=y*1e-3;
close all
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
la=96;

% ax(1)=subplot(2,1,1);
plot(1:100,HTEts+HSEts,'LineWidth',2)
hold on
plot(1:100,Hekmants,'LineWidth',2)
plot(1:100,Htotalts,'LineWidth',2)
plot(1:100,-Hsuf(la)*ones(size(1:100)),'LineWidth',2)
plot(101:130,HTEtsdw+HSEtsdw,'-.','color',cbar(1,:),'LineWidth',2)
plot(101:130,Hekmantsdw,'-.','color',cbar(2,:),'LineWidth',2)
plot(101:130,Htotaltsdw,'-.','color',cbar(3,:),'LineWidth',2)
% ylim([-50 100])
grid on
legend({'H_{eddy}','H_{Ekman}','H_{total}','H_{surface}',...
    'H_{eddy} double wind','H_{Ekman} double wind','H_{total} double wind'},'fontsize',12,'location','bestoutside')
set(gca,'fontsize',14)
ylabel('[TW]')
xlabel('Time [year]','fontsize',14)
print -dpng ~/heatransportdw.png





