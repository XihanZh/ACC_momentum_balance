%================ plots for the first draft ===============================

%--- fig1 model configurations 
x1=x*1e-3;
y1=y*1e-3;
[X,Y]=meshgrid(x1(2:end-1),y1(2:end-1));

close all
x0=10;
y0=10;
width=800;
height=1000;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,3,1);
ax(1).Position=[0.1 0.52 0.18 0.4];% [left bottom width height]
rho0=999.8;
plot(rho0*windstress,y1,'LineWidth',2)
hold on
plot(rho0*windstressdw,y1,'LineWidth',2)
grid on
xlabel('Wind stress [N m^{-2}]')
ylabel('Meridional distance [km]')
legend({'reference-wind','double-wind'},'fontsize',12,'Location','best')
title('(a)','position',[0.012 y1(end-1)+20])
set(gca,'fontsize',14,'ytick',200:200:1800)


ax(2)=subplot(2,3,2);
ax(2).Position=[0.32 0.52 0.18 0.4];% [left bottom width height]
plot(sflx(96,:),y1,'LineWidth',2)
grid on
xlabel('Heat flux [W m^{-2}]')
title('(b)','position',[-9 y1(end-1)+20])
set(gca,'fontsize',14,'ytick',[])

ax(3)=subplot(2,3,3);
ax(3).Position=[0.55 0.52 0.42 0.4];% [left bottom width height]
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

ax(4)=subplot(2,3,4);
ax(4).Position=[0.1 0.08 0.4 0.35];% [left bottom width height]
imagesc(x1(2:end-1),y1(2:end-1),sst(2:end-1,2:end-1)'),axis xy
colorbar
xlabel('Zonal distance [km]')
set(gca,'fontsize',14,'xtick',200:400:1800,'ytick',200:200:1800)
ylabel('Meridional distance [km]')
ylabel(colorbar,'[^{o}C]','fontsize',14)
% set(gca,'fontsize',14)
caxis([2 8])
title('(d)','position',[x1(5) y1(end-1)+1])

ax(5)=subplot(2,3,5);
ax(5).Position=[0.6 0.08 0.38 0.35];% [left bottom width height]
plot(tend,y1,'LineWidth',2)
hold on
plot(adv,y1,'LineWidth',2)
plot(cori,y1,'LineWidth',2)
plot(diss+vis,y1,'LineWidth',2)
% plot(diss,y1,'LineWidth',2)
plot(dphi,y1,'LineWidth',2)
plot(ext,y1,'LineWidth',2)
% plot(vis,y1,'LineWidth',2)
plot(res,y1,'k','LineWidth',2)
% legend({'tendency','advection','Coriolis','bottom friction','PGF(TFS)','wind stress','vertical viscosity','residual'},'fontsize',10,'location','northeast')
legend({'tendency','advection','Coriolis','friction','TFS','wind stress','residual'},'fontsize',10,'location','northeast')
grid on
xlabel('momentum terms [m^{3}s^{-2}]','fontsize',14)
set(gca,'fontsize',14,'ytick',200:200:1800)
xlim([-200 200])
title('(e)','position',[-195 y1(end-1)+1])
print -dpng ~/Desktop/fig1_modeldomain.png

%---fig3 stage0 and stage 1
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

raw=10420*10420*192*192;
rho0=999.8;
load mmtdaily.mat
load figureplot2.mat

close all
x0=10;
y0=10;
width=800;
height=750;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(3,1,1);
plot(t,surfref,'LineWidth',2)
hold on
plot(t,surfhom,'LineWidth',2)
grid on
plot(t,-(mean(tauLx)/1e6)*ones(size(t)),'k--','LineWidth',2)
ylim([1.2 1.55])
xlabel('Time [day]','fontsize',14)
legend({'stratified','homogenous','theoretical'},'location','best','fontsize',14)
set(gca,'fontsize',14)
ylabel('[Sv]')
title('(a)','position',[2 1.552])


ax(2)=subplot(3,1,2);
plot(t,refmean,'LineWidth',2)
hold on
plot(t,homean,'LineWidth',2)
grid on

set(ax(2),'fontsize',14)
ylabel('[Sv]')
% xlabel('Time [day]')
legend({'stratified','homogenous'},'fontsize',14,'location','best')
ylim([-0.6 1])
title('(b)','position',[2 1.01])

ax(3)=subplot(3,1,3);
plot(t,-rho0*tfsref/raw,'LineWidth',2)
hold on
plot(t,-rho0*tfshom/raw,'LineWidth',2)
plot(t,rho0*ext/raw,'k--','LineWidth',2)
grid on 
xlabel('Time [day]')
ylabel('[N m^{-2}]')
% ylabel('[m^{4}s^{-2}]')
legend({'-TFS stratified','-TFS homogenous','wind stress'},'location','best','fontsize',14)
set(gca,'fontsize',14,'Layer','top')
title('(c)','position',[2 0.107])
print -dpng stage1.png

%---fig4 barotropic flow 

load manuscriptvort_ideal.mat
load figureplot2.mat
load zonalu.mat
y1=y*1e-3;
t=1:120;
sth=39;
nth=160;
[T,Z]=meshgrid(t,zc);

close all
x0=10;
y0=10;
width=800;
height=1000;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(3,2,1);
% ax(1).Position=[0.1 0.65 0.38 0.3];% [left bottom width height]
ax(1).Position=[0.1 0.65 0.36 0.3];% [left bottom width height]
imagesc(t,y1(2:end),etameanhom(t,2:end)'),axis xy
caxis([-0.04 0.04])
hold on
plot(t,y1(sth)*ones(size(t)),'k--')
plot(t,y1(nth)*ones(size(t)),'k--')
set(ax(1),'ytick',0:500:2000,'fontsize',12)
ylabel('Meridional distance [km]')
title('(a)','position',[2 2001])


ax(2)=subplot(3,2,2);
% ax(2).Position=[0.55 0.65 0.4 0.3];% [left bottom width height]
ax(2).Position=[0.52 0.65 0.42 0.3];% [left bottom width height]
imagesc(t,y1(2:end),etameanref(t,2:end)'),axis xy
colorbar
caxis([-0.04 0.04])
hold on
plot(t,y1(sth)*ones(size(t)),'k--')
plot(t,y1(nth)*ones(size(t)),'k--')
set(ax(2),'ytick',0:500:2000,'fontsize',12)
ylabel(colorbar,'[m]','fontsize',14)
title('(b)','position',[2 2001])


ax(3)=subplot(3,2,3);
ax(3).Position=[0.1 0.32 0.36 0.3];% [left bottom width height]
homprof=nan(size(uprofmidlat_hom));
homprof(1:end-1,1:end-1)=uprofmidlat_hom(2:end,2:end);
pcolor(T,Z,homprof')
shading flat
% colorbar
caxis([0 1e-2])
set(ax(3),'fontsize',12,'layer','top')
ylim([zc(end) zc(1)])
ylabel('Depth [m]')
xlabel('Time [day]')
% ylabel(colorbar,'[m s^{-1}]','fontsize',14)
title('(c)','position',[2 0.01])


ax(4)=subplot(3,2,4);
ax(4).Position=[0.52 0.32 0.42 0.3];% [left bottom width height]
refprof=nan(size(uprofmidlat_ref));
refprof(1:end-1,1:end-1)=uprofmidlat_ref(2:end,2:end);
pcolor(T,Z,refprof')
shading flat
colorbar
caxis([0 1e-2])
set(ax(4),'fontsize',12,'layer','top')
ylim([zc(end) zc(1)])
% ylabel('Depth [m]')
ylabel(colorbar,'[m s^{-1}]','fontsize',14)
title('(d)','position',[2 0.01])
set(ax(4),'fontsize',12,'layer','top')
xlabel('Time [day]')


ax(5)=subplot(3,2,5);
ax(5).Position=[0.1 0.05 0.85 0.2];% [left bottom width height]
plot(t,etameanref(:,nth)-etameanref(:,sth),'LineWidth',2)
hold on
plot(t,etameanhom(:,nth)-etameanhom(:,sth),'LineWidth',2)
grid on
legend({'stratified','homogenous'},'fontsize',14)
ylabel('[m]')
xlabel('Time [day]')
set(ax(5),'fontsize',12)
title('(e)','position',[2 0.101])
print -dpng BTflow.png

%---fig5 transport separation druing first 4 months and long time
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
plot(6*ones(1,100),1:100,'k--','LineWidth',2)
legend({'total','baroclinic','barotropic','total double-wind','baroclinic double-wind','barotropic double-wind','Year 6'},'fontsize',10,'location','best')
% set(ax(1),'fontsize',14)
xlabel('Time [year]')
ylabel('[Sv]')
title('(b)','position',[2 100.1])
set(gca,'fontsize',14)
print -dpng transport.png


%---fig6 isotherm slopes (stage0-2a)
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


%---fig7 EP flux sections
load review.mat
load blue_red_saturated.mat
x1=x*1e-3;
% [X,Z]=meshgrid(x1,zc(2:end-1));
[X,Z]=meshgrid(x1,zc(2:end));
close all
x0=10;
y0=10;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,2,1);
epavg2anew=nan(size(epavg2a,1),size(epavg2a,2)+1);
epavg2anew(1:end-1,1:end-2)=epavg2a(2:end,2:end);
pcolor(X,Z,epavg2anew')
shading flat
hold on
plot(x1,-d(:,96),'-k','LineWidth',2)
ylim([-d(end,96),zc(1)])
colorbar
caxis([-5e-4 5e-4])
colormap(map)
set(ax(1),'fontsize',14,'Layer','top')
xlabel('Zonal distance [km]')
ylabel('Depth [m]')
ylabel(colorbar,'[m^{2}s^{-2}]','fontsize',14)
title('(a)','position',[50 0.1])

ax(2)=subplot(2,2,2);
% epavg2bnew=nan(size(epavg2b));
epavg2bnew=nan(size(epavg2b,1),size(epavg2b,2)+1);
% epavg2bnew(1:end-1,1:end-1)=epavg2b(2:end,2:end);
epavg2bnew(1:end-1,1:end-2)=epavg2b(2:end,2:end);
pcolor(X,Z,epavg2bnew')
shading flat
hold on
plot(x1,-d(:,96),'-k','LineWidth',2)
ylim([-d(end,96),zc(1)])
colorbar
caxis([-5e-4 5e-4])
colormap(map)
ylabel(colorbar,'[m^{2}s^{-2}]','fontsize',14)
set(ax(2),'fontsize',14,'Layer','top')
xlabel('Zonal distance [km]')
ylabel('Depth [m]')
title('(b)','position',[50 0.1])


ax(3)=subplot(2,2,3);
% epavgdwnew=nan(size(epavgdw));
epavgdwnew=nan(size(epavgdw,1),size(epavgdw,2)+1);
epavgdwnew(1:end-1,1:end-2)=epavgdw(2:end,2:end);
pcolor(X,Z,epavgdwnew')
shading flat
hold on
plot(x1,-d(:,96),'-k','LineWidth',2)
ylim([-d(end,96),zc(1)])
colorbar
caxis([-5e-4 5e-4])
colormap(map)
set(ax(3),'fontsize',14,'Layer','top')
xlabel('Zonal distance [km]')
ylabel('Depth [m]')
title('(c)','position',[50 0.1])
ylabel(colorbar,'[m^{2}s^{-2}]','fontsize',14)



ax(4)=subplot(2,2,4);
plot(eprofs2a(15:end),zc(16:end-1),'LineWidth',2)
hold on
plot(eprofs2b(15:end),zc(16:end-1),'LineWidth',2)
plot(eprofdw(15:end),zc(16:end-1),'LineWidth',2)
grid on
ylabel('Depth [m]')
% title('EP flux [m^{2}s^{-2}]')
xlabel('[m^{2}s^{-2}]','fontsize',14)
title('(d)','position',[2e-6 0.1])
set(gca,'fontsize',14)
legend({'Stage 2a','Stage 2b','double wind'},'fontsize',12,'location','best')
print -dpng EPfluxsections.png





%---fig 8 time series of heat transport after doubling wind (separate time-mean and time varying field)
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

%---fig9 (add one more panel to have the line plot)

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


%---fig10 isotherm slopes (doubling wind)
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
print -dpng fig9_isotherms.png

%---fig 11 hovomller of SSH and meridional SSH difference after doubling wind
load figureplot2.mat
y1=y*1e-3;
sth=39;
nth=160;
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

plot(t,y1(sth)*ones(size(t)),'k--','LineWidth',2)
plot(t,y1(nth)*ones(size(t)),'k--','LineWidth',2)
set(ax(1),'fontsize',14)
xlabel('Time [day]')
ylabel('Meridional distance [km]')
ylabel(colorbar,'[m]')
title('(a)','position',[2 2000.1])


ax(2)=subplot(2,1,2);
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


%---fig 12 time seires of TFS in both simulations after doubling wind
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

%======== Appendix
%---- experiment 1---------
load acc400mmtexp1.mat
t=1:10;
figure
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,2,1);
ax(1).Position=[0.1 0.12 0.38 0.8];% [left bottom width height]
plot(t,tendts,'LineWidth',2)
hold on
plot(t,advts,'LineWidth',2)
plot(t,corits,'LineWidth',2)
% plot(t,dissts,'LineWidth',2)
plot(t,dissts+vsflxtotts,'LineWidth',2)
plot(t,dphits,'LineWidth',2)
plot(t,extts,'LineWidth',2)
% plot(t,vsflxtotts,'LineWidth',2)
plot(t,rests,'k','LineWidth',2)
grid on
xlabel('Time [year]')
ylabel('[m^{4}s^{-2}]')
xlim([1 10])
legend({'tendency','advection','Coriolis','friction','TFS','wind stress','residual'},'fontsize',11,'Location','best')
set(ax(1),'fontsize',14,'Layer','top')
title('(a)','position',[2 1e9+1])


ax(2)=subplot(2,2,2);
ax(2).Position=[0.57 0.6 0.4 0.31];% [left bottom width height]
z1=9;
plot(-tend(1:z1),zc(1:z1),'LineWidth',2)
hold on
plot(adv(1:z1),zc(1:z1),'LineWidth',2)
plot(cori(1:z1),zc(1:z1),'LineWidth',2)
% plot(diss(1:z1),zc(1:z1),'LineWidth',2)
plot(diss(1:z1)+vsflxtot(1:z1),zc(1:z1),'LineWidth',2)
plot(dphi(1:z1),zc(1:z1),'LineWidth',2)
plot(ext(1:z1),zc(1:z1),'LineWidth',2)
% plot(vsflxtot(1:z1),zc(1:z1),'LineWidth',2)
plot(res(1:z1),zc(1:z1),'k','LineWidth',2)
grid on
xlim([-8e8 8e8])
ylim([zc(z1) 0])
ylabel('Depth [m]')
set(ax(2),'fontsize',14)
title('(b)','position',[-7.8e8 0.1])

ax(3)=subplot(2,2,4);
ax(3).Position=[0.57 0.12 0.4 0.31];
plot(-tend(z1+1:end),zc(z1+1:end),'LineWidth',2)
hold on
plot(adv(z1+1:end),zc(z1+1:end),'LineWidth',2)
plot(cori(z1+1:end),zc(z1+1:end),'LineWidth',2)
% plot(diss(z1+1:end),zc(z1+1:end),'LineWidth',2)
plot(diss(z1+1:end)+vsflxtot(z1+1:end),zc(z1+1:end),'LineWidth',2)
plot(dphi(z1+1:end),zc(z1+1:end),'LineWidth',2)
plot(ext(z1+1:end),zc(z1+1:end),'LineWidth',2)
% plot(vsflxtot(z1+1:end),zc(z1+1:end),'c','LineWidth',2)
plot(res(z1+1:end),zc(z1+1:end),'k','LineWidth',2)
xlim([-8e8 8e8])
ylim([zc(end) zc(z1+1)])
grid on
title('(c)','position',[-7.8e8 zc(z1+1)+0.1])
ylabel('Depth [m]')
xlabel('zonal momentum termmmtexp1.pngs [m^{3}s^{-2}]')
set(ax(3),'fontsize',14)
print -dpng ~/Desktop/mmtexp1.png

load newexp2zonalu.mat
load newexp.mat
load blue_red_saturated.mat
y1=y*1e-3;
x1=x*1e-3;

x0=10;
y0=10;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

ax(1)=subplot(2,2,1);
imagesc(x1,y1,squeeze(sst(:,:,1))'),axis xy
colorbar
xlabel('Zonal distance [km]')
ylabel('Meridional distance [km]')
set(ax(1),'fontsize',14,'layer','top')
ylabel(colorbar,'[^{o}C]','fontsize',14)
title('(a)','position',[50 4001])

ax(2)=subplot(2,2,2);
yled=1000:100:1300;
zled1=-2200;
zled2=-2400;


tc=1:7;
[Y,Z]=meshgrid(y1,zc);

contour(Y,Z,Tintmean',tc,'LineWidth',2,'ShowText','on')
set(ax(2),'XTick',[],'Ytick',[])
colormap(ax(2),[0 0 0])
h1=axes('position',get(ax(2),'position'),'color','none','fontsize',14);
hold on
contour(Y,Z,Tfinalmean',tc,'LineWidth',2,'ShowText','on')
colormap(h1,[1 0 0])
xlabel('Meridional distance [km]')
ylabel('Depth [m]')
title('zonally averaged isotherm')
title('(b)','position',[0.1e3 0])
set(gca,'fontsize',14)
hold on
plot(yled,zled1*ones(size(yled)),'color',[0,0,0],'LineWidth',2)
plot(yled,zled2*ones(size(yled)),'color',[1,0,0],'LineWidth',2)
text(1400,zled1,'month 1','fontsize',12)
text(1400,zled2,'year 30','fontsize',12)


ax(3)=subplot(2,2,3);
cbar=get(gca,'colororder');
t=31:90;
plot(t,tendts,'LineWidth',2)
hold on
plot(t,advts,'LineWidth',2)
plot(t,corits,'LineWidth',2)
% plot(t,dissts,'LineWidth',2)
plot(t,dissts+vsflxtotts,'LineWidth',2)
plot(t,dphits,'LineWidth',2)
plot(t,extts,'LineWidth',2)
% plot(t,vsflxtotts,'LineWidth',2)
plot(t,rests,'k','LineWidth',2)
plot(t,mean(corits)*ones(size(t)),'--','color',cbar(3,:))
plot(t,mean(dphits)*ones(size(t)),'--','color',cbar(5,:))
grid on
xlabel('Time [year]')
ylabel('[m^{4}s^{-2}]')
xlim([31 90])
legend({'tendency','advection','Coriolis','friction','TFS','wind stress','residual'...
    'time-mean Coriolis','time-mean TFS'},'fontsize',10,'Location','bestoutside')
set(ax(3),'fontsize',14,'layer','top')
title('(c)','position',[45 1e9+1])



ax(4)=subplot(2,2,4);
[Y,Z]=meshgrid(y1(2:end),zc(2:end-1));

pcolor(Y,Z,squeeze(mean(ep,1,'omitnan'))')
shading flat
colorbar
caxis([-5e-4 5e-4])
colormap(ax(4),map)
set(ax(4),'xtick',[],'ytick',[],'Layer','top')
ylabel(colorbar,'[m^{2}s^{-2}]','fontsize',14)
% title('EP flux [m^{2}s^{-2}]')

% [Y1,Z1]=meshgrid(y,zc);
uc=-0.05:0.01:0.05;
h1=axes('position',get(ax(4),'position'),'color','none','fontsize',14);
hold on
contour(Y,Z,umeanyz(2:end,2:end-1)',uc,'LineWidth',2,'Showtext','on')
colormap(h1,[0 0 0])
xlabel('Meridional distance [km]')
ylabel('Depth [m]')
set(ax(4),'fontsize',14)
title('(d)','position',[10 1])

print -dpng ~/Desktop/exp2.png

