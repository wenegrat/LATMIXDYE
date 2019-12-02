% Movie of Float in SI Boundary layer & Gradient
clear

% FLOAT

load '../Mar05_SI_2_Track.mat'
% F. lons lats yds
F.yds=F.yds-1;  % to Craig date
efile='/Users/eric/Documents/LateralMixing/2012 Experiment/Operations/Float Deployments/Mar05 SI_2/DATA/Env.mat';
E=load(efile);
E.yd=E.yd-1;
E.lon=interp1(F.yds,F.lons,E.yd);
E.lat=interp1(F.yds,F.lats,E.yd);

% Fake until I get float 
% E.yd=F.yds;
% E.lon=F.lons;
% E.lat=F.lats;
% E.P=F.Ps;
%   %%%%%

%MET DATA
uway=load('wind_speed_dir.mat');
% tmet        wind_dir    wind_speed
uway.yday=uway.tmet;

met.t=uway.yday;
g=find(uway.yday>E.yd(1) & uway.yday<E.yd(end) );
met.t=uway.yday(g);
met.wind_speed=uway.wind_speed(g);
met.wind_dir=uway.wind_dir(g);

% Load Triaxus
root='/Volumes/triaxus/LatMix_2012/SI_2a/';
D=dir([root 'CTD*.mat']);
yd=[];Sig0=[];lat=[];lon=[];X=[];S=[];T=[];Fl=[];P=[];depth=[];
for i=1:length(D)
    file=[root D(i).name];
    load(file,'ts');
    disp([ file(end-30:end) ' ' num2str(ts.t(1),3) ' -' num2str(ts.t(end),3) ' ' num2str(size(ts.sigma_t))]);
    Sig0=[Sig0;ts.sigma_t];
    S=[S;ts.S];
    T=[T;ts.T];
    Fl=[Fl;ts.fluor];
    yd=[yd; ts.t];
    X=[X; ts.x];
    lat=[lat; ts.lat];
    lon=[lon; ts.lon];
    depth=[depth;ts.depth];
    P=[P;sw_pres(ts.depth,36*ones(size(ts.depth)))];
end
dX=diff(X);
dX(dX<0)=0.7;
X=[0;cumsum(dX)];
ts.lon=lon;ts.lat=lat;ts.P=P;ts.x=X;ts.yd=yd;
ts.T=T;ts.S=S;ts.sigma_t=Sig0;ts.t=yd;ts.depth=depth;ts.F=Fl;
ts.F=ts.F-median(ts.F);

Dye=ts.F;
Dye(Dye<0.003)=0.;

% float data to triaxus time
lonf=interp1(E.yd,E.lon,ts.t,'linear');
latf=interp1(E.yd,E.lat,ts.t,'linear');
Pf=interp1(E.yd,E.P,ts.t);

% smooth a bit  (nyquist is 1 sec)
[B,A] = butter(2,0.00015,'low');
g=find(~isnan(lonf));
lonf(g)=filtfilt(B,A,lonf(g));
latf(g)=filtfilt(B,A,latf(g));

% Move lat/lon to time of Triaxus - later
Tlag=850/3; % seconds = cable out / speed
ts.lon=interp1(ts.yd,ts.lon,ts.yd+Tlag/86400);
ts.lat=interp1(ts.yd,ts.lat,ts.yd+Tlag/86400);

%% find cross track direction and rotate ship-float that way
lonx=lonf*111200*sin(36/180*pi);
latx=latf*111200;
dx= diff(lonx);  dx=[dx;dx(end)];
dy= diff(latx);   dy=[dy;dy(end)];
th=atan2(dx,dy);

X0=(ts.lon-lonf)*111200*sin(36/180*pi);  % ship in float coordinates
Y0=(ts.lat-latf)*111200;

X=X0.*cos(th)+Y0.*sin(th);
Y= -X0.*sin(th)  +Y0.*cos(th);

% test rotation
if 0
    dX0=10000.*ones(size(X0)); dY0=0;
    dX=  dX0.*cos(th)+dY0.*sin(th);
    dY= -dX0.*sin(th)+dY0.*cos(th);
    plot(lonx,latx,'k-');
    for i=1:1000:length(dx)
        plot([lonx(i) lonx(i)+dX(i)],[latx(i) latx(i)+dY(i)],'-');
        hold on
        axis equal
    end
end

L=5000;  % half width of plot
Zmax=120;  % depth of plot
tsm=800;  % distance smoothing
zsm=5;  % z smoothing
limit=5;
dsm=200; % dye dist smoothing
dlimit=20;

movie=1;
domat=0;
Nframe=50;
Toverlap=0.1; %0.15; % time average of each frame
Tstart=ts.t(1)+Toverlap/2;
days=linspace(Tstart,ts.t(end),Nframe);

for i=[1:Nframe] %Nframe]
    T1=max(days(i)-Toverlap,Tstart);
    T2=min(days(i)+Toverlap,ts.t(end));
    g=find(ts.t >T1 & ts.t <T2 & ts.depth<Zmax & ts.depth>0 & ~isnan(ts.S) & ts.S> 34.8 );
    g=g(1:3:end);
    gdye=find(ts.F>0.003 & ts.t >T1 & ts.t <T2 & ts.depth<Zmax & ts.depth>0 );
    if(~movie)
        figure
    end
    %disp([i length(g) T2-200]);
    a(1)=axes('position',[0.1 0.3 0.8 0.6]);
    %a(2)=axes('position',[0.5 0.3 0.4 0.6]);
    a(3)=axes('position',[0.1 0.06 0.5 0.13]);  % wind
    amap=axes('position',[0.65 0.06 0.3 0.13]);  % map
    tit=axes('position',[0.1 0.90 0.8 0.08]);

    % Sections along Y=0
    axes(a(1));  % DENSITY SALINITY
    %patch([-L L L -L],[-L -L L L],[1 1 1]*0.5); hold on;

    G=ts.sigma_t;
    gmin=24.;gmax=26.0;
    gthick=[gmin:0.25:gmax];
    icn=0;
    for ic=1:length(gthick)
        for ic2=1:4
            icn=icn+1;
            gthin(icn)=gthick(ic)+ic2*0.05;
        end
    end

    S=ts.S;
    smin=35.1;smax=36.8;
    sthick=[smin:0.25:smax];
    sthin=[smin:0.01:smax];

    y2x=zsm/tsm;
    Xgrid=[-L:100:L];
    Zgrid=[0:1:Zmax];
    [xg, zg]=meshgrid(Xgrid,Zgrid);

    % SALINITY DENSITY
    caxis([smin smax]);
    Smap=Ggrid(X(g),ts.depth(g),S(g),xg,zg,zsm,limit,y2x,0);
    % adjust X position to center map
    gcent=find(~isnan(Smap(:)));
    
    Xcent=nmean(xg(gcent));% center on section
    Xcent=0.;  % center on float

    %     [c,h]=contour(xg,-zg,Smap,sthin);hold on
    %     set(h,'LineWidth',0.5);hold on
    pcolor(xg-Xcent,zg,Smap);shading flat;hold on;
    caxis([smin smax]);

    axis([ -L L 0 Zmax]);
    hold on
    
    % test - input points
    %color_mark(X(g)-Xcent, -ts.depth(g), S(g) );

    Gmap=Ggrid(X(g),ts.depth(g),G(g),xg,zg,zsm,limit,y2x,0);
    [c,h]=contour(xg-Xcent,zg,Gmap,gthin,'k-');hold on
    set(h,'LineWidth',0.1,'Color',[1 1 1]*0.1);
    [c,h]=contour(xg-Xcent,zg,Gmap,gthick,'k-');set(h,'LineWidth',2);
    hl=clabel(c,h,'Color','k');
    fixclabel(hl );
    caxis([smin smax]);
    set(gcf,'Color','w')
    set(gca,'color',[1 1 1]*0.5);
    colorbarthin;
    
    % Dye
    dy2x=zsm/dsm;
    Dmap=Ggrid(X(g),ts.depth(g),Dye(g),xg,zg,zsm,dlimit,dy2x,0);
    [c,h]=contour(xg-Xcent,zg,Dmap,[0.003 0.01 0.03 0.1 0.3 1] ,'g-');set(h,'LineWidth',1);
    %plot(X(gdye)-Xcent,-ts.P(gdye),'g.','markerSize',3);

    % float - position is -Xcent 
    if i==1
        xf=Xcent*ones(size(g));
    else  % sort of a fake position base on trend from last offset
        xf=Xcentp+(Xcent-Xcentp)*(ts.t(g)-days(i-1))/(days(i)-days(i-1));
    end
    xf=xf-100+[1:length(xf)]'/length(xf)*200; % very slight motion to spread
    plot(-xf,Pf(g),'m-');hold on;
    gx=floor(length(xf)/2);
    plot(-xf(gx),Pf(g(gx)),'m.','LineWidth',1,'MarkerSize',2);
    plot(-xf(gx),Pf(g(gx)),'yo','LineWidth',4,'MarkerSize',2);
    plot(-xf(gx),Pf(g(gx)),'ko','LineWidth',1,'MarkerSize',7);

    hx=xlabel('Cross Frontal Distance / m');
    hy=ylabel('Depth / m');
    set([hx hy],'FontWeight','bold','FontSize',12);
    set(gca,'ydir','reverse');
%     YY=get(gca,'yticklabel');
%     YY(1:6,1)=' ';
    ht=text(L-800,Zmax-5,sprintf('%d\n%6.3f',i,days(i)));
    set(ht,'Fontsize',8,'HorizontalAlignment','Right');

    axes(a(3));
    %%  Met / Nav panel
    plot(met.t,met.wind_speed,'k-');
    hold on
    gw=find(met.t>T1 & met.t<T2);
    plot(met.t(gw),met.wind_speed(gw),'r-','LineWidth',5);
    %axis([215 216.6 0 20]);
    set(gca,'FontSize',9);
    hy= ylabel('Wind Spd./ ms^{-1}');
    set(hy,'FontWeight','bold','FontSize',10);
    aa=axis;
    %axdate(10)
    %  wind direction
    ht=text(62.6,4,sprintf('%3.0f^o',median(met.wind_dir(gw))));
    set(ht,'Color','r','FontSize',12,'FontWeight','Bold');
    
    %%  Nav panel
    axes(amap)
    plot(lon,lat,'k-')
    hold on
    plot(lon(g),lat(g),'r-','LineWidth',2);
    axis tight
    set(gca,'dataaspectratio',[1 cos(36/180*pi) 1]);
    set(gca,'FontSize',10);
    % add wind vector
    Umean= -met.wind_speed(gw).*sin(met.wind_dir(gw)/180*pi);
    Vmean= -met.wind_speed(gw).*cos(met.wind_dir(gw)/180*pi);
    Umean=nmean(Umean);
    Vmean=nmean(Vmean);
    lonmean=-64.5;latmean=39.3;
    S=0.08;
    dlat= Vmean*S;
    dlon= Umean*S/cos(36/180*pi);
    plot([lonmean lonmean+dlon],[latmean latmean+dlat],'r-');
    plot(lonmean,latmean,'ro');
    
    axes(tit);
    ht=text(0,0.8,'LATMIX: SI Turbulence');
    set(ht,'FontSize',14,'horizontalalignment','Center','FontWeight','Bold');
    ht=text(0,0.2,'Potential Density Contours    Salinity Colors');
    set(ht,'FontSize',9,'horizontalalignment','Center');
    set(gca,'visible','off','box','off');
    datetick('x',15);

    if domat
        file=['Mat/R' sprintf('%04d.mat',i)];
        Tmap=Ggrid(X(g),ts.depth(g),ts.T(g),xg,zg,sm,limit,y2x,0);
        readme='Variables in the map: Smap=salinity, Gmap=Sigma0 Tmap=T (untested)';
        save(file,'xg','zg','Smap','Xcent','Gmap','T1','T2','Tmap','readme');
    end
    if movie
        file=['Frames/R' sprintf('%04d',i)];
        print('-dpng',file);
        clf
    end
    Xcentp=Xcent;
end

