clear

flist = {'mor01','mor01_a35b05','mor01_a35b34','mor02','mor02_a35b05','mor02_a35b34',...
         'mor05','mor05_a35b05','mor05_a35b34','mor1','mor1_a35b05','mor1_a35b34'};
for ii=4%:length(flist)
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
nc = ncgeodataset([url,flist{ii},'.nc'])
%%
run_name=flist{ii};
% read vertical grid parameters
Vtransform = nc{'Vtransform'}(:);
Vstretching = nc{'Vstretching'}(:);

s_rho = nc{'s_rho'}(:);
s_w = nc{'s_w'}(:);
Cs_r = nc{'Cs_r'}(:);
Cs_w = nc{'Cs_w'}(:);
N = length(s_rho);
Np = length(s_w);

theta_s = nc{'theta_s'}(:)
theta_b = nc{'theta_b'}(:)
depth_c = nc{'hc'}(:)

a=theta_s;
b=theta_b;
sr = s_rho;
C = (1-b)*sinh(a*sr)/sinh(a) + b*[tanh(a*(sr+0.5))/(2*tanh(0.5*a)) - 0.5];

%% read water depth
h = nc{'h'}(3,4);
hc = nc{'hc'}(:);
zeta = nc{'zeta'}(:,3,4);

z=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                1, h, 0,0))';
zw=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                 5, h, 0,0))';
             
time = nc{'ocean_time'}(:);
nt = length(time);
nz = length(z);
nzw = length(zw);
dz = diff(zw)

dz2d = repmat(dz,[nt,1]);
%%
fdiam = 1e6*nc{'Sd50'}(1)
ws = 1e3*nc{'Wsed'}(1)

% check mass conservation of suspended NCS classes only
% (last class in these runs is sand)
clear mud
ncs = length(ws)-1;
for n=1:ncs
   ncname = sprintf('mud_%02d',n)
   mud(n,:,:)=squeeze(nc{ncname}(:,:,3,4));
end
muds = squeeze(sum(mud));
mmud = muds.*dz2d;

initial_mass = sum(mmud(1,:));
final_mass = sum(mmud(nt-1,:));
summud_ts = sum(mmud,2);
max_mud_change = max( abs( summud_ts-summud_ts(1) ))

%%
% make 2D arrays of time and depths
t2d = repmat(time,[1,nz]);
tw2d = repmat(time,[1,nz+1]);
z2d = repmat(h+z,[nt,1]);
zw2d = repmat(h+zw,[nt,1]);

%%
% calculate size- and settling-velocity weighted floc dists
ws_av = squeeze(sum((mud+eps).*repmat(ws(1:ncs),[1,nt,nz]),1)./sum((mud+eps),1));
fdiam_av = squeeze(sum((mud+eps).*repmat(fdiam(1:ncs),[1,nt,nz]),1)./sum((mud+eps),1));

mconc = sum(mud,3);


%%

clear ho
%%Calculate and plot turbulence and mixing info
tke = nc{'tke'}(:,:,3,4);
gls = nc{'gls'}(:,:,3,4);
akv_bak = nc{'Akv_bak'}(1);
akt_bak = nc{'Akt_bak'}(1);
akv = nc{'AKv'}(:,:,3,4);
nueau =      1.5e-6;
gls_p =     -1.000;  %gls_p           GLS stability exponent.
gls_m =      0.500;  %gls_m           GLS turbulent kinetic energy exponent.
gls_n =     -1.000;  %gls_n           GLS turbulent length scale exponent.
gls_cmu0 =   5.4770e-01; %            GLS stability coefficient.
 
exp1 = 3.0+gls_p/gls_n;
exp2 = 1.5+gls_m/gls_n;
exp3 = -1.0/gls_n;
diss = gls_cmu0^exp1.*tke.^exp2.*gls.^exp3;
Gval=sqrt(diss/nueau);


%% Calculate G on rho 
Gc = 0.5*(Gval(:,1:end-1)+Gval(:,2:end));
if size(ws_av,1)>=11
    nts = 10;
else
    nts = size(ws_av,1)-1;
end

tconc = squeeze(sum(mud,1));

nf = 2.;
m = .2;
q = (nf-1.)/(2.*m);
%% colormaps
cb1=colormap('pink(32)');cb2=flipud(cb1);
cb3=colormap('hsv');cb3=flipud(cb3);
cb4=colormap('copper(32)');cb4=flipud(cb4);
ccm1 = contrast(squeeze(sum(mud,1)));ccm1 = flipud(ccm1);
ccm2 = contrast(fdiam_av);
ccm3 = contrast(ws_av);ccm3=flipud(ccm3);
my_map = [.90 .90 .55;
          .85 .85 .45;
          .80 .80 .40;
          .75 .75 .35;
          .70 .65 .25;
          .65 .50 .15;
          .55 .40 .1;
          .45 .30 .05;
          .50 .25 .0;
          .40 .20 .0;
          .30 .15 .0;
          .25 .10 .0];
for i=1:ncs
    if i<6;
        cl(i,:)=[0,0,(i)*50]/255;
    elseif i>=6&&i<10
        cl(i,:)=[(i-5)*50,0,(i-5)*50]/255;
    else
        cl(i,:)=[(i-10)*50,(i-10)*50,0]/255;
    end
    if ncs>16
        if i<7;
            cl(i,:)=[0,0,(i)*40]/255;
        elseif i>=7&&i<14
            cl(i,:)=[0,(i-5)*30,(i-5)*30]/255;
        else
            cl(i,:)=[(i-10)*16,(i-10)*16,0]/255;
        end
    end 
end
%%     Figure creation
figure(1);clf
set(gcf,'PaperPosition',[.5,.5,10,8]);
subplot ('Position',[0.1100    0.7093    0.2012    0.2157])
ab=get(gca,'Position');
set(gca,'FontSize',11,'FontWeight','bold')
dum=squeeze(sum(mud,1));
pcolor(t2d/3600.,z2d,dum);
cran=[median(dum(:))-2*iqr(dum(:)),median(dum(:))+2*iqr(dum(:))];
if cran(1)<0;cran(1)=0;end
caxis(cran);

shading interp;set(gca,'XTick',[0:2:12])
title(['(a) Conc. (kg m^{-3})'],'FontSize',13,'FontWeight','bold')
colormap(gca,my_map); 
hcb=colorbar; set(hcb,'Position',[0.3224    0.7093    0.019    0.2157],'FontSize',12,'FontWeight','bold')
set(gca,'Position',ab)
subplot('Position',[0.1100    0.4096    0.2012    0.2157])
ab=get(gca,'Position');
set(gca,'FontSize',11,'FontWeight','bold')
pcolor(t2d/3600.,z2d,fdiam_av);
shading interp;set(gca,'XTick',[0:2:12])
ylabel('Elevation (m)')
title('(b)  Diameter (\mum)','FontSize',13,'FontWeight','bold')
colormap(gca,ccm2); 
hcb=colorbar;set(hcb,'Position',[0.3224    0.4096    0.019    0.2157],'FontSize',12,'FontWeight','bold'); 
set(gca,'Position',ab)
subplot('Position',[0.1100    0.1100    0.2012    0.2157])
ab=get(gca,'Position');
set(gca,'FontSize',11,'FontWeight','bold')
pcolor(t2d/3600.,z2d,ws_av);
shading interp;set(gca,'XTick',[0:2:12]);
cran=[median(ws_av(:))-2*iqr(ws_av(:)),median(ws_av(:))+2*iqr(ws_av(:))];
if cran(1)<0;cran(1)=0;end
caxis(cran);
colormap(gca,cb2); 
hcb=colorbar; set(hcb,'Position',[0.3224    0.1100    0.019    0.2157],'FontSize',12,'FontWeight','bold');
xlabel('Time (h)');
title('(c)  Settling velocity (mm s^{-1})','FontSize',13,'FontWeight','bold')
set(gca,'Position',ab)
subplot ('Position',[0.4208    0.7093    0.2134    0.2157])
set(gca,'FontSize',11,'FontWeight','bold')
line(tconc(end,:),h+z,'Color','k','LineWidth',1.5); %axtt
set(gca,'YLim',[0,12])
title(['(d) Total Conc. (kg m^{-3})'],'FontSize',13,'FontWeight','bold')
subplot ('Position',[0.4208    0.4096    0.2134    0.2157])
set(gca,'FontSize',11,'FontWeight','bold')
line(fdiam_av(end,:),h+z,'Color','k','LineWidth',1.5);%axtt
set(gca,'YLim',[0,12])
title('(e) Diameter (\mum)','FontSize',13,'FontWeight','bold')
subplot ('Position',[0.4208    0.1100    0.2134    0.2157])
set(gca,'FontSize',11,'FontWeight','bold')
line(ws_av(end,:),h+z,'Color','k','LineWidth',1.5);%axtt
set(gca,'YLim',[0,12])
title('(f) Avg. Ws (mm s^{-1})','FontSize',13,'FontWeight','bold')
subplot (1,3,3)
set(gca,'FontSize',11,'FontWeight','bold')
for i=1:ncs
    ho(i)=semilogx(squeeze(mud(i,nt,:))+eps,h+z,'Color',cl(i,:),'LineWidth',1.5); hold on
end
ho(i+1)=semilogx(squeeze(sum(mud(:,nt-1,:),1)),h+z,'Color','k','LineStyle','--','LineWidth',2);
set(gca,'XLim',[1e-5,3e0],'XTick',[1e-5,1e-3,1e-1,1])
hle=legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tot.'],'Location','BestOutside');
set(hle,'Position',[0.8965    0.2815    0.0558    0.6240]);
title('(g)  Concentration (kg m^{-3})','FontSize',13,'FontWeight','bold')
ylabel('Elevation (m)')
xlabel('Concentration (kg m^{-3})')
% Figure 1 output
eval(['print -dpng Fig4_',run_name,'_redone.png'])
eval(['print -depsc -zbuffer Fig10_',run_name,'_redone.eps'])

%% Figure 2 plots of response
% load the alpha_v
a=load('alpha_v_nf2.txt')
alpha_v(i) 
figure(2); clf
subplot(4,2,1)

end


