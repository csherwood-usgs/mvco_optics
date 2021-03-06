% p_fit_periods - Average profiles for selected two-day periods
% Refactored from p_times_mo_new
close all
clear all

ifnorm = 0; % divide all of the values by the global mean
Ca_norm = 0; % normalize profiles by Ca
lglg = 0; % make log-log plots
load ustar_av % made by plot_ustar.m
burst_length = 3480/(3600*24);
dus = [ustar_av.dn]+burst_length/2.;
speed = sqrt(ustar_av.u.^2+ustar_av.v.^2);
us = [ustar_av.us];
ew = sign([ustar_av.u]);
ewus = ew.*us;
ewspeed = ew.*speed;

% 20-minute bin averages
% load suspsed_ba20_20rstrim_crs_cen % (missing small LISST bins)
load suspsed_ba20_20rstrim_crs_cen_june_25_2018.mat
% use bin centers when available
inst = {...
   'LISSTtc','LISSTzc','LISSTvtot','LISST Tot. Vol.';...
   'LISSTtc','LISSTzc','LISSTfinesv','LISST Fines';...
   'LISSTtc','LISSTzc','LISSTmicrov','LISST Micro. Vol.';...
   'LISSTtc','LISSTzc','LISSTmacrov','LISST Macro. Vol.';...
   'LISSTtc','LISSTzc','LISSTattn','LISST Attenuation';...
   'YSIdn','YSIz','YSIturb','YSI Turbidity';... % no bin centers?
   'UMEtc','UMEzc','UMEattn650','c_{p}_{650}';...
   'UMEtc','UMEzc','UMEchl','Chl';...
   'UMEtc','UMEzc','UMEbs650','bs_{650}';...
   'trandn','tranz','tranattn','Trans. Attn.';...
   'ADVdn','OBSz','OBS','OBS Backscatter';...
   'ADVdn','ADVz','ADVagc','ADV AGC';...
   'absstc','absszc','abss1','1 MHz ABSS';...
   'absstc','absszc','abss2','2.5 MHz ABSS';...
   'absstc','absszc','abss3','4 MHz ABSS';...
   'UMEtc','UMEzc','UMEgamma','gamma';...
   };

za = 0.1; % (m) Ca estimates standardized to this elevation
nav = 0;  % number of profiles to left and right to be averaged
%list of target year days

%% Changes from Emmanuel's version in this section make all of the plots for
% one period end up on the same plot

% do all of them
%  ydt = [261:(1/3)/24:288];
%  plotdir = 'new'

ydt = [261:(1/3)/24:263];
plotdir = 'maria'
%
% ydt = [268:(1/3)/24:278];
% plotdir = 'spgtides'

% ydt = [274:(1/3)/24:276];
% plotdir = 'ophelia'
%
% ydt = [280:(1/3)/24:282];
% plotdir = 'calm'
%
% ydt = [286:(1/3)/24:288];
% plotdir = 'noreaster'
%%
yday_off = datenum('1-Jan-2011 00:00:00');
dnt = yday_off+ydt;
datestr(dnt);

% info box
%subplot(4,3,3)
xlim([0 1])
ylim([0 1])
speed = interp1(dus,ewspeed,dnt,'nearest'); speed_mn=nanmean(abs(speed)), speed_sd=nanstd(speed)
%
ustar = interp1(dus,us,dnt,'nearest'); ustar_mn=nanmean(ustar), ustar_sd = nanstd(ustar)
ubr = interp1(dus,ustar_av.ubr,dnt,'nearest'); ubr_mn=nanmean(ubr), ubr_sd = nanstd(ubr)
zo = interp1(dus,ustar_av.zoa,dnt,'nearest'); zo_mn=nanmean(zo), zo_sd=nanstd(zo)

figure(1); clf;
text(0.1,.9,datestr(nanmean(dnt)))
ydayc = nanmean(dnt)-datenum('1-Jan-2011 00:00:00');
text(0.1,.8,sprintf('Day: %6.2f',ydayc))
text(0.1,.7,sprintf('Speed: %5.2f (%6.3f) m/s',speed_mn,speed_sd))
text(0.1,.6,sprintf('ubr: % 4.3f (%5.3f) m/s',ubr_mn,ubr_sd))
text(0.1,.5,sprintf('u*c: % 4.3f (%5.3f) m/s',ustar_mn,ustar_sd))
text(0.1,.4,sprintf('zo : % 6.4f (%6.4f) m',zo_mn,zo_sd))
set(gca,'Visible','off')


%% Cp650
figure(2); hold on
plotnum=1;
ino=[7]
disp('UME Cp650')
allt=ba.(inst{ino,1});
allz=ba.(inst{ino,2});
allc=ba.(inst{ino,3});
alld=ba.UMEbs650;
allvc=vv.(inst{ino,3});
allvd=vv.UMEbs650;


% find target time
i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));

c = nanmean(allc(:,i)')';
d = nanmean(alld(:,i)')'*40;
varc = allvc(:,i);
vard = allvd(:,i)*40;
var = allvc(:,i);
z = nanmean(allz(:,i)')';
t = nanmean(allt(:,i));
%fprintf(1,'Target time: %s; mean profile time: %s\n',datestr(dnt(ii)),datestr(t));
%ustar = interp1(dus,us,t,'nearest')
%%
ok = (~isnan(c+z));

if(sum(ok)>3)
   subplot(2,3,1)
   hold on
   plot_snippet_new
   xlabel('c_p(650) and 40xb_{bp}(650) m^{-1} )');
%    ttext1 = sprintf('ws_cp= %4.2f  ws_bbp=%4.2f',-1000*0.41*pfnl.p*ustar,-1000*0.41*pfnl2.p*ustar)
%    title(ttext1);
   %xlabel('Attenuation ( m^{-1} )');
   %      ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl2.p*ustar,pfnl2.r2)
   %       hold off
end

%% UME gamma
plotnum=2;
ino=[16]
disp('UME gamma')
allt=ba.(inst{ino,1});
allz=ba.(inst{ino,2});
allc=ba.(inst{ino,3});
alld=ba.UMEbs532;
alle=ba.UMEbs650;

% find target time
i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));
c = nanmean(allc(:,i)')';
d = nanmean(alld(:,i)')';
e = nanmean(alle(:,i)')';
gamma_bb=(650/532)*log(d./e);
z = nanmean(allz(:,i)')';
t = nanmean(allt(:,i)')';
tmin = nanmin(nanmin(allt(:,i)));
tmax = nanmax(nanmax(allt(:,i)));
ok = (~isnan(c+z));

if(sum(ok)>3)
   datestr(t)
   subplot(2,3,2)
   hold on
   %plot_snippet
   hp=plot(c(ok),z(ok),'ok');
   set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
   hold on
   hp2=plot(gamma_bb(ok)+0.2,z(ok),'ok');
   set(hp2,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0])
   xlabel('\gamma_{cp}(R), \gamma_{bbp}(K)+0.3');
   ylabel('Elevation [m]')
   title(plotdir)
   %          hold off
   ylim([0.1 2])
   %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
   %title(ttext)
end


%% LISST D50
plotnum=3;
ino=[1]
disp('LISST D50')
allt=ba.(inst{ino,1});
allz=ba.(inst{ino,2});
allc=ba.LISSTD50a;

% find target time
i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));
c = nanmean(allc(:,i)')';
z = nanmean(allz(:,i)')';
t = nanmean(allt(:,i)')';
ok = (~isnan(c+z));

if(sum(ok)>3)
   datestr(t)
   %         datestr(tmin)
   %         datestr(tmax)
   subplot(2,3,3)
   hold on
   %plot_snippet
   hp=plot(c(ok),z(ok),'ok');
   ylabel('Elevation [m]')
   xlabel('Lisst D_{50} [\mum]');
   %       ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
   %       title(ttext)
end
ylim([0.1 2])

%% Backscattering ratio
plotnum=4;
ino=[8]
disp('UME bbr')
plotnum=plotnum+1;
allt=ba.(inst{ino,1});
allz=ba.(inst{ino,2});
allc=(ba.UMEbs650./ba.UMEattn650);

%       if(ino==8 || ino==16) % no variance for Chl or gamma
%          allvc=ones(size(allc));
%       else
%          allvc=vv.(inst{ino,3});
%       end

% find target time
i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));
c = nanmean(allc(:,i)')';
z = nanmean(allz(:,i)')';
t = nanmean(allt(:,i)')';
ok = (~isnan(c+z));

if(sum(ok)>3)
   datestr(t)
   subplot(2,3,4)
   hold on
   hp=plot(c(ok),z(ok),'ok');
   ylabel('Elevation [m]')
   xlabel('b_{bp}/b_p(650)');
   %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
   %title(ttext)
   ylim([0.1 2])
end
%% UME chlorophyll over attenuation
plotnum=5;
for ino=[8]
   disp('UME chlorophyll over attenuation')
   plotnum=plotnum+1;
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=(ba.(inst{ino,3})./(ba.(inst{7,3})));
   
   %       if(ino==8 || ino==16) % no variance for Chl or gamma
   %          allvc=ones(size(allc));
   %       else
   %          allvc=vv.(inst{ino,3});
   %       end
   
   % find target time
   i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));
   c = nanmean(allc(:,i)')';
   z = nanmean(allz(:,i)')';
   t = nanmean(allt(:,i)')';
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,5)
      hold on
      hp=plot(c(ok),z(ok),'ok');
      ylabel('Elevation [m]')
      xlabel('Chl/c_p(650)');
      %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
      %title(ttext)
      ylim([0.1 2])
   end
end
%% LISST density - special case
plotnum=6;
ino=[1]
disp('LISST density')
allt=ba.(inst{ino,1});
allz=ba.(inst{ino,2});
allc=ba.(inst{5,3})./ba.(inst{1,3});

% find target time
i = find(allt(6,:)>=dnt(1) & allt(6,:)<dnt(end));
c = nanmean(allc(:,i)')';
z = nanmean(allz(:,i)')';
t = nanmean(allt(:,i)')';

ok = (~isnan(c+z));

if(sum(ok)>3)
   datestr(t)
   subplot(2,3,6)
   hold on
   %plot_snippet
   hp=plot(c(ok),z(ok),'ok');
   ylabel('Elevation [m]')
   xlabel('LISST Density Proxy');
   ylim([0.1 2])
   %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
   %title(ttext)
   
end
%%

shg
pfn = sprintf('./%s.png',plotdir)
print(pfn,'-dpng')
