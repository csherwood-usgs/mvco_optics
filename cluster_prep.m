% cluster_prep.m
% make matrices suitable for k-means clustering

% Script to plot profiles. Version 3 was made by EB and the plotting was
% modified by CRS. This version is changed as follow:
% Parameters to be plotted are extracted from the database with names at
% the top. This is a waste of memory, but makes the plotting code easier to
% read and understand.
close all
clear all

load ustar_av % made by plot_ustar.m
burst_length = 3480/(3600*24);
dus = [ustar_av.dn]+burst_length/2.;
speed = sqrt(ustar_av.u.^2+ustar_av.v.^2);
us = [ustar_av.us];
us_wave = [ustar_av.us];
us_wc=[ustar_av.ustrr];
ubr = [ustar_av.ubr];
zoa = [ustar_av.zoa];
ew = sign([ustar_av.u]);
ewus = ew.*us;
ewspeed = ew.*speed;

% 20-minute bin averages
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
% a few of these are negative. Not sure why.
ineg = ba.UMEbs650 <0.;
ba.UMEbs650(ineg)=1e-8;
%% UMEattn650 profiles
allt=ba.UMEtc;
allz=ba.UMEzc;
allc=ba.UMEattn650;
alld=ba.UMEbs532;
alle=ba.UMEbs650;
t1 = allt(3:9,1:2031);
z1 = allz(3:9,1:2031);
c1 = allc(3:9,1:2031); %attn650
c4 = alld(3:9,1:2031); %bs532
c5 = alle(3:9,1:2031); %bs650

allf = ba.UMEchl;
c6= allf(3:9,1:2031); % Chl

% LISST D50v
allt=ba.LISSTtc;
allz=ba.LISSTzc;
allc=ba.LISSTD50v;
t2 = allt(3:9,1:2031);
z2 = allz(3:9,1:2031);
c2 = allc(3:9,1:2031);

% abss
allt=ba.absstc;
allz=ba.absszc;
allc=ba.abss2;
t3 = allt(3:9,1:2031);
z3 = allz(3:9,1:2031);
c3 = allc(3:9,1:2031);

good_columns = find(~isnan(sum([t1;t2;t3;z1;z2;z3;c1;c2;c3;c4;c5;c6])));
c1 = c1(:,good_columns);
c2 = c2(:,good_columns);
c3 = c3(:,good_columns);
c4 = c4(:,good_columns);
c5 = c5(:,good_columns);
c6 = c6(:,good_columns);
c7 =  (650/532)*log(c4./c5); % gammabb


% calculate mean depths and mean times
t = mean([t1(:,good_columns);t2(:,good_columns);t3(:,good_columns)],1);
z = mean([z1(:,good_columns),z2(:,good_columns),z3(:,good_columns)],2);
tstd = std([t1(:,good_columns);t2(:,good_columns);t3(:,good_columns)],0,1);
zstd = std([z1(:,good_columns),z2(:,good_columns),z3(:,good_columns)],0,2);
size(c1)
size(c2)
size(c3)
size(z)
size(t)

%% standardize each array
c1mean = mean(c1(:));
c2mean = mean(c2(:));
c3mean = mean(c3(:));
c4mean = mean(c4(:));
c5mean = mean(c5(:));
c6mean = mean(c6(:));
c7mean = mean(c7(:));

c1std = std(c1(:));
c2std = std(c2(:));
c3std = std(c3(:));
c4std = std(c4(:));
c5std = std(c5(:));
c6std = std(c6(:));
c7std = std(c7(:));


c1s = (c1-c1mean)/c1std;
c2s = (c2-c2mean)/c2std;
c3s = (c3-c3mean)/c3std;
c4s = (c4-c4mean)/c4std;
c5s = (c5-c5mean)/c5std;
c6s = (c6-c6mean)/c6std;
c7s = (c7-c7mean)/c7std;

save profiles.mat t tstd z zstd c1mean c1std c2mean c2std c3mean c3std c1s c2s c3s

% interpolate
uwi = interp1(dus,us_wave,t);
usi = interp1(dus,ewus,t);

%% cp605, d50, abss
c123s = [c1s;c2s;c3s]';
nk = 5
clr = linspecer(nk);
idx = kmeans(c123s,nk);
figure(1); clf
ip=1
for k=1:nk
   subplot(nk,3,ip)
   plot( mean(c1(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c1(:,find(idx==k)),2)+std(c1(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c1(:,find(idx==k)),2)-std(c1(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 11])
   if(k==1),title('cp650'),end
   
   subplot(nk,3,ip+1)
   plot( mean(c2(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c2(:,find(idx==k)),2)+std(c2(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c2(:,find(idx==k)),2)-std(c2(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([100 350])
   if(k==1),title('LISST d50'),end
   
   subplot(nk,3,ip+2)
   if(k==1),title('abss 2.5'),end
   plot( 1000*mean(c3(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( 1000*mean(c3(:,find(idx==k)),2)+1000*std(c3(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( 1000*mean(c3(:,find(idx==k)),2)-1000*std(c3(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 5])
   if(k==1),title('abss 2.5'),end
   ip=ip+3;
end
print -dpng -r300 optics_acoustics_classes.png
%
figure(2)
for k=1:nk
   h=scatter(t(find(idx==k)),uwi(find(idx==k)),12,'filled');
   set(h,'markerfacecolor',clr(k,:),'markeredgecolor',clr(k,:))
   hold on
end
datetick('x','keeplimits')
xlabel('2011')
ylabel('u*cw (m/s)')
print -dpng -r300 optics_acoustics_time_series.png


%% all optics
c12456s = [c1s;c2s;c4s;c5s;c6s]';
nk = 5
clr = linspecer(nk);

idx = kmeans(c12456s,nk);
%
figure(3); clf
ip=1
for k=1:nk
   
   subplot(nk,5,ip)
   plot( mean(c1(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c1(:,find(idx==k)),2)+std(c1(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c1(:,find(idx==k)),2)-std(c1(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 11])
   if(k==1),title('cp650'),end
   
   subplot(nk,5,ip+1)
   plot( mean(c2(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c2(:,find(idx==k)),2)+std(c2(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c2(:,find(idx==k)),2)-std(c2(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([100 350])
   if(k==1),title('LISST d50'),end
   
   subplot(nk,5,ip+2)
   plot( mean(c4(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c4(:,find(idx==k)),2)+std(c4(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c4(:,find(idx==k)),2)-std(c4(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 .3])
   if(k==1),title('bs 532'),end
   
   subplot(nk,5,ip+3)
   plot( mean(c5(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c5(:,find(idx==k)),2)+std(c5(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c5(:,find(idx==k)),2)-std(c5(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 .3])
   if(k==1),title('bs 650'),end
   
      subplot(nk,5,ip+4)
   plot( mean(c6(:,find(idx==k)),2),z, '-','linewidth',3,'color',clr(k,:))
   hold on
   plot( mean(c6(:,find(idx==k)),2)+std(c6(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   plot( mean(c6(:,find(idx==k)),2)-std(c6(:,find(idx==k)),0,2),...
      z, '--','linewidth',2,'color',clr(k,:))
   xlim([0 20])
   if(k==1),title('chl'),end
   
   ip=ip+5;

end
print -dpng -r300 optics_classes.png


figure(4)
for k=1:nk
   h=scatter(t(find(idx==k)),uwi(find(idx==k)),12,'filled');
   set(h,'markerfacecolor',clr(k,:),'markeredgecolor',clr(k,:))
   hold on
end
datetick('x','keeplimits')
xlabel('2011')
ylabel('u*cw (m/s)')
print -dpng -r300 optics_time_series.png