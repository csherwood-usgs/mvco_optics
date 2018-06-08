% p_fit_mo_times - Fit settling velocities for several quantities at
% specific times
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
load suspsed_ba20_20rstrim_crs_cen
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
% ydt = [261:(1/3)/24:263];
% plotdir = 'maria'
% 
% ydt = [268:(1/3)/24:270];
% plotdir = 'spgtides'

% ydt = [274:(1/3)/24:276];
% plotdir = 'ophelia'
% 
% ydt = [280:(1/3)/24:282];
% plotdir = 'calm'
% 
ydt = [286:(1/3)/24:288];
plotdir = 'noreaster'

yday_off = datenum('1-Jan-2011 00:00:00');
dnt = yday_off+ydt;
datestr(dnt);
%%
for ii=1:length(dnt)
   figure(1); clf
   %% info box
   subplot(4,3,3)
   xlim([0 1])
   ylim([0 1])
   speed = interp1(dus,ewspeed,dnt(ii),'nearest');
   ustar = interp1(dus,ewus,dnt(ii),'nearest');
   ubr = interp1(dus,ustar_av.ubr,dnt(ii),'nearest');
   zo = interp1(dus,ustar_av.zoa,dnt(ii),'nearest');
   text(0.1,.9,datestr(dnt(ii)))
   ydayc = dnt(ii)-datenum('1-Jan-2011 00:00:00');
   text(0.1,.8,sprintf('Day: %6.2f',ydayc))
   text(0.1,.7,sprintf('Speed: %5.2f m/s',speed))
   text(0.1,.6,sprintf('ubr: % 4.2f m/s',ubr))
   text(0.1,.5,sprintf('u*c: % 4.2f m/s',ustar))
   text(0.1,.4,sprintf('zo : % 6.4f m',zo))
   set(gca,'Visible','off') 
   
   %% Cp650
   plotnum=1;
   ino=[7]
   disp('UME Cp605')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});
   allvc=vv.(inst{ino,3});
   if(ifnorm)
      allc = allc/(nanmean(allc(:)));
   end
   
   % find target time
   i = find(allt(6,:)>=dnt(ii),1,'first');
   
   c = allc(:,i);
   varc = allvc(:,i);
   z = allz(:,i);
   t = nanmean(allt(:,i));
   fprintf(1,'Target time: %s; mean profile time: %s\n',datestr(dnt(ii)),datestr(t));
   ustar = interp1(dus,us,t,'nearest')
   
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      subplot(4,3,1)
      plot_snippet
      xlabel('Attenuation ( m^{-1} )');
      ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2)
      title(ttext);
   end
   
   %% LISST attenuation
   plotnum=2;
   ino=5;
   disp('LISST Attenuation')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});
   
   if(ifnorm)
      allc = allc/(nanmean(allc(:)));
   end
   
   % find target time
   i = find(allt(6,:)>dnt(ii),1,'first');
   fprintf(1,'Target time: %s; mean profile time: %s\n',datestr(dnt(ii)),datestr(t));
   % get closest u* value
   ustar = interp1(dus,us,t,'nearest')
   c = allc(:,i);
   varc = allvc(:,i);
   z = allz(:,i);
   t = nanmean(allt(:,i));
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      subplot(4,3,plotnum)
      plot_snippet
      xlabel('Attenuation ( m^{-1} )');
      ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2)
      title(ttext);
   end
   
   %% LISST Vol. conc - micro and macro flocs
   plotnum=3;
   for ino=[3:4]
      plotnum = plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=ba.(inst{ino,3});
      
      if(ifnorm)
         allc = allc/(nanmean(allc(:)));
      end
      
      % find target time
      i = find(allt(6,:)>dnt(ii),1,'first');
      
      c = allc(:,i);
      varc = allvc(:,i);
      z = allz(:,i);
      t = nanmean(allt(:,i))
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         %         datestr(tmin)
         %         datestr(tmax)
         subplot(4,3,plotnum)
         plot_snippet
         xlabel('Volume Concentration (m^3/m^3)');
         ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         title(ttext)
      end
   end
   
   %% LISST Vol. conc total
   plotnum=6;
   ino=[1]
   disp('LISST volume concentration')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});
   
   if(ifnorm)
      allc = allc/(nanmean(allc(:)));
   end
   
   % find target time
   i = find(allt(6,:)>dnt(ii),1,'first');
   
   c = allc(:,i);
   allvc=vv.(inst{ino,3});
   varc = allvc(:,i);
   z = allz(:,i);
   disp([z,c,varc])
   c(c<1e-8)=NaN;
   
   t = nanmean(allt(:,i));
   tmin = nanmin(nanmin(allt(:,i)));
   tmax = nanmax(nanmax(allt(:,i)));
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      %         datestr(tmin)
      %         datestr(tmax)
      subplot(4,3,plotnum)
      plot_snippet
      xlabel('Volume Concentration (m^3/m^3)');
      ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
      title(ttext)
   end
   
   
   %% ABSS profiles
   plotnum=6;
   for ino=[13:15]
      disp('ABSS')
      plotnum = plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=ba.(inst{ino,3});
      
      % find target time
      i = find(allt(6,:)>dnt(ii),1,'first');
      
      c = allc(:,i);
      varc = allvc(:,i);
      z = allz(:,i);
      t = nanmean(allt(:,i));
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(4,3,plotnum)
         plot_snippet
         xlabel('Backscatter Intensity (rel) )');
         ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         title(ttext)
      end
   end
   %% UME gamma
   plotnum=9;
   for ino=[16]
      disp('UME gamma')
      plotnum=plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=ba.(inst{ino,3});
      
      % find target time
      i = find(allt(6,:)>dnt(ii),1,'first');
      
      c = allc(:,i);
      varc = allvc(:,i);
      z = allz(:,i);
      t = nanmean(allt(:,i));
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(4,3,plotnum)
         plot_snippet
         xlabel('Gamma value');
         ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         title(ttext)
      end
   end
   
   %% UME chlorophyll
   plotnum=10;
   for ino=[8]
      disp('UME chlorophyll')
      plotnum=plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=ba.(inst{ino,3});
      
      if(ino==8 || ino==16) % no variance for Chl or gamma
         allvc=ones(size(allc));
      else
         allvc=vv.(inst{ino,3});
      end
      
      % find target time
      i = find(allt(6,:)>dnt(ii),1,'first');
      
      c = allc(:,i);
      varc = allvc(:,i);
      z = allz(:,i);
      t = nanmean(allt(:,i));
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(4,3,plotnum)
         plot_snippet
         xlabel('Chlorophyll value');
         ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         title(ttext)
      end
   end
   %% LISST density - special case
   plotnum=12;
   ino=[1]
   disp('LISST density')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{1,3})./ba.(inst{5,3});
   
   % find target time
   i = find(allt(6,:)>dnt(ii),1,'first');
   
   c = allc(:,i);
   allv = sqrt(vv.(inst{1,3}).^2+vv.(inst{5,3}).^2);
   varc = allv(:,i);
   
   z = allz(:,i);
   t = nanmean(allt(:,i));
   tmin = nanmin(nanmin(allt(:,i)));
   tmax = nanmax(nanmax(allt(:,i)));
   ok = (~isnan(c+z));
   disp([z c varc])
   
   if(sum(ok)>3)
      datestr(t)
      subplot(4,3,plotnum)
      plot_snippet
      xlabel('LISST Density Proxy');
      ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
      title(ttext)
      
   end
   %%
   shg
   pfn = sprintf('./%s/p%d.png',plotdir,fix(ydayc*100))
   print(pfn,'-dpng')
end
