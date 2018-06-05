% p_fit_mo_times - Fit settling velocities for several quantities at
% specific times
ifnorm = 0; % divide all of the values by the global mean
Ca_norm = 0; % normalize profiles by Ca
lglg = 0; % make log-log plots
load ustar_av % made by plot_ustar.m
burst_length = 3480/(3600*24);
dus = [ustar_av.dn]+burst_length/2.;
us = [ustar_av.us];
ew = sign([ustar_av.u]);
ewus = ew.*us;

% 20-minute bin averages
load suspsed_ba20_20rstrim_crs_cen
% use bin centers when available
inst = {...
   'LISSTtc','LISSTzc','LISSTvconc','LISST Volume';...
   'LISSTtc','LISSTzc','LISSTfinesv','LISST Fines';...
   'LISSTtc','LISSTzc','LISSTmicrov','LISST Micro';...
   'LISSTtc','LISSTzc','LISSTmacrov','LISST Macro';...
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
%
ydt = 276.2:(1/3)/24:285; % list of target year days
dnt = datenum(2011,1,1)+ydt-1;
datestr(dnt);
%%
for ii=1% 1:length(dnt)
   figure(1); clf
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
   % get closest u* value
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
   plotnum=5;
   for ino=[1]
      disp('LISST volume concentration')
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
         pf = pfit( c(ok), z(ok),0, za);
         loglog(c(ok)./pf.Ca,z(ok),'ok')
         hold on
         zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
         Cest = (zest./za).^pf.p;
         loglog(Cest,zest,'--k','linewidth',2)
         pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
         % normalize the nl fit by the linear Ca
         Cest_nl = pfnl.Ca./pf.Ca*(zest./za).^pfnl.p;
         loglog(Cest_nl,zest,'--r','linewidth',2)
         
         %xlim( [.5 500] )
         ylim( [.1 2.5] )
         xlim( [.01, 1.5])
         ylabel('Elevation (mab)')
         xlabel('Vol. Concentration (m^3/m^3)');
         title(inst{ino,4});
         ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
            datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
         ht = text(.015,.3,ts);
         set(ht,'fontsize',9);
      end
   end
   
   %% ABSS profiles
   plotnum=6;
   for ino=[13:15]
      disp('ABSS')
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
      t = nanmean(allt(:,i));
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         %         datestr(tmin)
         %         datestr(tmax)
         subplot(4,3,plotnum)
         pf = pfit( c(ok), z(ok),0, za);
         loglog(c(ok)./pf.Ca,z(ok),'ok')
         hold on
         zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
         Cest = (zest./za).^pf.p;
         loglog(Cest,zest,'--k','linewidth',2)
         pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
         % normalize the nl fit by the linear Ca
         Cest_nl = pfnl.Ca./pf.Ca*(zest./za).^pfnl.p;
         loglog(Cest_nl,zest,'--r','linewidth',2)
         
         %xlim( [.5 500] )
         ylim( [.1 2.5] )
         xlim( [.01, 1.5])
         ylabel('Elevation (mab)')
         xlabel('Backscatter Intensity (rel) )');
         title(inst{ino,4});
         ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
            datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
         ht=text(.015,.3,ts);
         set(ht,'fontsize',9);
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
      
      if(ifnorm)
         allc = allc/(nanmean(allc(:)));
      end
      
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
         %         datestr(tmin)
         %         datestr(tmax)
         subplot(4,3,plotnum)
         pf = pfit( c(ok), z(ok),0, za);
         loglog(c(ok)./pf.Ca,z(ok),'ok')
         hold on
         zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
         Cest = (zest./za).^pf.p;
         loglog(Cest,zest,'--k','linewidth',2)
         pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
         % normalize the nl fit by the linear Ca
         Cest_nl = pfnl.Ca./pf.Ca*(zest./za).^pfnl.p;
         loglog(Cest_nl,zest,'--r','linewidth',2)
         
         %xlim( [.5 500] )
         ylim( [.1 2.5] )
         xlim( [.01, 5])
         ylabel('Elevation (mab)')
         xlabel('Gamma value');
         title(inst{ino,4});
         ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
            datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
         ht = text(.015,.3,ts);
         set(ht,'fontsize',9);
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
      
      if(ifnorm)
         allc = allc/(nanmean(allc(:)));
      end
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
         %         datestr(tmin)
         %         datestr(tmax)
         subplot(4,3,plotnum)
         pf = pfit( c(ok), z(ok),0, za);
         loglog(c(ok)./pf.Ca,z(ok),'ok')
         hold on
         zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
         Cest = (zest./za).^pf.p;
         loglog(Cest,zest,'--k','linewidth',2)
         pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
         % normalize the nl fit by the linear Ca
         Cest_nl = pfnl.Ca./pf.Ca*(zest./za).^pfnl.p;
         loglog(Cest_nl,zest,'--r','linewidth',2)
         
         %xlim( [.5 500] )
         ylim( [.1 2.5] )
         xlim( [.01, 5])
         ylabel('Elevation (mab)')
         xlabel('Chlorophyll value');
         title(inst{ino,4});
         ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
            datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
         ht=text(.015,.3,ts);
         set(ht,'fontsize',9);
      end
   end
   %% LISST density - special case
   plotnum=11;
   for ino=[1]
      disp('LISST density')
      plotnum = plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=ba.(inst{1,3})./ba.(inst{5,3});
      
      if(ifnorm)
         allc = allc/(nanmean(allc(:)));
      end
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
         %         datestr(tmin)
         %         datestr(tmax)
         subplot(4,3,plotnum)
         plot_snippet
         xlabel('LISST Density Proxy');
         title('LISST Attenuation / Volume Conc.');
         
      end
   end
   %%
   shg
   pfn = sprintf('./plots/p%d.png',ii)
   %print(pfn,'-dpng')
end
