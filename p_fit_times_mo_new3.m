close all
clear all
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
us_wave = [ustar_av.us];
us_wc=[ustar_av.ustrr];
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
%compute the Sauter diameter

for i=1:12
   for j=1:2442
      Lisst_sumV(i,j)=nansum(ba.LISSTvconc(i,j,:));
      Lisst_sumV_fines(i,j)=nansum(ba.LISSTvconc(i,j,4:12));
      Lisst_sumV_micro(i,j)=nansum(ba.LISSTvconc(i,j,13:23));
      Lisst_sumV_macro(i,j)=nansum(ba.LISSTvconc(i,j,24:29));
   end
end
D=ba.LISSTsz;
for i=1:12
   for j=1:2442
      for z=4:29
         Lisst_A(i,j,z)=ba.LISSTvconc(i,j,z)*1.5/D(z);
      end
      Lisst_sumA(i,j)=nansum(Lisst_A(i,j,:));
   end
end
D_S=3/2*Lisst_sumV./Lisst_sumA;
Agg_Dens=ba.(inst{5,3})./ba.(inst{1,3});
f_over_m=Lisst_sumV_fines./Lisst_sumV_micro;
f_over_ma=Lisst_sumV_fines./Lisst_sumV_macro;
za = 0.1; % (m) Ca estimates standardized to this elevation
nav = 0;  % number of profiles to left and right to be averaged
%list of target year days
% ydt = [261:(1/3)/24:263];
% plotdir = 'maria'
%
%ydt = [261:(1/3)/24:261.3];
%ydt = [268+(18/3)/24:(1/3)/24:269];
%ydt = [275:(1/3)/24:276];
ydt =[281+(18/3)/24:(1/3)/24:282];
%,268+(2/3)/24:(1/3)/24:269,275:(1/3)/24:276];%,281+(2/3)/24:(1/3)/24:282,286:(1/3)/24:287];
% plotdir = 'spgtides'

% ydt = [274:(1/3)/24:276];
% plotdir = 'ophelia'
%
%ydt = [280:(1/3)/24:282];
plotdir = 'new';
%
% ydt = [286:(1/3)/24:288];
% plotdir = 'noreaster'

yday_off = datenum('1-Jan-2011 00:00:00');
dnt = yday_off+ydt;
datestr(dnt);
%%
ii=1
N=3;

for ii=1:length(dnt)/N
   
   figure(1); clf
   %% info box
   %subplot(4,3,3)
   %    xlim([0 1])
   %    ylim([0 1])
   speed = interp1(dus,ewspeed,dnt(1+(ii-1)*N:ii*N),'nearest'); speed=nanmean(speed);
   ustar = interp1(dus,us,dnt(1+(ii-1)*N:ii*N),'nearest'); ustar=nanmean(ustar);
   ustar_wc = interp1(dus,us_wc,dnt(1+(ii-1)*N:ii*N),'nearest'); ustar_wc=nanmean(ustar_wc);
   ubr = interp1(dus,ustar_av.ubr,dnt(1+(ii-1)*N:ii*N),'nearest'); ubr=nanmean(ubr);
   zo = interp1(dus,ustar_av.zoa,dnt(1+(ii-1)*N:ii*N),'nearest'); zo=nanmean(zo);
   J=floor(nanmedian(1+(ii-1)*N:ii*N));
   text(0.1,.9,datestr(dnt(J)))
   ydayc = dnt(J)-datenum('1-Jan-2011 00:00:00');
   text(0.1,.8,sprintf('Day: %6.2f',ydayc))
   text(0.1,.7,sprintf('Speed: %5.2f m/s',speed))
   text(0.1,.6,sprintf('ustar_wc: % 4.3f m/s',ustar_wc))
   text(0.1,.5,sprintf('u*c: % 4.3f m/s',ustar))
   text(0.1,.4,sprintf('zo : % 6.4f m',zo))
   set(gca,'Visible','off')
   %pause
   
   figure(2); clf
   input('Resize Figure 2 and hit return')
   %% Cp650
   plotnum=1;
   ino=[7]
   disp('UME Cp650')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});
   alld=ba.UMEbs650;
   allvc=vv.(inst{ino,3});
   allvd=vv.UMEbs650;
   %    if(ifnorm)
   %       allc = allc/(nanmean(allc(:)));
   %    end
   
   % find target time
   i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));

   % cp 650 attenuation
   c = nanmean(allc(:,i)')';
   varc = allvc(:,i);
   
   % bs 650 backscatter
   % scale to fit on same x axis
%    d = nanmean(alld(:,i)')'*40;
%    vard = allvd(:,i)*40;
   % no scaling for dual axes
   d = nanmean(alld(:,i)')';
   vard = allvd(:,i);
   
   z = allz(:,floor(mean(i)));
   t = nanmean(allt(:,floor(mean(i))));
   fprintf(1,'Target time: %s; mean profile time: %s\n',datestr(dnt(ii)),datestr(t));
   %ustar = interp1(dus,us,t,'nearest')
   
   ok = (~isnan(c+z));
   
   % colors for cp and bpp
   cpcol = [.2 .2 1]
   bpcol = [1 .2 .2]
   
   if(sum(ok)>3)
      subplot(2,3,1)
      % script plot_snippet is replaced by custom text here
      pf = pfit( c(ok), z(ok),0, za);
      zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
      pf2 = pfit( d(ok), z(ok),0, za);
      zest2 = logspace( log10(pf2.za),log10(max(z(ok))), 20);
      %Cest = pf.Ca*(zest./za).^pf.p;
      %plot(Cest,zest,'--k','linewidth',2)
      %hold on
      pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
      pfnl2 = pfit_nlp( d(ok), z(ok), vard(ok),0, za);
      % normalize the nl fit by the linear Ca
      Cest_nl = pfnl.Ca*(zest./za).^pfnl.p;
      Cest_nl2 = pfnl2.Ca*(zest2./za).^pfnl2.p;
      % plot 
      plot(Cest_nl,zest,'--b','linewidth',2);
         hold on
      hp=plot(c(ok),z(ok),'ob');
      set(hp,'MarkerFaceColor',[0 0 1])
      set(gca,'YLim',[.1 2])
      set(gca,'XColor',[0 0 1])

      ylabel('Elevation [m]')
      ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
         datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)

      hx=xlabel('c_p(650) and b_{bp}(650) m^{-1} )');
      set(hx,'Color','k')
      
      % get axis info
      ax1 = gca;
      ax1_pos = ax1.Position;
      ax1_YLim = ax1.YLim;
            
      % create second axis on top
      ax2 = axes('Position',ax1_pos,...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none')
      line(Cest_nl2,zest2,'Parent',ax2,'Color','r')
      hold on
      hp2=plot(d(ok),z(ok),'or');
      set(hp2,'MarkerFaceColor',[1 0 0])

      set(gca,'YLim',ax1_YLim);
      set(gca,'XColor',[1 0 0 ])
      set(gca,'YTickLabel',[])
      
      ttext1 = sprintf('ws_{cp}= %4.2f  ws_{bbp}=%4.2f',-1000*0.41*pfnl.p*ustar,-1000*0.41*pfnl2.p*ustar)
      ht=title(ttext1);

   end
   
   nanmean(c)
   nanstd(c)
   nanmean(d/40)
   nanstd(d/40)
   
   
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
   i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));
   c = nanmean(allc(:,i)')';
   d = nanmean(alld(:,i)')';
   e = nanmean(alle(:,i)')';
   gamma_bb=(650/532)*log(d./e);
   z = allz(:,floor(mean(i)));
   t = nanmean(allt(:,floor(mean(i))));
   tmin = nanmin(nanmin(allt(:,i)));
   tmax = nanmax(nanmax(allt(:,i)));
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,2)
      %plot_snippet
      hp=plot(c(ok),z(ok),'ok');
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      hp2=plot(gamma_bb(ok)+0.2,z(ok),'ok');
      set(hp2,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0])
      xlabel('\gamma_{cp}(B), \gamma_{bbp}+0.2(R)');
      ylabel('Elevation [m]')
      hold off
      ylim([0.1 2])
      %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
      %title(ttext)
   end
   nanmean(c(ok))
   nanstd(c(ok))
   nanmean(gamma_bb(ok))
   nanstd(gamma_bb(ok))
   
   
   
   %% LISST D50
   plotnum=3;
   ino=[1]
   disp('LISST D50')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=1./Agg_Dens;
   %allc=ba.LISSTD50a;
   alld=D_S;
   
   %    if(ifnorm)
   %       allc = allc/(nanmean(allc(:)));
   %    end
   
   % find target time
   i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));
   c = nanmean(allc(:,i)')';
   d = nanmean(alld(:,i)')';
   z = allz(:,floor(mean(i)));
   t = nanmean(allt(:,floor(mean(i))));
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      %         datestr(tmin)
      %         datestr(tmax)
      subplot(2,3,3)
      %     plot_snippet
      hp=plot(c(ok)*8,z(ok),'ok');
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      hp2=plot(d(ok),z(ok),'ok');
      set(hp2,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0])
      hold off
      
      ylabel('Elevation [m]')
      xlabel('8\times\rho_a^{-1}(B), D_s(R)');
      %       ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
      %       title(ttext)
   end
   ylim([0.1 2])
   nanmean(c(ok))
   nanstd(c(ok))
   nanmean(d(ok))
   nanstd(d(ok))
   
   
   %% Backscattering ratio
   plotnum=4;
   for ino=[8]
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
      i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));
      c = nanmean(allc(:,i)')';
      z = allz(:,floor(mean(i)));
      t = nanmean(allt(:,floor(mean(i))));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(2,3,4)
         hp=plot(c(ok),z(ok),'ok');
         set(hp,'markerfacecolor','k','markeredgecolor','k')
         ylabel('Elevation [m]')
         xlabel('b_{bp}/b_p(650)');
         %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         %title(ttext)
         ylim([0.1 2])
      end
   end
   
   nanmean(c(ok))
   nanstd(c(ok))
   
   
   
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
      i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));
      c = nanmean(allc(:,i)')';
      z = allz(:,floor(mean(i)));
      t = nanmean(allt(:,floor(mean(i))));
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(2,3,5)
         hp=plot(c(ok),z(ok),'ok');
         set(hp,'markerfacecolor','k','markeredgecolor','k')
         ylabel('Elevation [m]')
         xlabel('Chl/c_p(650)');
         %ttext = sprintf('%s w_s: %4.2f r^2: %4.2f',inst{ino,4},-1000*0.41*pfnl.p*ustar,pfnl.r2);
         %title(ttext)
         ylim([0.1 2])
      end
   end
   
   nanmean(c(ok))
   nanstd(c(ok))
   
   
   %% LISST density - special case
   plotnum=6;
   ino=[1]
   %disp('LISST density')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=f_over_m;
   alld=f_over_ma*3;
   %ba.(inst{5,3})./ba.(inst{1,3});
   
   % find target time
   i = find(allt(6,:)>=dnt(1+(ii-1)*N) & allt(6,:)<dnt(ii*N));
   c = nanmean(allc(:,i)')';
   d = nanmean(alld(:,i)')';
   z = allz(:,floor(mean(i)));
   t = nanmean(allt(:,floor(mean(i))));
   
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,6)
      %plot_snippet
      hp=plot(c(ok),z(ok),'ok');
      ylabel('Elevation [m]')
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      hp2=plot(d(ok),z(ok),'ok');
      set(hp2,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0])
      hold off
      ylim([0.1 2])
      xlabel('V_f/V_m (B) 3xV_f/V_M(R)');
      nanmean(c(ok))
      nanstd(c(ok))
      nanmean(d(ok)/3)
      nanstd(d(ok)/3)
      
      
      
   end
   %%
   %pause
   shg
   pfn = sprintf('./%s/p%d.png',plotdir,fix(ydayc*100))
   print(pfn,'-dpng')
end
