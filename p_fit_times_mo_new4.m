% p_fit_times_mo_new4.m
% Script to plot profiles. Version 3 was made by EB and the plotting was
% modified by CRS. This version is changed as follow:
% Parameters to be plotted are extracted from the database with names at
% the top. This is a waste of memory, but makes the plotting code easier to
% read and understand.
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

%compute the Sauter diameter and some other stuff we will plot
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

%% Here is where we set up cases to plot


% interpolate wave stresses to profile times
ti = ba.ADVtc(5,:);
uwi = interp1(dus,us_wave,ti);

% indices into wave conditions
big_wave = uwi>.01;
small_wave = uwi<=0.01;

% indices into current conditions
flood = ba.u(5,:)>5;
ebb = ba.u(5,:)<-5;
slack = ba.u(5,:)>=-5 & ba.u(5,:)<=5;

% Named periods. Adjust these, or add to them
cases(1).name = 'Maria';
cases(1).ts = 261.;
cases(1).te = 262.;
cases(2).name = 'SpgTide';
cases(2).ts = 267.1;
cases(2).te = 268.1;
cases(3).name = 'Ophelia';
cases(3).ts = 275.1;
cases(3).te = 276.1;
cases(4).name = 'Calm';
cases(4).ts = 280;
cases(4).te = 282;
cases(5).name = 'Noreaster';
cases(5).ts = 286;
cases(5).te = 287;

if(1) % replace the times with my best guess of Emmanuel's times
   % little function to convert hour to decimal day and subtract two
   % minutes
   dch =@(h) h/24.- 2./60./24.;
   % 62 minutes later (three profiles)
   min62 = 1./24.+(2./60.)/24.
   cases(1).name = 'Fig3';
   cases(1).ts = 261.+dch(2.)
   cases(1).te = cases(1).ts+min62;
   cases(2).name = 'Fig4';
   cases(2).ts = 268.+dch(14.);
   cases(2).te = cases(2).ts+2*min62;
   cases(3).name = 'Fig5';
   cases(3).ts = 275.+dch(4.)
   cases(3).te = cases(3).ts+2*min62;
   cases(4).name = 'Fig6';
   cases(4).ts = 281.+dch(22.)
   cases(4).te = cases(4).ts+min62;
   cases(5).name = 'Fig6';
   cases(5).ts = 286.+dch(2.)
   cases(5).te = cases(5).ts+min62;
end
% conversion from datenum time to year day
yday_off = datenum('1-Jan-2011 00:00:00');
tiyd = ti-yday_off;

plotdir = 'test'
% open a file to record stats
sfname = ['./',plotdir,'/stats.txt']
fid = fopen(sfname,'w');

% fix that xlimits for all plots?
fixx = 1
for icase = 1:length(cases)
   ba_tindex = tiyd >= cases(icase).ts & tiyd <= cases(icase).te;
   us_tindex = dus-yday_off >= cases(icase).ts & dus-yday_off <= cases(icase).te;
   disp(cases(icase).name)
   
   % indices
   ius = us_tindex;
   i = ba_tindex;
   disp(sum(i))
   
   % info box
   figure(1); clf
   sp_t=nanmean(abs(speed(ius)));
   us_t=nanmean(us(ius));
   uswc_t=nanmean(us_wc(ius));
   ubr_t=nanmean(ubr(ius));
   zo_t=nanmean(zoa(ius));
   J=nanmean(ti(i));
   
   text(0.1,.9,datestr((J)))
   ydayc = nanmean(tiyd(i));
   text(0.1,.8,sprintf('Days: %6.2f - %6.2f',min(tiyd(i)),max(tiyd(i))))
   text(0.1,.7,sprintf('Speed: %5.2f m/s',sp_t))
   text(0.1,.6,sprintf('u*wc: % 4.3f m/s',uswc_t))
   text(0.1,.5,sprintf('u*c: % 4.3f m/s',us_t))
   text(0.1,.4,sprintf('zo : % 6.4f m',zo_t))
   set(gca,'Visible','off')
   %pause
   fprintf(fid,'%5d, %6.2f, %6.2f, %5.2f, %5.3f, %5.3f, %6.4f\n',fix(100*ydayc),min(tiyd(i)),max(tiyd(i)),...
      sp_t,uswc_t,us_t,zo_t);
   
   figure(2); clf
   fig = gcf;
   fig.Position =[397 288 913 1027];
   fig.PaperPositionMode = 'auto';
   
   %% Cp650 and Bbp650
   plotnum=1;
   ino=[7]
   disp('UME Cp650')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});
   alld=ba.UMEbs650;
   allvc=vv.(inst{ino,3});
   allvd=vv.UMEbs650;
   
   % cp 650 attenuation
   c = nanmean(allc(:,i)')';
   varc = nanvar(allvc(:,i)')';
   sdc = nanstd(allc(:,i)')';
   
   d = nanmean(alld(:,i)')';
   vard =nanvar(allvd(:,i)')';
   sdd = nanstd(alld(:,i)')';
   
   z = nanmean(allz(:,i)')';
   t = nanmean(allt(:,i)')';
   fprintf(1,'Mean profile time: %s\n',datestr(t));
   
   ok = (~isnan(c+z));
   
   % colors for cp and bpp
   cpcol = [.2 .2 1]
   bpcol = [1 .2 .2]
   %
   if(sum(ok)>3)
      % script plot_snippet is replaced by custom text here
      pf = pfit( c(ok), z(ok),0, za);
      zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
      pf2 = pfit( d(ok), z(ok),0, za);
      zest2 = logspace( log10(pf2.za),log10(max(z(ok))), 20);
      
      pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
      pfnl2 = pfit_nlp( d(ok), z(ok), vard(ok),0, za);
      % normalize the nl fit by the linear Ca
      Cest_nl = pfnl.Ca*(zest./za).^pfnl.p; % cp
      Cest_nl2 = pfnl2.Ca*(zest2./za).^pfnl2.p; % bs
      
      
      % plot cp in blue
      subplot(2,3,1)
      plot(Cest_nl,zest,'--b','linewidth',2);
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+b')
      end
      hp=plot(c(ok),z(ok),'ob');
      set(hp,'MarkerFaceColor',[0 0 1])
      set(gca,'YLim',[.1 2])
      set(gca,'XColor',[0 0 1])
      % if(fixx), set(gca,'XLim',[0 40]); end
      set(gca,'Xscale','log')
      set(gca,'XLim',[.1, 50])
      set(gca,'Xtick',[.1, 1, 10])
      ylabel('Elevation [m]','fontsize',14)
      set(gca,'fontsize',11)
      text(0.05, 0.97,'a','units','normalized','fontsize',14)
      ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
         datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
      
      hx=xlabel('{\color{blue}{\itc_{p(650)}}} and {\color{red}{\itb_{bp(650)}}} [m^{-1}]',...
         'fontsize',14);
      set(hx,'Color','k')
      
      % get axis info
      ax1 = gca;
      ax1_pos = ax1.Position;
      ax1_YLim = ax1.YLim;
      
      % create second axis on top
      ax2 = axes(...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none')
      
      % plot bp in red
      hl=line(Cest_nl2,zest2,'Parent',ax2,'Color','r')
      set(hl,'Linewidth',2)
      hold on
      for ik = 1:length(ok)
         plot([max(0., d(ik)-sdd(ik)), d(ik)+sdd(ik)],[z(ik) z(ik)],'+r')
      end
      hp2=plot(d(ok),z(ok),'or');
      set(hp2,'MarkerFaceColor',[1 0 0])
      
      set(gca,'YLim',ax1_YLim);
      set(gca,'XColor',[1 0 0 ])
      %if(fixx), set(gca,'XLim',[0 0.7]); end
      set(gca,'Xscale','log')
      set(gca,'XLim',[.01, 1])
      Yax=get(gca,'Yaxis')
      set(Yax,'Visible','off')
      set(gca,'Position',ax1_pos)
      set(ax1,'Position',ax1_pos)
      set(gca,'fontsize',11)
      
      ttext1 = sprintf('ws_{cp}=%4.2f ws_{bbp}=%4.2f',-1000*0.41*pfnl.p*us_t,-1000*0.41*pfnl2.p*us_t)
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
   allvd=vv.UMEbs532;
   alle=ba.UMEbs650;
   allve=vv.UMEbs650;
   
   c = nanmean(allc(:,i)')';
   sdc = nanstd(allc(:,i)')';
   
   d = nanmean(alld(:,i)')';
   e = nanmean(alle(:,i)')';
   
   gamma_bb=(650/532)*log(d./e);
   sd_gbb = nanstd( (650/532)*log(alld(:,i)./alle(:,i))')';
   
   z = nanmean(allz(:,i)')';
   t = nanmean(allt(:,i)')';
   tmin = nanmin(nanmin(allt(:,i)));
   tmax = nanmax(nanmax(allt(:,i)));
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,2)
      
      % plot gamma_cp in blue
      hp=plot(c(ok),z(ok),'ob');
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+b')
      end
      set(gca,'YLim',[.1 2])
      if(fixx), set(gca,'XLim',[0 0.6]); end
      text(0.05, 0.97,'b','units','normalized','fontsize',14)
      
      set(gca,'XColor',[0 0 1])
      set(gca,'fontsize',11)
      set(gca,'Yticklabel',[])
      
      %ylabel('Elevation [m]')
      ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
         datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
      hx=xlabel('{\it{\color{blue}\gamma_{cp}}} and {\it{\color{red}\gamma_{bbp}}} [ ]',...
         'fontsize',14);
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
      
      % plot gamma bb in red
      % (for some reason, the plot command won't work...so use line, but
      % make it invisble below. WTF)
      hp2=line(gamma_bb(ok),z(ok),'Parent',ax2,'Color','r');
      hold on
      hp3=plot(gamma_bb(ok),z(ok),'or')
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,gamma_bb(ik)-sd_gbb(ik)), gamma_bb(ik)+sd_gbb(ik)],[z(ik) z(ik)],'+r')
      end
      set(hp3,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0],'Parent',ax2)
      set(hp2,'Color','none')
      
      set(gca,'YLim',ax1_YLim);
      set(gca,'XColor',[1 0 0 ])
      set(gca,'fontsize',11)
      if(fixx), set(gca,'XLim',[0 0.4]); end
      Yax=get(gca,'Yaxis')
      set(Yax,'Visible','off')
      set(gca,'Position',ax1_pos)
      set(ax1,'Position',ax1_pos)
      
      
   end
   nanmean(c(ok))
   nanstd(c(ok))
   nanmean(gamma_bb(ok))
   nanstd(gamma_bb(ok))
   
   %% LISST D50 and rho-1
   plotnum=3;
   ino=[1]
   disp('LISST D50')
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=1./Agg_Dens;
   %allc=ba.LISSTD50a;
   alld=D_S;
   
   % find target time
   c = nanmean(allc(:,i)')';
   sdc=  nanstd(allc(:,i)')';
   d = nanmean(alld(:,i)')';
   sdd = nanstd(alld(:,i)')';
   z = nanmean(allz(:,i)')';
   t = nanmean(allt(:,i)')';
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,3)
      
      % plot rho-1 in blue
      hp=plot(c(ok),z(ok),'ob');
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+b')
      end
      set(gca,'YLim',[.1 2])
      set(gca,'XColor',[0 0 1])
      set(gca,'YTickLabel',[])
      if(fixx), set(gca,'XLim',[0 30]); end
      text(0.05, 0.97,'c','units','normalized','fontsize',14)
      
      set(gca,'fontsize',11)
      
      %ylabel('Elevation [m]')
      hx=xlabel('{\color{blue}{\it\rho_a^{-1}} [\mum^3 m^{-1}]} and {\color{red}{\itD_s}[\mum]}',...
         'fontsize',14);
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
      
      % plot Ds in red
      % (for some reason, the plot command won't work...so use line, but
      % make it invisble below. WTF)
      hp2=line(d(ok),z(ok),'Parent',ax2,'Color','r');
      hold on
      for ik = 1:length(ok)
         plot([max(0., d(ik)-sdd(ik)), d(ik)+sdd(ik)],[z(ik) z(ik)],'+r')
      end
      hp3=plot(d(ok),z(ok),'or')
      set(hp3,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0],'Parent',ax2)
      set(hp2,'Color','none')
      set(gca,'YLim',ax1_YLim);
      if(fixx), set(gca,'XLim',[0 250.]); end
      
      set(gca,'XColor',[1 0 0 ])
      set(gca,'fontsize',11)
      Yax=get(gca,'Yaxis')
      set(Yax,'Visible','off')
      set(gca,'Position',ax1_pos)
      set(ax1,'Position',ax1_pos)
   end
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
      
      % find target time
      c = nanmean(allc(:,i)')';
      sdc = nanstd(allc(:,i)')';
      
      z = nanmean(allz(:,i)')';
      t = nanmean(allt(:,i)')';
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(2,3,4)
         hp=plot(c(ok),z(ok),'ok');
         hold on
         set(hp,'markerfacecolor','k','markeredgecolor','k')
         % error bars
         for ik = 1:length(ok)
            plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+k')
         end
         ylabel('Elevation [m]','fontsize',14)
         text(0.05, 0.97,'d','units','normalized','fontsize',14)
         
         xlabel('{\itb_{bp(650)}/b_{p(650)}} [ ]','fontsize',14);
         ylim([0.1 2])
         if(fixx), set(gca,'XLim',[0 0.05]); end
         
         set(gca,'fontsize',11)
      end
   end
   
   
   
   %% UME chlorophyll over attenuation
   plotnum=5;
   for ino=[8]
      disp('UME chlorophyll over attenuation')
      plotnum=plotnum+1;
      allt=ba.(inst{ino,1});
      allz=ba.(inst{ino,2});
      allc=(ba.(inst{ino,3})./(ba.(inst{7,3})));
      
      c = nanmean(allc(:,i)')';
      sdc = nanstd(allc(:,i)')';
      z = nanmean(allz(:,i)')';
      t = nanmean(allt(:,i)')';
      ok = (~isnan(c+z));
      
      if(sum(ok)>3)
         datestr(t)
         subplot(2,3,5)
         hp=plot(c(ok),z(ok),'ok');
         hold on
         % error bars
         for ik = 1:length(ok)
            plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+k')
         end
         set(hp,'markerfacecolor','k','markeredgecolor','k')
         xlabel('{\itChl/c_{p(650)}} [mg m^{-2}]','fontsize',14);
         if(fixx), set(gca,'XLim',[0 5.5]); end
         text(0.05, 0.97,'e','units','normalized','fontsize',14)
         
         ylim([0.1 2])
         set(gca,'fontsize',11)
         set(gca,'Yticklabel',[])
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
   alld=f_over_ma;
   
   c = nanmean(allc(:,i)')';
   sdc = nanstd(allc(:,i)')';
   
   d = nanmean(alld(:,i)')';
   sdd = nanstd(alld(:,i)')';
   
   z = nanmean(allz(:,i)')';
   t = nanmean(allt(:,i)')';
   ok = (~isnan(c+z));
   
   if(sum(ok)>3)
      datestr(t)
      subplot(2,3,6)
      
      % plot rho-1 in blue
      hp=plot(c(ok),z(ok),'ob');
      set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,c(ik)-sdc(ik)), c(ik)+sdc(ik)],[z(ik) z(ik)],'+b')
      end
      set(gca,'YLim',[.1 2])
      set(gca,'XColor',[0 0 1])
      set(gca,'Yticklabels',[])
      if(fixx), set(gca,'XLim',[0 0.4]); end
      text(0.05, 0.97,'f','units','normalized','fontsize',14)
      
      set(gca,'fontsize',11)
      hx=xlabel('{\color{blue}{\itV_f / V_m}} and {\color{red}{\itV_f / V_M}} [ ]',...
         'fontsize',14);
      
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
      
      % plot Ds in red
      % (for some reason, the plot command won't work...so use line, but
      % make it invisble below. WTF)
      hp2=line(c(ok),z(ok),'Parent',ax2,'Color','b');
      hold on
      
      
      hp3=plot(d(ok),z(ok),'or')
      hold on
      % error bars
      for ik = 1:length(ok)
         plot([max(0,d(ik)-sdd(ik)), d(ik)+sdd(ik)],[z(ik) z(ik)],'+r')
      end
      set(hp3,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0],'Parent',ax2)
      set(hp2,'Color','none')
      set(gca,'YLim',ax1_YLim);
      if(fixx), set(gca,'XLim',[0 0.4]); end
      
      set(gca,'XColor',[1 0 0 ])
      set(gca,'fontsize',11)
      Yax=get(gca,'Yaxis')
      set(Yax,'Visible','off')
      set(gca,'Position',ax1_pos)
      set(ax1,'Position',ax1_pos)
      
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
fclose(fid);
