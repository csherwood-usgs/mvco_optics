% ws_fit_mo - Estmate ws and C from profile data

ifnorm = 1; % divide all of the values by the global mean

avinterval = '20'

load suspsed_ba20_20rstrim_crs_cen
use_bin_centers=1 % use the bin centers when available
if(use_bin_centers)
inst = {...
   'LISSTtc','LISSTzc','LISSTvconc','LISST Volume Conc.';...
   'LISSTtc','LISSTzc','LISSTfinesv','LISST Fines Vol. Conc.';...
   'LISSTtc','LISSTzc','LISSTmicrov','LISST Microfloc Vol. Conc.';...
   'LISSTtc','LISSTzc','LISSTmacrov','LISST Macrofloc Vol. Conc.';...
   'LISSTtc','LISSTzc','LISSTattn','LISST Attenuation';...
   'YSIdn','YSIz','YSIturb','YSI Turbidity';... % no bin centers?
   'UMEtc','UMEzc','UMEattn650','c_{p}_{650}';...
   'UMEtc','UMEzc','UMEchl','Chl';...
   'UMEtc','UMEzc','UMEbs650','bs_{650}';...
   'trandn','tranz','tranattn','Trans. Attn.';...
   'ADVdn','OBSz','OBS','OBS Backscatter';...
   'ADVdn','ADVz','ADVagc','ADV AGC';...
   'absstc','absszc','abss1','ABSS 1';...
   'absstc','absszc','abss2','ABSS 2';...
   'absstc','absszc','abss3','ABSS 3';...
   };
else
inst = {...
   'LISSTdn','LISSTz','LISSTvconc';...
   'LISSTdn','LISSTz','LISSTfinesv';...
   'LISSTdn','LISSTz','LISSTmicrov';...
   'LISSTdn','LISSTz','LISSTmacrov';...
   'LISSTdn','LISSTz','LISSTattn';...
   'YSIdn','YSIz','YSIturb';...
   'UMEdn','UMEz','UMEattn650';...
   'UMEdn','UMEz','UMEchl';...
   'UMEdn','UMEz','UMEbs650';...
   'trandn','tranz','tranattn';...
   'ADVdn','OBSz','OBS';...
   'ADVdn','ADVz','ADVagc';...
   'abssdn','abssz','abss1';...
   'abssdn','abssz','abss2';...
   'abssdn','abssz','abss3';...
   };
end

za = 0.1; % (m) Ca estimates standardized to this elevation
nav = 0;  % number of profiles to left and right to be averaged
%%
ino = 2;
for ino=[13:15]
   allt=ba.(inst{ino,1});
   allz=ba.(inst{ino,2});
   allc=ba.(inst{ino,3});

   if(ifnorm)
      allc = allc/(nanmean(allc(:)));
   end
   if(ino==5) % convert LISST to mass (roughly)
      % allc = (allc*(2.7/1.5)-1.3); % g/m3
      allc = (allc*(2.7/1.5)); % g/m3
   end
   if(ino==7) % convert attn to Mass (forget offset...can lead to neg. Concs)
      % allc = (allc*2.7-1.3); %g/m3
      allc = (allc*2.7); %g/m3
   end
   
   if(ino==8) % no variance for Chl (why not?)
      allvc=ones(size(allc))
   else
      allvc=vv.(inst{ino,3});
   end
   dz = diff(allz(1:2,1));
   [ nz nt ]= size( allz )
   
   ic = 1;
   for i=1:nt
      c = allc(:,i);
      varc = allvc(:,i);
      z = allz(:,i);
      t = nanmean(allt(:,i))
      tmin = nanmin(nanmin(allt(:,i)));
      tmax = nanmax(nanmax(allt(:,i)));
      ok = (~isnan(c+z))
      
      if(sum(ok)>3)
         datestr(t)
         datestr(tmin)
         datestr(tmax)
         figure(1); clf
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
         xlim( [.01, 1.2])
         ylabel('Elevation (mab)')
         xlabel('Attenuation ( m^{-1} )');
         title(inst{ino,3});
         ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
            datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
         text(.015,.3,ts)
         shg
         pause(.2)
         pf.i = i;
         pf.t = t;
         pfa(ic)=pf
         pfanl(ic)=pfnl;
         ic = ic+1;
      end
   end
   pfaname = [inst{ino,3},avinterval,'_pfa_mo']
   pfanamenl = [inst{ino,3},avinterval,'_pfnl_mo']
   eval([pfaname,'= pfa;'])
   eval([pfanamenl,'= pfanl;'])
   eval(['save ',pfanamenl,' ',pfaname,' ',pfanamenl])
end