% script snippet invoked by p_fit_times_mo.m
pf = pfit( c(ok), z(ok),0, za);
zest = logspace( log10(pf.za),log10(max(z(ok))), 20);
Cest = pf.Ca*(zest./za).^pf.p;
plot(Cest,zest,'--k','linewidth',2)
hold on
pfnl = pfit_nlp( c(ok), z(ok), varc(ok),0, za);
% normalize the nl fit by the linear Ca
Cest_nl = pfnl.Ca*(zest./za).^pfnl.p;
if(Ca_norm)
   Cest    = Cest./pf.Ca
   Cest_nl = Cest_nl./pfnl.Ca
end
plot(Cest_nl,zest,'--r','linewidth',2);
if(Ca_norm)
   hp=plot(c(ok)./pf.Ca,z(ok),'ok');
else
   hp=plot(c(ok),z(ok),'ok');
end
set(hp,'markerfacecolor',[.4 .4 .4],'markeredgecolor',[.2,.2,.2])
%xlim([0 1.5])
ylim([0 2.5])
if(lglg)
   ylim([.1 2.5] )
   xlim([.01, 1.5])
end
ylabel('Elevation (m)')
ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
   datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
% ht=text(1,1.5,ts);
% set(ht,'fontsize',9);
