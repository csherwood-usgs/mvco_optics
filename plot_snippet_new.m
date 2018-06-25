% script snippet invoked by p_fit_times_mo.m
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
plot(Cest_nl,zest,'--b','linewidth',2);
hold on
plot(Cest_nl2,zest2,'--r','linewidth',2);
hp=plot(c(ok),z(ok),'ok');
set(hp,'markerfacecolor',[0 0 1],'markeredgecolor',[0,0,1])
hp2=plot(d(ok),z(ok),'ok');
set(hp2,'markerfacecolor',[1 0 0],'markeredgecolor',[1,0,0])
%xlim([0 1.5])
ylim([0.1 2])
if(lglg)
   ylim([.1 2] )
   %xlim([.01, 1.5])
end
ylabel('Elevation [m]')
ts = sprintf('%s\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\nCa_n=%7.2f\np_n=% 5.2f\nr_n^2=%06.4f',...
   datestr(t,'dd-mmm-yyyy HHMM'),pf.N,pf.Ca,-pf.p,pf.r2,pfnl.Ca,-pfnl.p,pfnl.r2)
   %datestr(t,'dd-mmm-yyyy HHMM'),pf2.N,pf2.Ca,-pf2.p,pf2.r2,pfnl2.Ca,-pfnl2.p,pfnl2.r2)
% ht=text(1,1.5,ts);
% set(ht,'fontsize',9);
