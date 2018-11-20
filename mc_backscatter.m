% mc_backscatter
%
% model backscattered area with occlusion

np = 1e3;
xp = 1:np;
yp = 1:np;
Vv = np^3
iplot=1;
nit = 100;

rlist = [12, 20, 50, 100];
nr = length(rlist)

Vc = logspace(log10(.001),log10(.1),20);
%Vc = .001

nc = length(Vc);
N = NaN*ones(nc,nr);
Nr = NaN*ones(nc,nr);
Ap = NaN*ones(nc,nr);
Asm_mean = zeros(nc,nr);
Asm_std = zeros(nc,nr);
Apm_mean = zeros(nc,nr);
Apm_std = zeros(nc,nr);

% loop through the range of particle radii
for ir=1:nr
   r = rlist(ir);
   fprintf(1,'r = %4.0f\n',r)
   fprintf(1,'Vc,  N, Ap, amean, astd, amean/Ap\n')
   
   % loop through the range of volume concentrations
   for m=1:nc
      Vp = (4./3.)*pi*r^3;
      Nr(m,ir) = Vc(m)*Vv/Vp;
      N(m,ir) = max(1,round(Vc(m)*Vv/Vp));
      Ap(m,ir) = Nr(m,ir)*pi*r^2;
      Asi = NaN*ones(nit,1);
      Api = NaN*ones(nit,1);
      Asi_mean = NaN*ones(nit,1);
      Api_std = NaN*ones(nit,1);
      
      % do a bunch of simulations for each case
      for it=1:nit
         % pick random locations for the particle centers
         xolist = np*rand(N(m,ir),1);
         yolist = np*rand(N(m,ir),1);
         
         % zero out pixel arrays to sum particle areas
         zp = zeros(round(np),round(np));
         zps = zeros(round(np),round(np));
         
         % loop through all of the particles
         for k = 1:N(m,ir)
            xo = xolist(k);
            yo = yolist(k);
            for i = -r:r
               for j=-r:r
                  % increment pixel arrays overlapped by each particles 
                  if( ...
                        ( (i^2 <= r^2-j^2 ) && ...
                        (j^2 <= r^2-i^2 ) ) && ...
                        ( round(i+yo) > 0 && round(i+yo) <= np ) && ...
                        ( round(j+xo) > 0 && round(j+xo) <= np ) )
                     % this array is set to one if any particle overlaps it
                     zp(round(i+yo), round(j+xo) )=1.;
                     % this array is increments for each particle that
                     % overlaps it...so it counts all particle areas
                     zps(round(i+yo), round(j+xo) ) = 1.+zps (round(i+yo), round(j+xo) );
                  end
               end
            end
         end
         Asi(it) = sum( zp(:) );
         Api(it) = sum( zps(:) );
         Asi_mean(it) = nanmean(Asi);
         Asi_std(it) = nanstd(Asi);
         
         % plot the first four simulations
         if(it==1), figure(1); clf, end
         if(it<=4)
            subplot(2,2,it)
            pcolor(zps)
            shading interp
            axis equal
            axis([0 1000 0 1000])
            caxis([0 4])
            set(gca,'xticklabels',[],'yticklabels',[])
            ts = sprintf('V_c=%6.4f, r=%4.0f, N=%d, A_s/A_p=%5.2f',Vc(m), r, N(m,ir), Asi(it)/Api(it));
            title(ts,'fontsize',10)
         end
         if(it==5)
            figure(1)
            colormap hot
            shg
            pfname = sprintf('fVc%03d_r%02d.png',round(1000*Vc(m)),round(r));
            print(pfname, '-dpng', '-r200')
         end
         
      end
      Asm_mean(m,ir) = mean(Asi);
      Apm_mean(m,ir) = mean(Api);
      Asm_std(m,ir) = std(Asi);
      Apm_std(m,ir) = std(Api);
      
      % plot to test for convergence of iterations
      if(1)
         figure(2); clf
         subplot(211)
         plot([1:100],ones(nit,1),'--k')
         hold on
         plot(Asi_mean./Asm_mean(m,ir),'linewidth',2)
         ylabel('Mean As(1:i) / Mean As')
         ts = sprintf('V_c = %6.4f, r = %4.0f, N = %d, A_s/A_p = %5.2f',Vc(m), r, N(m,ir), Asm_mean(m,ir)/Apm_mean(m,ir));
         title(ts)
         set(gca,'fontsize',14)
         
         subplot(212)
         plot([1:100],ones(nit,1),'--k')
         hold on
         plot(Asi_std./Asm_std(m,ir),'linewidth',2)
         ylabel('Std. As(1:i) / Std(As)')
         xlabel('Iteration')
         set(gca,'fontsize',14)
         shg
         pfname = sprintf('convergeVc%03d_r%02d.png',round(1000*Vc(m)),round(r));
         print(pfname, '-dpng', '-r200')
      end
      fprintf(1,'%f, %d, %f, %f, %f, %f\n', Vc(m), N(m,ir), Ap(m,ir), Apm_mean(m,ir), Asm_mean(m,ir), Asm_mean(m,ir)/Apm_mean(m,ir))
   end
end
%% plot results
col = linspecer(4);
figure(3);clf
subplot(211)
plot([0 .1],[1 1],'--k')
hold on
for i=1:4
   h(i)=plot(Vc,1e-6*Apm_mean(:,i),'--','linewidth',2,'color',col(i,:));
   hold on
end
for i=1:4
  h(i+4) = plot(Vc,1e-6*Asm_mean(:,i),'-','linewidth',2,'color',col(i,:))
end
legend([h(5);h(6);h(7);h(8)],' 12 \mum',' 20 \mum',' 50 \mum','100 \mum','location','northwest')
ylim([0 3])
ylabel('Particle area A_p, Sampled area A_s (mm^2)')
subplot(212)
for i=1:4
hh(i)=plot(Vc,Asm_mean(:,i)./Apm_mean(:,i),'-','linewidth',2,'color',col(i,:));
hold on
end
ylabel('Ratio A_p / A_m')
xlabel('Volume concentration V_c')
legend(hh,' 12 \mum',' 20 \mum',' 50 \mum','100 \mum','location','southwest')

%%
figure(4); clf
subplot(211)
hold on
for i=1:4
   h(i)=plot(Vc,(1-1e-6*Asm_mean(:,i)),'-','linewidth',2,'color',col(i,:));
   hold on
end
legend([h(1);h(2);h(3);h(4)],' 12 \mum',' 20 \mum',' 50 \mum','100 \mum','location','northeast')

ylabel('Transmissivity (1/mm)')
subplot(212)
%plot([0 .1],[1 1],'--k')
hold on
for i=1:4
   h(i)=plot(Vc,1e-6*Asm_mean(:,i).*(1-1e-6*Asm_mean(:,i)),'-','linewidth',2,'color',col(i,:));
   hold on
end
ylabel('Backscatter response (mm)')
xlabel('Volume concentration V_c')






%ts = sprintf('V_c = %6.4f, r = %4.0f, N = %d, A_s/A_p = %5.2f',Vc(m), r, N(m), Asm_mean(m)/Apm_mean(m));





