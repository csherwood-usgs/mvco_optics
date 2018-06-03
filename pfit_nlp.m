function pf=pfit_nlp(c,z,varc,iplot,za)
% % Non-linear fit to Rouse profile
% % Based on fit_lambda by Pat Dickhudt
% 
if(exist('iplot','var')~=1),iplot=0,end
if(exist('za','var')~=1), za=z(1); end
ok = (~isnan(c+z+varc)&c>0&z>0&varc>0) 
y = c(ok);
x = z(ok)./za;
vary = varc(ok);
disp([x y vary])
if(length(y)<3), error('Need at least 3 points in pfit_nl'), end

% fit exponential function
[efit,model] = expfit(x,y,[x(1) -4]);

[sse, yhat] = model(efit)
sst = sum( (yhat- mean(yhat) ).^2);
r2 = 1-sse./sst;

zest = logspace( log10(min(z)),log10(max(z)), 20);
Cest = efit(1).*(zest./za).^efit(2);
if(iplot)
   hh=loglog(c,z,'ob');
   set(hh,'color',[ .4 .4 .4])
   hold on
   for i=1:length(x)
      hh = loglog([y(i)-sqrt(vary(i)) y(i)+sqrt(vary(i))],za*[x(i) x(i)],'-k')
   set(hh,'color',[.4 .4 .4],'linewidth',2)
   end
   loglog(Cest,zest,'-r')
   hh=loglog(y,za*x,'ob');
   set(hh,'markerfacecolor',[1 .4 .4],'markersize',10)
end
pf.N = length(c)
pf.p = efit(2);
pf.Ca = efit(1);
pf.sa = NaN;
pf.sb = NaN
pf.r2 = r2;

   function [estimates, model] = expfit(x,y,start_point)
      % exponential fit function using fminsearch
      % copied from Matlab help: Curve Fitting via Optimization
      
      % Call fminsearch with a random starting point.
      %start_point = [.005 .0005];
      model = @expfun;
      options = optimset('TolFun',1E-6,'TolX',1E-6,'MaxFunEvals',1000,'MaxIter',1000);
      [estimates,fval,exitflag] = fminsearch(model,start_point,options);
      
      % expfun accepts curve parameters as inputs, and outputs sse,
      % the sum of squares error for
      % and the FittedCurve. FMINSEARCH only needs sse, but we want
      % to plot the FittedCurve at the end.
      
      function [sse, yhat] = expfun(params)
         yhat = params(1) .* x.^params(2);
         sse = sum( (yhat - y).^2 );
      end
   end
end
