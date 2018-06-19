% rouse - Integrate a Rouse profile
rhow = 1000 % water density kg m-3
nu = 1.36e-6 % kinematic viscosity m2 s-1
mu = rhow*nu % dynamic viscosity
g = 9.81; % m s-2
vk = 0.41;

% Set ustar
ustr = 0.02 % m/s

% number of size classes
nc = 1;

% elevations (m)
za = 0.02;
z1 = 2.5;
z = logspace( log10(za), log10(z1), 20 )';

% diameters (m)
nf = 2;
Dp0 = 4e-6;
Df = logspace(log10(4e-6), log10(1500e-6), nc)
rhos = 2650*ones(1,nc);
rhof = NaN*ones(1,nc);
for i=1:length(nc)
rhof(i) = max(1020,min(rhos(i),rhow + (rhos(i)-rhow)*(Df(i)./Dp0).^(nf-3)));
end

ws = -(2.*g.*(rhof-rhow).*(Df./2).^2)./(9*mu)
p = ws/(vk*ustr)
%% Calc log profiles
Ca = ones(size(Df))./nc;
Cm = NaN*ones(length(z),length(Df));

for i=1:length(nc)
   Cm(:,i) = Ca(i) .* (z./za).^p(i);
end
Cmt = sum(Cm,2);
%% plot
figure(1); clf
h1=loglog(Cm,z);
hold on
h2=loglog(Cmt,z,'linewidth',2);
xlabel('C')
ylabel('z')