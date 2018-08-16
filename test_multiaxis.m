% test dual axis


%%
figure(3);clf
subplot(221)
x1 = 0:0.1:40;
y1 = 4.*cos(x1)./(x1+2);
line(x1,y1,'Color','r')
hold on
plot(x1,y1,'or')
ax1 = gca; % current axes
YLim1 = ax1.YLim
ax1.XColor = 'r';
ax1.YColor = 'r'
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none')
 
x2 = 1:0.2:20;
y2 = x2.^2./x2.^3;
line(x2,y2,'Parent',ax2,'Color','k')
hold on
plot(x2,y2,'+k')
set(gca,'YLim',YLim1);
