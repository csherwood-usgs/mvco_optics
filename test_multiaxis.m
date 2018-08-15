% test dual axis

z = linspace(0, 5, 10);
a = z.^2;
b = z.*1.5 + 7;

figure(1); clf
plot(z,a)
plot(z,a,'ob')
ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
plot(z,b,'Parent',ax2)
plot(z,b,'or')