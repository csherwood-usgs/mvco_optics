% p_fit_mo_times - Fit settling velocities for several quantities at
% specific times
ifnorm = 0; % divide all of the values by the global mean

avinterval = '20'

load suspsed_ba20_20rstrim_crs_cen
% use bin centers when available
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
        'UMEtc','UMEzc','UMEgamma','gamma';...
        };

za = 0.1; % (m) Ca estimates standardized to this elevation
nav = 0;  % number of profiles to left and right to be averaged
%%
ino = 2;

ydt = 276.; % list of target year days
dnt = datenum(2011,1,1)+ydt-1
datestr(dnt)

figure(1); clf
% Cp650
plotnum=0;
for ino=[7]
    plotnum = plotnum+1
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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

% LISST Vol. conc - micro and macro flocs
plotnum=3;
for ino=[3:4]
    plotnum = plotnum+1;
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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

% LISST Vol. conc total
if(0)
plotnum=5;
for ino=[1]
    plotnum = plotnum+1;
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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
end
% ABSS profiles
plotnum=6;
for ino=[13:15]
    plotnum = plotnum+1
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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
% LISST Vol. conc
plotnum=3;
for ino=[2:4]
    plotnum = plotnum+1;
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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
% LISST attenuation
plotnum=2;
for ino=[5]
    plotnum=plotnum+1;
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
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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
% UME gamma
plotnum=9;
for ino=[16]
    plotnum=plotnum+1;
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
    
    if(ino==8 || ino==16) % no variance for Chl or gamma
        allvc=ones(size(allc))
    else
        allvc=vv.(inst{ino,3});
    end
    dz = diff(allz(1:2,1));
    [ nz nt ]= size( allz )
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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
% UME chlorophyll
plotnum=10;
for ino=[8]
    plotnum=plotnum+1;
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
    
    if(ino==8 || ino==16) % no variance for Chl or gamma
        allvc=ones(size(allc))
    else
        allvc=vv.(inst{ino,3});
    end
    dz = diff(allz(1:2,1));
    [ nz nt ]= size( allz )
    
    % find target time
    i = find(allt(6,:)>dnt,1,'first');
    
    ic = 1;
    
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
        subplot(4,3,plotnum)
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