% cluster_prep.m
% make matrices suitable for k-means clustering

% Script to plot profiles. Version 3 was made by EB and the plotting was
% modified by CRS. This version is changed as follow:
% Parameters to be plotted are extracted from the database with names at
% the top. This is a waste of memory, but makes the plotting code easier to
% read and understand.
close all
clear all

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
%% UMEattn650 profiles
allt=ba.UMEtc;
allz=ba.UMEzc;
allc=ba.UMEattn650;
t1 = allt(3:9,1:2031);
z1 = allz(3:9,1:2031);
c1 = allc(3:9,1:2031);

% LISST D50v
allt=ba.LISSTtc;
allz=ba.LISSTzc;
allc=ba.LISSTD50v;
t2 = allt(3:9,1:2031);
z2 = allz(3:9,1:2031);
c2 = allc(3:9,1:2031);
% abss
allt=ba.absstc;
allz=ba.absszc;
allc=ba.abss2;
t3 = allt(3:9,1:2031);
z3 = allz(3:9,1:2031);
c3 = allc(3:9,1:2031);

good_columns = find(~isnan(sum([t1;t2;t3;z1;z2;z3;c1;c2;c3])));
c1 = c1(:,good_columns);
c2 = c2(:,good_columns);
c3 = c3(:,good_columns);

t = mean([t1(:,good_columns);t2(:,good_columns);t3(:,good_columns)],1);
z = mean([z1(:,good_columns),z2(:,good_columns),z3(:,good_columns)],2);
size(c1)
size(c2)
size(c3)
size(z)
size(t)
pcolorjw(t,z,c3)
% 