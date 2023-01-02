clear; close all;
prefix = '3dvar'; %['4dvar_out2_obs2']; %ensrf_m3_infl08_oberr5 letkf_m3_infl08_oberr5
%xlim = [12000 16000];
xlim = [3000 4000];
ylim = [0 25];

xa = load([prefix '.10040']);
xt = load([prefix '.10020']);
dat = xa - xt;
rmsea = sqrt( sum(dat'.^2)/3 );
tmp = rmsea;
nsmooth = 20;
for i = 1+nsmooth:length(tmp)-nsmooth
    rmsea(i) = sum(tmp(i-nsmooth:i+nsmooth))/(2*nsmooth+1);
end

prefix = ['4dvar_out3_obs2']; %ensrf_m3_infl08_oberr5 letkf_m3_infl08_oberr5
xa = load([prefix '.10040']);
xt = load([prefix '.10020']);
dat = xa - xt;
rmsea2 = sqrt( sum(dat'.^2)/3 );
tmp = rmsea2;
for i = 1+nsmooth:length(tmp)-nsmooth
    rmsea2(i) = sum(tmp(i-nsmooth:i+nsmooth))/(2*nsmooth+1);
end

prefix = ['4dvar_out3_obs3']; %ensrf_m3_infl08_oberr5 letkf_m3_infl08_oberr5
xa = load([prefix '.10040']);
xt = load([prefix '.10020']);
dat = xa - xt;
rmsea3 = sqrt( sum(dat'.^2)/3 );

tmp = rmsea3;
for i = 1+nsmooth:length(tmp)-nsmooth
    rmsea3(i) = sum(tmp(i-nsmooth:i+nsmooth))/(2*nsmooth+1);
end

prefix = ['4dvar_out3_obs4']; %ensrf_m3_infl08_oberr5 letkf_m3_infl08_oberr5
xa = load([prefix '.10040']);
xt = load([prefix '.10020']);
dat = xa - xt;
rmsea4 = sqrt( sum(dat'.^2)/3 );

tmp = rmsea4;
for i = 1+nsmooth:length(tmp)-nsmooth
    rmsea4(i) = sum(tmp(i-nsmooth:i+nsmooth))/(2*nsmooth+1);
end





figure;
hold on;
%et(gca,'xlim',[1000 1500]);
plot(rmsea,'k','linestyle','-','marker','.','linewidth',2);
plot(rmsea2,'color',[ 1 0 0],'linestyle','-','marker','.','linewidth',2);
plot(rmsea3,'color',[ 0.7 0 0],'linestyle','-','marker','.','linewidth',2);
plot(rmsea4,'color',[ 0 1 0],'linestyle','-','marker','.','linewidth',2);

set(gca,'box','on','fontsize',20);
title('RMSE');
%set(gca,'ylim',[0 20]);

mean(rmsea(3000:end))
mean(rmsea2(3000:end))
mean(rmsea3(3000:end))
mean(rmsea4(3000:end))


% LETKF: 0.3469
% 3DVAR: 0.7833
% hybrid: 
   
