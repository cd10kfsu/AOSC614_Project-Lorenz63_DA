clear; close all;
prefix = ['fort']; %ensrf_m3_infl08_oberr5 letkf_m3_infl08_oberr5
inflation = 0.8;
nmember   = 4;
method        = '';

xlim = [2900 3500];
ylim = [0 25];


xb = load([prefix '.10030']);
xo = load([prefix '.10010']);
xa = load([prefix '.10040']);
xt = load([prefix '.10020']);
%gr = load([prefix '.10070']);

%sprd2_his = load([prefix '.10050']);
%sprd2_his = load([prefix '.10065']);

% rmse
dat = xa - xt;
dbt = xb - xt;

rmseb = sqrt( sum(dbt'.^2)/3 );
rmsea = sqrt( sum(dat'.^2)/3 );


figure;
%--------------------------------------------------------------------------
% x:
% true trajectory (green), backround (red), analysis (blue)
%--------------------------------------------------------------------------
subplot(4,1,1);
hold on;
set(gca,'xlim',xlim);
plot( xt(:,1), 'g','linewidth',2);
plot( xb(:,1), 'r','linewidth',2);
plot( xa(:,1), 'b','linewidth',2);

set(gca,'box','on','fontsize',20);
%set(gca,'ylim',[-30 30]);



h=title([method ' : x']);
set(h,'fontsize',22);


%--------------------------------------------------------------------------
% y:
% true trajectory (green), backround (red), analysis (blue)
%--------------------------------------------------------------------------
subplot(4,1,2);
hold on;
set(gca,'xlim',xlim);
plot( xt(:,2), 'g','linewidth',2);
plot( xb(:,2), 'r','linewidth',2);
plot( xa(:,2), 'b','linewidth',2);

set(gca,'box','on','fontsize',20);
%set(gca,'ylim',[-30 30]);

h=title([method ' : y']);
set(h,'fontsize',22);


%--------------------------------------------------------------------------
% z:
% true trajectory (green), backround (red), analysis (blue)
%--------------------------------------------------------------------------
subplot(4,1,3);
hold on;
set(gca,'xlim',xlim);
plot( xt(:,3), 'g','linewidth',2);
plot( xb(:,3), 'r','linewidth',2);
plot( xa(:,3), 'b','linewidth',2);

set(gca,'box','on','fontsize',20);
%set(gca,'ylim',[0 50]);


h=title([method ' : z']);
set(h,'fontsize',22);

%--------------------------------------------------------------------------
% RMSE of background & analysis
% the equation is RMSE = sqrt[ (x-xt)^2 + (y-yt)^2 +(z-zt)^2 ]
%--------------------------------------------------------------------------


subplot(4,1,4);
hold on;
set(gca,'xlim',xlim);
plot(rmseb,'r','linestyle','-','marker','.','linewidth',2);
plot(rmsea,'b','linestyle','-','marker','.','linewidth',2);

set(gca,'box','on','fontsize',20);
title('RMSE');
set(gca,'ylim',[0 20]);

mean(rmsea(3000:end))

% LETKF: 0.3469
% 3DVAR: 0.7833
% hybrid: 
   
