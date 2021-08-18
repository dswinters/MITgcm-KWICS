function [aflxBC aflx] = plotReson(kTopo)

fs = 14; fn = 'times';
lam = 2*pi./kTopo;
dx_inner = 5e3; %2.e3;
dy_inner = 5e3; %2.e3;
params = gendata_params();
om = params.om;

datt = [];
aflxBC = [];
aflx = [];
yf = [];

for j = 1:length(kTopo)
    % Load run data
    kTopoPrefix = sprintf('kTopo%.8f',kTopo(j)); % File prefix for kTopo
    rdir = '..';
    rname = sprintf('run_%s',kTopoPrefix);
    load(fullfile(rdir,'setup','corrugation_params.mat'),'xSin1')
    datt = rdmnc(fullfile(rdir,'runs',rname,'outs_sn.*'),'T','iter');
    fluxFile = fullfile(rdir,'analysis','calcFluxBC',[rname '_flux_wave_avg_end.mat']);

    if exist(fluxFile,'file')
        fl = load(fluxFile);
        if j==1
            gridm = rdmnc(fullfile(rdir,'runs',rname,'grid*'));
            yf = fl.yc(find(gridm.Depth(1,:)==0,1,'last')) + [0 100e3];
        end
        indy = fl.yc >= yf(1) & fl.yc <= yf(2);

        % X location where we computed fluxes (70km after end of sine topog)
        [~,indx] = min(abs(fl.xc-xSin1+70e3));

        aflxBC(j) = squeeze(nansum(fl.upbcNz(indx,indy,:)*dy_inner,2));
        aflx(j) = squeeze(nansum(fl.upz(indx,indy,:)*dy_inner,2));
    else
        aflxBC(j) = nan;
        aflx(j) = nan;
    end
    fprintf('\r%s [%d of %d]',fluxFile,j,length(kTopo))
end

fprintf('\n')
save plotReson.mat aflxBC aflx indx datt yf

figure('position',[0 0 1065 400],'color','w','paperpositionmode','auto');

% subplot(2,2,3)
pl1 = [];
pl2 = [];
cols = get(groot,'DefaultAxesColorOrder');
clf

subplot(1,2,1); hold on
pl1 = plot(lam/1e3,aflxBC(end,:),'-','linewidth',1,'color',cols(1,:)); grid on
xlabel('Corrugation lengthscale (km)')
ylabel('Along-shelf baroclinic flux to 100 km (W/m)')
set(gca,'fontsize',8)

% subplot(2,2,4)
subplot(1,2,2); hold on
pl2 = plot(lam,aflx(end,:),'-','linewidth',1,'color',cols(1,:)); grid on
xlabel('Corrugation lengthscale (km)')
ylabel('Along-shelf flux to 100 km (W/m)')
set(gca,'fontsize',8)


print('-dpng','kwics_resonance.png')
save plotReson.mat aflxBC aflx kTopo
