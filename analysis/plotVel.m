function plotVel(kTopo)
lTopo = 2*pi/kTopo;
params = gendata_params();
om = params.om;
fs = 8; fn = 'times';
kTopoPrefix = sprintf('kTopo%.8f',kTopo); % File prefix for kTopo
rname = sprintf('run_%s',kTopoPrefix);
froot = fullfile('..','runs',rname); % FIXME revert this
% froot = '/home/dw/Work/MITgcm/simulations/run_kTopo0.00002720/';
fig_dir = sprintf('plotVel/kTopo%.8f',kTopo);
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% Load grid and time info
gridm = rdmnc(fullfile(froot, 'grid*'));
datt = rdmnc(fullfile(froot,'outs_sn.*'),'T','iter');
files = dir(fullfile(froot,'outs_sn.*.nc'));         % all files
fids = extractBetween({files.name},'outs_sn.','.t'); % time identifiers (ignore tile suffixes)
fids = unique(fids);
yf = gridm.Y(find(gridm.Depth(1,:)==0,1,'last')) + [0 100e3];

% Load corrugation parameters (contains location of flux line)
load('../setup/corrugation_params.mat','xSin1');
[~,nx] = min(abs(gridm.Xp1-xSin1));


figure('position',[0 0 1000 800])
cm = cmocean('balance',101);

u0 = 0.1; % forcing amplitude
uclim = 0.5*u0;
vclim = 0.1*u0;

aspect_1to1 = true;
xl = [30e6 35e6];
xt = 30e6:.5e6:35e6;

% Load flux file
ffile = fullfile('calcFluxBC', [rname '_flux_wave_avg_end.mat']);
flux = load(ffile);
flux_nskip = 5;
flux_ix = find(flux.xc > 31e6,1,'first'):flux_nskip:find(flux.xc < 33e6,1,'last');
flux_iy = 1:flux_nskip:find(flux.yc < 600e3,1,'last');

for i = 1:length(datt.iter)
    dat = load_data(froot,fids,datt.iter(i),{'UVEL','VVEL'});
    clf
    tl = tiledlayout(3,2);

    % First plot: u plan-view
    nexttile([1,2])
    pcolor(gridm.Xp1,gridm.Y,squeeze(dat.UVEL(:,:,1))'); hold on
    title('u (m/s)')
    caxis([-1 1]*uclim);
    %
    colorbar
    shading flat
    colormap(gca,cm)
    plot(xSin1*[1 1]+70e3,yf,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    ylim([0 1e6]);
    xticklabels(xticks/1e3);
    yticklabels(yticks/1e3);
    contour(gridm.XC',gridm.YC',gridm.Depth',[0:500:4e3],'k');
    contour(gridm.XC',gridm.YC',gridm.Depth',[0 0],'color','k','linewidth',2);
    if aspect_1to1
        axis equal
        xlim(xl)
        ylim([0 1e6])
        xticks(xt);
        xticklabels(xticks/1e3);
    end
    hq = quiver(flux.xc(flux_ix), flux.yc(flux_iy),...
                flux.upbcNz(flux_ix, flux_iy)',...
                flux.vpbcNz(flux_ix, flux_iy)',...
                'k');

    % Second plot: v plan-view
    nexttile([1,2])
    pcolor(gridm.X,gridm.Yp1,squeeze(dat.VVEL(:,:,1))'); hold on
    title('v (m/s)')
    caxis([-1 1]*vclim);
    %
    colorbar
    shading flat
    colormap(gca,cm)
    plot(xSin1*[1 1]+70e3,yf,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    ylim([0 1e6]);
    xticklabels(xticks/1e3);
    yticklabels(yticks/1e3);
    contour(gridm.XC',gridm.YC',gridm.Depth',[0:500:4e3],'k');
    contour(gridm.XC',gridm.YC',gridm.Depth',[0 0],'color','k','linewidth',2);
    if aspect_1to1
        axis equal
        xlim(xl)
        ylim([0 1e6])
        xticks(xt);
        xticklabels(xticks/1e3);
    end
    hq = quiver(flux.xc(flux_ix), flux.yc(flux_iy),...
                flux.upbcNz(flux_ix, flux_iy)',...
                flux.vpbcNz(flux_ix, flux_iy)',...
                'k');


    % Third plot: u slice at flux calculation line
    nexttile
    pcolor(gridm.Y,gridm.Z,squeeze(dat.UVEL(nx,:,:))'); hold on
    title('u (m/s)')
    caxis([-1 1]*uclim);
    %
    plot(yf(end)*[1 1],ylim,'k--')
    colorbar
    shading flat
    colormap(gca,cm)
    plot(gridm.Y,-gridm.Depth(end,:),'k-');
    xlim([0 500e3]);
    xlabel('y (km)')
    ylabel('z (km)')
    xticklabels(xticks/1000)
    yticklabels(yticks/1000)

    % Fourth plot: v slice at flux calculation zone
    nexttile
    pcolor(gridm.Yp1,gridm.Z,squeeze(dat.VVEL(nx,:,:))'); hold on
    title('v (m/s)')
    caxis([-1 1]*vclim);
    %
    plot(yf(end)*[1 1],ylim,'k--')
    colorbar
    shading flat
    colormap(gca,cm)
    plot(gridm.Y,-gridm.Depth(end,:),'k-');
    xlim([0 500e3]);
    xlabel('y (km)')
    ylabel('z (km)')
    xticklabels(xticks/1000)
    yticklabels(yticks/1000)

    ttxt = sprintf('kTopo=%.8f, lTopo=%.2f | T=%.1f cycles | [nx ny nz]=[%d %d %d]',...
                   kTopo,lTopo,datt.T(i)/(2*pi/om),...
                   length(gridm.X),length(gridm.Y),length(gridm.Z));
    hax = axes('visible','off','position',[0 0 1 1]);
    xlim(hax,[-1 1]); ylim(hax,[-1 1]);
    text(hax,0,1,ttxt,...
         'fontweight','bold',...
         'verticalalignment','top',...
         'horizontalalignment','center');

    fout=fullfile(fig_dir,sprintf('vel_%04d',i));
    print('-djpeg90','-r300',fout)
    fprintf('\rSaved %s [%0.2f%%]',fout,100 * i / length(datt.iter));
end




function [dat] = load_data(froot,fids,iter,vars)
    % If the simulation is too large for a single .nc file, rdmnc fails when
    % trying to load "outs_sn.*.nc". This workaround keeps trying each file
    % until the one containing the desired iteration is found.
    % TODO: There's probably a better way to do this.
    fileFound = false;
    nf = 0;
    while ~fileFound & nf < length(fids)
        nf = nf + 1; % try next file
        try
            dat = rdmnc(fullfile(froot,['outs_sn.' fids{nf},'.t*.nc']),vars{:},iter);
            fileFound = true;
        catch err
        end
    end
    if ~fileFound
        error('Failed loading data from iter %d',iter)
    end
