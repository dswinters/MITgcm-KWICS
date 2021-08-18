%% calculate baroclinic energy flux
%% RCM Feb 2020
%% Modified by DW 2021
function outname = calcFluxBC(kTopo)

fs = 8; fn = 'times';
kTopoPrefix = sprintf('kTopo%.8f',kTopo); % File prefix for kTopo
rname = sprintf('run_%s',kTopoPrefix);
froot = fullfile('..','runs',rname);
fig_dir = 'calcFluxBC';
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

outname = fullfile(fig_dir, [rname '_flux_wave_avg_end.mat']);
if exist(outname,'file'); return; end % comment this to overwrite output

% Define some constants
descr = ['calcFluxBC ' rname];
rhoNil = 999.8;
tAlpha = 2e-4;
params = gendata_params();
g = params.g;
f = params.f;
om = params.om;

% Read grid data from model output
gridm = rdmnc(fullfile(froot, 'grid*'));
xc2d = gridm.XC; xc = squeeze(xc2d(:,1)); nx = length(xc);
xg2d = gridm.XG; xg = squeeze(xg2d(:,1)); Lx = max(xg);
dxc2d = gridm.dxC; dxc = squeeze(dxc2d(:,1));
dxg2d = gridm.dxG;  dxg = squeeze(dxg2d(:,1));
yc2d = gridm.YC; yc = squeeze(yc2d(1,:)); ny = length(yc);
yg2d = gridm.YG; yg = squeeze(yg2d(1,:)); Ly = max(yg);
dyc2d = gridm.dyC; dyc = squeeze(dyc2d(1,:));
dyg2d = gridm.dyG; dyg = squeeze(dyg2d(1,:));
rc = gridm.RC; nz = length(rc);
rf = gridm.RF;
drc = gridm.drC;
drf = gridm.drF;
hfacc = squeeze(gridm.HFacC);
hfacs = squeeze(gridm.HFacS);
hfacw = squeeze(gridm.HFacW);
raw = gridm.rAw;
ras = gridm.rAs;
rac = gridm.rA;
raz = gridm.rAz;
dpth = gridm.Depth;

% find unique output times
files = dir(fullfile(froot,'outs_sn.*.nc'));         % all files
fids = extractBetween({files.name},'outs_sn.','.t'); % time identifiers (ignore tile suffixes)
fids = unique(fids);
[~,fidx] = unique(fids);

% Read time data
datt = rdmnc(fullfile(froot,'outs_sn.*'),'T','iter');
nt = length(datt.T);

% get background pressure
dat = rdmnc(fullfile(froot,'outs_sn.0000000000.*.nc'),'PHIHYD',datt.iter(1));
p0 = dat.PHIHYD*rhoNil; % kgm^-3 * m^2/s^2 = pa

% initialize summed flux terms
up = zeros(nx,ny,nz);
vp = zeros(nx,ny,nz);
upbt_ = zeros(nx,ny);
vpbt_ = zeros(nx,ny);
upbcN = zeros(nx,ny,nz);
vpbcN = zeros(nx,ny,nz);
upbcK = zeros(nx,ny,nz);
vpbcK = zeros(nx,ny,nz);

% alias for vertically integrating an arbitrary variable
% This is equivalent to something like:
%    for jj = 1:nz
%        vint = vint + x(:,:,jj).*hfacc(:,:,jj)*gridm.drF(jj);
%    end
vint=@(x) sum(x.*hfacc.*permute(drf,[1,3,2]),3);

% number of timesteps in a wave period
nperiod = sum(datt.T >= datt.T(end)-2*pi/om);

% initialize depth- and time-averaved flux terms
upbt = zeros(nx,ny);
vpbt = zeros(nx,ny);
upz = zeros(nx,ny);
vpz = zeros(nx,ny);
upbcNz = zeros(nx,ny);
vpbcNz = zeros(nx,ny);
upbcKz = zeros(nx,ny);
vpbcKz = zeros(nx,ny);

for ii = [-nperiod+1:0] + length(datt.T) % sum over the final forcing period
    % load data and compute fluxes
    [pbt,ubt,vbt,pbc,ubc,vbc,ufull,vfull,panom] = ...
        compute_fluxes(froot,fids,datt.iter(ii),rhoNil,p0,hfacc,drf,dpth);

    % add to summed terms
    up = up + ufull;
    vp = vp + vfull;
    upbt_ = upbt_ + ubt.*pbt.*dpth; %  [W/m^2] this is K&F
    vpbt_ = vpbt_ + vbt.*pbt.*dpth;
    upbcN = upbcN + ubc.*pbc; %  Nash 05
    vpbcN = vpbcN + vbc.*pbc;
    upbcK = upbcK + ubc.*panom; %  K&F
    vpbcK = vpbcK + vbc.*panom;
end

% Vertically integrate and time-average the summed fluxes
upz = vint(up/nperiod);
vpz = vint(vp/nperiod);
upbt = upbt_/nperiod; % already depth-integrated
vpbt = vpbt_/nperiod; % already depth-integrated
upbcNz = vint(upbcN/nperiod);
vpbcNz = vint(vpbcN/nperiod);
upbcKz = vint(upbcK/nperiod);
vpbcKz = vint(vpbcK/nperiod);

save(outname,'xc','yc','upz','vpz','upbcNz','vpbcNz','upbt','vpbt')
fprintf('\rSaved %s\n',outname);

%% TODO: put figure generation back in here

% Load data and compute fluxes. Putting this in its own sub-function makes the
% code more flexible.
function [pbt,ubt,vbt,pbc,ubc,vbc,ufull,vfull,panom] = ...
    compute_fluxes(froot,fids,iter,rhoNil,p0,hfacc,drf,dpth)

    % If the simulation is too large for a single .nc file, rdmnc fails when
    % trying to load "outs_sn.*.nc". This workaround keeps trying each file
    % until the one containing the desired iteration is found.
    % TODO: There's probably a better way to do this.
    fileFound = false;
    nf = 0;
    vars = {'UVEL','VVEL','PHIHYD'};
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
    panom = dat.PHIHYD*rhoNil-p0;

    % cell-centered velocities
    uC = 0.5*(dat.UVEL(2:end,:,:)+dat.UVEL(1:end-1,:,:));
    vC = 0.5*(dat.VVEL(:,2:end,:)+dat.VVEL(:,1:end-1,:));

    % full flux
    ufull = uC.*panom;
    vfull = vC.*panom;

    % BC/BT fluxes
    % use BT/BC decomposition on stratified fields
    pbt = sum(panom.*hfacc.*permute(drf,[1,3,2]),3)./dpth;
    ubt = sum(uC.*hfacc.*permute(drf,[1,3,2]),3)./dpth;
    vbt = sum(vC.*hfacc.*permute(drf,[1,3,2]),3)./dpth;
    pbc = bsxfun(@minus,panom,pbt);  % this is Nash 05
    ubc = bsxfun(@minus,uC,ubt);
    vbc = bsxfun(@minus,vC,vbt);
