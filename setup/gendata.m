function params = gendata(kTopo, rdir, overwrite)

params = gendata_params();
lprof = params.lprof; % depth of corrugations
shelf_offset = 2.1*abs(lprof);

%% setup modified from MITgcm-PWICS for Kelvin wave simulation

rname = sprintf('run_kTopo%.8f',kTopo);
prec='real*4';
ieee='b';
fs = 12; fn = 'times';

% % use the setup from Jim Lerczak's model to set bathymetry & stratification
% jl = load('linearModel/subcritTopo_modes.mat'); % warning: in sigma coordinates!

g = params.g;
f = params.f; % gsw_f(43.29);
om = params.om;
np = [2 2]; % processors [nx, ny]


% sponge cells on north, east and west boundaries
nsponge = 15;

% vertical grid parameters
Lz = 4000; % [m]

% this is a crude estimate for the KW neglecting the shelf
cp = sqrt(g*Lz);
k = om/cp;
lam = 2*pi/k;

Ld = cp/f; % offshore length scale

% for this setup, and sig = 1.36e-4; I get l = 7.008e-7 from Jim's model (look in
% /Users/ruth/research/software/lerczakScatterInvisc_ForPublic/dispFreeSurface/smallL for the
% resonance scan)
% The "crude" estimate works surprisingly well - 
% the modal structure of the homogeneous KW has a similar offshore lengthscale to that
% from the (smoothed) stratified model. It is much simpler to use this unstratified approximation
% to force the model, and allow it to evolve into something dynamically consistent within MITgcm.


% vertical grid, exponential high res near surface
nzc = 30;
mindz = 25; maxdz = 500;
dz = smooth([ones(1,floor(700/mindz))*mindz logspace(log10(mindz),log10(maxdz),nzc)],5)';
dz = smooth(dz,5)';
tmp = [cumsum(dz)];
ind = find(tmp<Lz);
dze = Lz-max(abs(tmp(ind)));  % make sure that depth is Lz (I am not totally sure this is ok)
dz = [dz(ind) dze];
zf = -[0 cumsum(dz)]; % this is RF
z = 0.5*(zf(1:end-1)+zf(2:end));
nzc = length(dz);


%horizontal grid parameters
dx_outer = lam/10; % set the outer grid resolution to resolve the forcing wave
dy_outer = dx_outer;
dx_inner = 5e3; %2.e3;
dy_inner = 5e3; %2.e3;

Lx = 2.*lam + 2*nsponge*dx_outer;
Ly = 2.5*Ld; % in m

LF = 1200e3; % nominal length of corrugated region

high_res_pad = 400e3;

x0inner = 0.7*Lx - LF/2 - high_res_pad; % let's start the high res region just before the sine topo
x1inner = 0.7*Lx + LF/2 + high_res_pad; % how long after sine topo should we have high res?

% now make horizontal grids
% x
dx = [dx_outer*ones(1,ceil(x0inner/dx_outer)) ...
      dx_inner*ones(1,ceil((x1inner-x0inner)/dx_inner))...
      dx_outer*ones(1,ceil((Lx-x1inner)/dx_outer))];
% add cells on right to make domain divisible by # of processors
n_add = np(1) - mod(length(dx),np(1));
dx = [dx, dx_outer*ones(1,n_add)];

% smooth dx
dx = smooth(smooth(dx,5),5)';
xg = [0 cumsum(dx)];
xc = 0.5*(xg(2:end)+xg(1:end-1));

% y
dy = [dy_inner*ones(1,floor((0.25*Ly)/dy_inner)) ...
      dy_outer*ones(1,floor((0.75*Ly)/dy_outer))];
% add cells offshore to make domain divisible by # of processors
n_add = np(2) - mod(length(dy),np(2));
dy = [dy, dy_outer*ones(1,n_add)];

% smooth dy
dy = smooth(smooth(dy,5),5)';
yg = [0 cumsum(dy)];
yc = 0.5*(yg(2:end)+yg(1:end-1));

nxc = length(xc); nyc = length(yc);
Lx = max(xc); Ly = max(yc);

%% Topography

lTopo = 2*pi/kTopo;

xSin1 = x1inner - 500e3; % end corrugations at some constant offset from high-res region
nLambdas = floor(LF/lTopo); % number of full wavelengths of corrugations
xSin0 = xSin1 - lTopo*nLambdas; % start sine wave around here
if overwrite
    save corrugation_params.mat xSin1 xSin0 nLambdas
end

hsh = 250; ysh = 75e3; dysl = 25e3;
prof = -hsh -0.5*(Lz-hsh)*(1+tanh((yc-ysh)/dysl));
%prof(1) = 0; % vertical wall
PROF = NaN(nxc,nyc);

cutoff = 0*xc;
cutoff(xc >= xSin0 & xc <= xSin1) = 1;
for ii = 1:nxc
    ycshift = yc + shelf_offset + lprof*(cutoff(ii)* (0.5 - 0.5*cos((xc(ii)-xSin1)*kTopo)));
    tmp = interp1(ycshift,prof,yc);
    PROF(ii,:) = tmp;
end

bads = find(isnan(PROF)); % these are where min(ycshift) was > 0
PROF(bads) = 0;


%% now stratification - from Jim's linear code
rho0 = 999.8; g = 9.81; alpha = 2e-4;

r1 = 992; r2 = 995;
r0 = (r1+r2)/2; dr = r2-r1;
N2back = (2*pi/(0.5*3600))^2;
mupyc = 400;
Zpyc = -400;
r = r2 - 0.5*dr*(1+tanh((z-Zpyc)/mupyc)) - z*N2back*r0/g;
n2 = 0.5*(dr/r0)*(g/mupyc)*sech((z-Zpyc)/mupyc).^2 + N2back;

t = (1-r/rho0)/alpha+5;


%% forcing region

x0mask = dx_outer*nsponge;
x1mask = x0mask + lam;

% rbcs mask: forcing will be applied within this region
% [~,xi0] = min(abs(xc-x0mask)); 
% [~,xi1] = min(abs(xc-x1mask));
% xind = xi0:xi1;
% maskxy = zeros(nxc,nyc);
% maskxy(xi0:xi1,:) = 1;
% mask = repmat(maskxy,[1 1 nzc]);


maskxy = zeros(nxc,nyc);
maskxy(xc>x0mask&xc<x1mask,:) = 1;
maskxy = conv2(maskxy,ones(4,1)/4,'same');
maskxy = conv2(maskxy,ones(4,1)/4,'same');
maskxy = conv2(maskxy,ones(4,1)/4,'same');


mask = repmat(maskxy,[1 1 nzc]);

%% rbcs forcing: modal structures to be read into rbcs_fields_load

% old way
% kb = load('linearModel/mode0_uvwrp.mat'); % KW from Ken's model, high res but only to 400km offshore
% kb2 = load('linearModel/mode0_uvwrp3.mat'); % as above, but lower res extend to 5000km offshore

% % merge the two model fields together
% indx = find(kb2.xpl>kb.xpl(end));
% xpl = cat(2,kb.xpl,kb2.xpl(indx));
% vg = cat(1,kb.vg,0.2395*kb2.vg(indx,:)); % youch!

% % interpolate onto MITgcm grid
% umode = interp2(xpl'*1e3,kb.z,vg',yc',z);
% umode = -umode./max(abs(umode(:))); % normalize
% bads = find(isnan(umode));
% umode(bads) = 0;

% %umode = repmat(exp(-(yc-yc(nsponge))/Ld),nzc,1); % modal structure for a vertical wall
% UMODE = permute(repmat(umode,[1 1 nxc]),[3 2 1]);

% stupid but better new way
vg = exp(-yc/Ld);
UMODE = repmat(vg,[nxc 1 nzc]);

diX = xSin0 - x1mask; % distance between forcing region and corrugation start
Tend = diX/cp + 5*(2*pi/om); 

%% output some parameters
params.t_end = Tend;
params.t_chkpt = Tend - 1.5*(2*pi/om); % time of checkpoint file
params.nxc = nxc;
params.nyc = nyc;
params.nzc = nzc;
params.np = np;


%% strat done, grids done, topo done

% initial fields
T = permute(repmat(t',[1 nxc nyc]),[2 3 1]);

U = 0*T;
V = 0*T;

% boundary sponge region fields

Uzonal = zeros(nyc,nzc,2); % [ny nz nt]
Vzonal = zeros(nyc,nzc,2);
Tzonal = repmat(squeeze(T(1,:,:)),[1 1 2]);

Umerid = zeros(nxc,nzc,2); % [nz nx nt]
Vmerid = zeros(nxc,nzc,2);
Tmerid = repmat(squeeze(T(:,end,:)),[1 1 2]);


%% write some files
if overwrite
    fid = fopen('inFiles/Uinit.bin','w',ieee);
    fwrite(fid,U,prec);
    fclose(fid);

    fid = fopen('inFiles/Vinit.bin','w',ieee);
    fwrite(fid,V,prec);
    fclose(fid);

    fid = fopen('inFiles/Tinit.bin','w',ieee);
    fwrite(fid,T,prec);
    fclose(fid);

    fid=fopen('inFiles/delZ.bin','w',ieee);
    fwrite(fid,dz,prec);
    fclose(fid);

    fid=fopen('inFiles/delX.bin','w',ieee);
    fwrite(fid,dx,prec);
    fclose(fid);

    fid=fopen('inFiles/delY.bin','w',ieee);
    fwrite(fid,dy,prec);
    fclose(fid);

    fid=fopen('inFiles/Umerid.bin','w',ieee);
    fwrite(fid,Umerid,prec);
    fclose(fid);

    fid=fopen('inFiles/Vmerid.bin','w',ieee);
    fwrite(fid,Vmerid,prec);
    fclose(fid);

    fid=fopen('inFiles/Tmerid.bin','w',ieee);
    fwrite(fid,Tmerid,prec);
    fclose(fid);

    fid=fopen('inFiles/Uzonal.bin','w',ieee);
    fwrite(fid,Uzonal,prec);
    fclose(fid);

    fid=fopen('inFiles/Vzonal.bin','w',ieee);
    fwrite(fid,Vzonal,prec);
    fclose(fid);

    fid=fopen('inFiles/Tzonal.bin','w',ieee);
    fwrite(fid,Tzonal,prec);
    fclose(fid);

    fid=fopen('inFiles/mask.bin','w',ieee);
    fwrite(fid,mask,prec);
    fclose(fid);

    fid=fopen('inFiles/umode.bin','w',ieee);
    fwrite(fid,UMODE,prec);
    fclose(fid);
end

% Write topography every time
fid = fopen(fullfile(rdir,'topog.bin'),'w',ieee);
fwrite(fid,PROF,prec);
fclose(fid);

% Make setup figure
figure('visible','on')
set(gcf,'position',[39 523  1736  439],'color','w')
subplot('position',[0.08 0.1 0.15 0.8])
plot(yc/1e3,prof)
title('Bathymetry')
set(gca,'fontsize',fs,'fontname',fn)

subplot('position',[0.3 0.1 0.15 0.8])
plot(t,z)
grid on; hold on
plot(t,z,'.')
title('Temperature')
set(gca,'fontsize',fs,'fontname',fn)

subplot('position',[0.48 0.3 0.05 0.6])
plot(dy/1e3,yc/1e3)
grid on; hold on
plot(dy/1e3,yc/1e3,'.')
ylim([0 Ly/1e3])
set(gca,'fontsize',fs,'fontname',fn,'layer','top')
ylabel('y [km]')

subplot('position',[0.55 0.3 0.15 0.6])
set(gcf,'paperpositionmode','auto','color','w')
pcolor(xc,yc,mask(:,:,1)'); shading flat;
cb = colorbar;
set(cb,'position',[0.71 0.3 0.015 0.6])
hold on
contour(xc,yc,PROF',[-4000:500:0],'k')
plot([xc(nsponge) xc(nsponge)],[yc(1) yc(end)],'k')
plot([xc(nxc-nsponge) xc(nxc-nsponge)],[yc(1) yc(end)],'k')
plot([xc(1) xc(end)],[yc(nyc-nsponge) yc(nyc-nsponge)],'k')
contour(xc,yc,PROF',[0 0],'linewidth',2,'color','k');
yf = shelf_offset + [0 100e3]; % y-coords of flux calc region
plot((xSin1*[1 1]+70e3),yf,'r-','linewidth',1);
caxis([-1 1])
grid on
xlim([0 Lx])
ylim([0 Ly])
xticklabels(xticks/1e3)
yticklabels([])
title('MASK')
set(gca,'fontsize',fs,'fontname',fn,'layer','top','yticklabel',{})

subplot('position',[0.55 0.1 0.15 0.1])
plot(xc/1e3,dx/1e3)
hold on
plot(xc/1e3,dx/1e3,'.')
plot(xSin0/1e3*[1 1],ylim,'k--')
plot(xSin1/1e3*[1 1],ylim,'k--')
xlim([0 Lx/1e3])
grid on
set(gca,'fontsize',fs,'fontname',fn)
xlabel('x [km]')
ylabel('dx [km]')

subplot('position',[0.77 0.1 0.2 0.8])
set(gcf,'paperpositionmode','auto','color','w')
pcolor(yc/1e3,z,squeeze(UMODE(1,:,:,1))'); shading flat; colorbar
hold on
plot(yc/1e3,prof,'k')
caxis([-1 1])
title('along shore modal structure')
set(gca,'fontsize',fs,'fontname',fn)
xlabel('y [km]')

print('-dpng',fullfile(rdir,'setup.png'))
