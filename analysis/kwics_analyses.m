addpath ../setup
addpath ../../../MITgcm/utils/matlab
params = gendata_params();
kTopo = params.kTopo;
addpath(fullfile(getenv('HOME'),'MATLAB')); % for cmocean

for i = 1:length(kTopo)
    % calcFluxBC(kTopo(i));
end
% plotReson(kTopo)

load plotReson.mat aflxBC aflx kTopo
[~,idx] = max(aflx);
% plotVel(kTopo(idx));

% Possible mode 2 peak?
% plotVel(kTopo(3));

% At ~600km corrugations
plotVel(kTopo(26));

% At 700km corrugations
plotVel(kTopo(30));
