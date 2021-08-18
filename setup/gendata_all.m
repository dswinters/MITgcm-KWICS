clear all, close all

params = gendata_params();
kTopo = params.kTopo;

% Some flags so we can generate files on e.g. the first theta of every kTopo,
% only once per group of runs, etc...
output_freq = 0; % don't write output
niter0 = 0;

overwrite = true;
for j = 1:length(kTopo) % kTopo index
        kTopoPrefix = sprintf('kTopo%.8f_',kTopo(j)); % File prefix for kTopo

        % Set up run directories
        rname = sprintf('run_%s',kTopoPrefix(1:end-1));
        rdir = fullfile('..','runs',rname);
        if ~exist(rdir,'dir'); mkdir(rdir); end

        % Generate run data
        params = gendata(kTopo(j), rdir, overwrite);
        overwrite = false;
end
