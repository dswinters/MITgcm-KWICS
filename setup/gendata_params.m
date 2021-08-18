% Use this function to get setup parameters to avoid hard-coding parameters in
% multiple places.

function params = gendata_params()

params.lprof = 20e3; % depth (negative) or extent (positive) of corrugations
params.lTopo = linspace(20e3,700e3,30);
params.kTopo = 2*pi./params.lTopo;
params.f = 1e-4;
params.g = 9.81;
params.om = 1.36*params.f;
params.deltaT = 500;

% params.lTopo = params.lTopo(7);
% params.kTopo = params.kTopo(7);
