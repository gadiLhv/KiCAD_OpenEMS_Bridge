clc
clear
close all

load('simParams.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize openEMS environment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_CSX = 'openEMS_sim.xml';

load(fullfile(Sim_Path,'model_data.mat'));
load(fullfile(Sim_Path,'model_ports.mat'));
load(fullfile(Sim_Path,'mesh_data.mat'));
load(fullfile(Sim_Path,'openEMS_data.mat'));


pkg load geometry
pkg load miscellaneous

addpath(modelerPath);
addpath(miscPath);
addpath(openEMSpath);

physical_constants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate radiation addition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Currently add wl/2
fmin = f_center - f_cutoff;
wl_max = c0/fmin;

bBoxAdd = wl_max*airbox_pad_wl/unit;


%%%%%%%%%%%%%%%%%%%%%%%%
% Add air box and mesh %
%%%%%%%%%%%%%%%%%%%%%%%%
step_max_air = (1/unit)*airboxMaxStep_wl*(c0/fmax);

meshLeft = AutoSmoothMeshLines([(allMesh.x(1) - [bBoxAdd 0]) allMesh.x(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshRight = AutoSmoothMeshLines([allMesh.x(end-1) (allMesh.x(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshDown = AutoSmoothMeshLines([(allMesh.y(1) - [bBoxAdd 0]) allMesh.y(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshUp = AutoSmoothMeshLines([allMesh.y(end-1) (allMesh.y(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshBot = AutoSmoothMeshLines([(allMesh.z(1) - [bBoxAdd 0]) allMesh.z(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshTop = AutoSmoothMeshLines([allMesh.z(end-1) (allMesh.z(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);

allMesh.x = [meshLeft(1:(end-2)) allMesh.x meshRight(3:end)];
allMesh.y = [meshDown(1:(end-2)) allMesh.y meshUp(3:end)];
allMesh.z = [meshBot(1:(end-2)) allMesh.z meshTop(3:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine mesh limitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = DefineRectGrid( CSX, unit, allMesh );

% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [allMesh.x(4)     allMesh.y(4)     allMesh.z(4)];
stop  = [allMesh.x(end-3) allMesh.y(end-3) allMesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX,'--dump-statistics','--numThreads=2');

save(...
   fullfile(Sim_Path,'sim_data.mat'),...
   '-mat7-binary',...
   'f_center',...
   'f_cutoff',...
   'maxTimeSteps',...
   'portZ0',...
   'unit',...
   'portObjs');
