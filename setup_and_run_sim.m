clc
clear
close all

%
% Simulation path
%%%%%%%%%%%%%%%%%%

Sim_Path = '.\Single_Layer_PIFA_1p6mm';
dxfPath = '.\Example_DXFs';

% Library folders
%%%%%%%%%%%%%%%%%%

modelerPath = 'C:\Repositories\EM2D_Solver\Modeler_2D';
miscPath = 'C:\Repositories\EM2D_Solver\misc';
openEMSpath = 'C:\Opt\openEMS\matlab';

% Model Preparation
%%%%%%%%%%%%%%%%%%%%

% Remove points closer than this resolution (in system units)
minPointRes = 0.1;

% Resolution (given in system unit) that is considered to be the same
% point\mesh line
samePointRes = 1e-6;

% For now, manually define stackup

%layers = [  struct('name','SM_Bottom','type','prepreg','material','SR','thickness',15e-3)           ; ...
%            struct('name','B_Cu','type','conductor-bottom','material','copper','thickness',30e-3)   ; ...
%            struct('name','Core','type','laminate','material','FR-4','thickness',1.2)           ; ...
%            struct('name','F_Cu','type','conductor-top','material','copper','thickness',30e-3)      ; ...
%            struct('name','SM_Front','type','prepreg','material','SR','thickness',15e-3)];

layers = [  struct('name','B_Cu','type','conductor-bottom','material','copper','thickness',35e-3)   ; ...
            struct('name','Core','type','laminate','material','FR-4','thickness',1.51)           ; ...
            struct('name','F_Cu','type','conductor-top','material','copper','thickness',35e-3)];


metals = [  struct('name','perfect','sigma',[]) ; ...
            struct('name','copper','sigma',56e6)];

materials = [  struct('name','FR-4','epsR',4.5,'tanD',0.02) ; ...
               struct('name','SR','epsR',3.3,'tanD',0.03)];

% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%

f_center = 2.5e9;
f_cutoff = 1e9;      % Simulation bandwidth is approx. f_cutoff*2

% Maximal number of time steps allowed in simulation
maxTimeSteps = 5000*100;

% Port impedance
portZ0 = 50;

unit = 1e-3; % all length in mm

viaMaterial = 'perfect';

% How much airbox to add, in wavelength terms
airbox_pad_wl = 0.25;

endSimulationEnergy_dB = -60;

% Mesh parameters
%%%%%%%%%%%%%%%%%%
minStep_wl = [1 1 1]*1/150;

% Maximal step, in wavelength terms, inside model_data
maxStep_wl = [1/20 1/20 1/20];

% Maximal step, in wavelength terms, for airbox
airboxMaxStep_wl = [1/20 1/20 1/20];

% Minimal number of steps to consider between each "hard edge" pair
minEdge2EdgeCount_abs = [2 2 2];

% Forces a global minimum for cells sizes. Dimensions with 0 will be skipped
%minCellSize_abs = [0.1 0.1 0];
minCellSize_abs = [0 0 0];

% Cell size ratio
gradRatio = 1.5;


% Store all parameters
%%%%%%%%%%%%%%%%%%%%%%%%

save('simParams.mat','-mat7-binary');

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%% Run Simulation %%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% None of these stages has to be run multiple times. All data is saved
KiCAD_dxf_2_mat;

KiCAD_mat_define_port_faces;

mat_2_openEMS_setup_sim;

openEMS_init_mesh;

openEMS_airbox_and_sim;

openEMS_results_postProcess;
