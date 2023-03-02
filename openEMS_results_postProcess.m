% openEMS_results_postProcess.m
clc
clear
close all

load('simParams.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize openEMS environment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(Sim_Path,'model_data.mat'));
load(fullfile(Sim_Path,'model_ports.mat'));
load(fullfile(Sim_Path,'sim_data.mat'));

pkg load geometry
pkg load miscellaneous

addpath(modelerPath);
addpath(miscPath);
addpath(openEMSpath);

physical_constants;

port = portObjs{1};
%% postprocessing & do the plots
freq = linspace( f_center - f_cutoff, f_center + f_cutoff, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = real(0.5 * port.uf.tot .* conj( port.if.tot )); % antenna feed power

figure;
plot(freq/1e9,20*log10(abs(s11(:))),'-b','linewidth',3);
xlabel('f [GHz]','fontsize',16);
ylabel('S_{11} [dB]','fontsize',16);
set(gca,'fontsize',14);

grid on;


