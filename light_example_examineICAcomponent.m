%% initiliase - clear all variables and scripts
clear all
close all

%% set-up
run localdef_ligthEEG.m
addpath(path_fieldtrip); % Fieldtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% plot the components for visual inspection
% subj = 'Ife_ft_DLT001.mat'
% load('Ife_ft_DLT001.mat')
File_Name = 'Ife_ft_DLT001.mat'
load([data_path filesep File_Name]);

figure
cfg = [];
cfg.component = 1:size(comp.topo,1);       % specify the component(s) that should be plotted
cfg.layout    = 'cain_elecloc_32ch_layout.mat'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

cfg = [];
cfg.layout = 'cain_elecloc_32ch_layout.mat'; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
ft_databrowser(cfg, comp)

% cfg = [];
% cfg.component = [1:12]; % to be removed component(s)
% data2 = ft_rejectcomponent(cfg, comp, data);