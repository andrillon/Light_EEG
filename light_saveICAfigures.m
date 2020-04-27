%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox
mkdir([data_path filesep 'fig_Comp']);

%% specify data file

List_Subj=dir([data_path filesep 'Ife_ft_*.mat']);

%% plot the components for visual inspection
for nS=1:length(List_Subj)
    
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    figure
    cfg = [];
    cfg.component = 1:size(comp.topo,1);      % specify the component(s) that should be plotted
    cfg.layout    = 'cain_elecloc_32ch_layout.mat'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    print(gcf, [data_path filesep 'fig_Comp' filesep  File_Name '.png'], '-dpng');
end
