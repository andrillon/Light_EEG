%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep 'e_*.mat']);

%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    %%% Re-reference to average THOMAS ADDS DETREND
    cfg=[];
    cfg.reref      = 'yes';
    cfg.refchannel = 'all';
    
    cfg.demean        = 'yes';
    
    cfg.hpfilter       = 'yes';        % enable high-pass filtering
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord         = 4;
    cfg.hpfreq         = 0.1;
    
    cfg.lpfilter       = 'yes';        % enable high-pass filtering
    cfg.lpfilttype     = 'but';
    cfg.lpfiltord         = 4;
    cfg.lpfreq         = 40;

    cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq        = [50]; % set up the frequencies for notch filtering
    
    data = ft_preprocessing(cfg,data);
    
    %%% run ICA
    cfg        = [];
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
    comp = ft_componentanalysis(cfg, data);
    save([data_path filesep 'If' File_Name],'data','comp');
end
