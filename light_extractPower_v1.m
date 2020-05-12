%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% List files and retrieve layout`
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep 'CIfe_*.mat']);

%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    
    %%% extract Power
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 0.5:0.1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
    cfg.toi          = [-4*60:0.5:0];                         % time
    cfg.keeptrials  = 'yes';
    TFRhann = ft_freqanalysis(cfg, data);
    save([data_path filesep 'TF' File_Name],'TFRhann','cfg');

%     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
end
