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
List_Subj=dir([data_path filesep 'CIfre_ft_*.mat']);

%% Loop across participants to extract power
reference_method='laplacian';
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    if strcmp(reference_method,'mastoids')
        cfg=[];
        cfg.reref      = 'yes';
        cfg.refchannel = {'TP9','TP10'};
        data = ft_preprocessing(cfg,data);
        
    elseif strcmp(reference_method,'laplacian')
        cfg=[];
        cfg.method       = 'spline';
        elec = ft_read_sens('cain_elecloc_32ch.sfp');
        elec.label{match_str(elec.label,'Tp10')}='TP10';
        cfg.elec         = elec;
        cfg.trials       = 'all';
        cfg.order        = 4;
        cfg.degree       = 14;
        [data3] = ft_scalpcurrentdensity(cfg, data);
    end
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
    if strcmp(reference_method,'mastoids')
        save([data_path filesep 'mastTF' File_Name],'TFRhann','cfg');
    elseif strcmp(reference_method,'laplacian')
        save([data_path filesep 'laplTF' File_Name],'TFRhann','cfg');
    else
        save([data_path filesep 'TF' File_Name],'TFRhann','cfg');
    end
    %     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
end
