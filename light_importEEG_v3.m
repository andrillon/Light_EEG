%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m
run light_badChannels.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep '*/*.eeg']);

%% Loop across participants to extract power
duration_epoch=[-4.5 0.5]; % in minutes
nc=0;
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    nc=nc+1;
    fprintf('... processing %s (%g/%g)\n',File_Name,nS,length(List_Subj))
    hdr=ft_read_header([File_Path filesep File_Name]);
    
    
    
    %%% Define epochs
    cfg=[];
    cfg.trialfun             = 'lightEEG_trialfun';
    cfg.dataset             = [File_Path filesep File_Name];
    cfg.trialdef.eventtype  = '';
    cfg.trialdef.eventvalue = {'FG1','CT1','CT2','CT3','CT4','fg1','ct1','ct2','ct3','ct4'};
    cfg.trialdef.prestim    = 4.5*60;
    cfg.trialdef.poststim   = 0.5*60;
    cfg = ft_definetrial(cfg);
    data                   = ft_preprocessing(cfg); % read raw data
    data.label(1:32)={'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','Fp2'};

    this_subj=find(~(cellfun(@isempty,regexpi(BadChannels(:,1),File_Name(1:end-4)))));
    thesebadch = BadChannels{this_subj,2};  % take from list
    badchannel_labels = data.label(thesebadch);
    if ~isempty(badchannel_labels)
        fprintf('... ... interpolating %g channels\n',length(badchannel_labels))
        % find neighbours
        cfg=[];
        cfg.method        = 'triangulation';
        cfg.layout        = layout;
        cfg.feedback      = 'no';
        cfg.channel = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','Fp2'};
        [neighbours] = ft_prepare_neighbours(cfg);
        
        % interpolate channels
        cfg=[];
        cfg.method         = 'weighted';
        cfg.badchannel     = badchannel_labels;
        cfg.missingchannel = [];
        cfg.neighbours     = neighbours;
        cfg.trials         = 'all';
        cfg.layout         = layout;
        cfg.channel = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','Fp2'};
        [data] = ft_channelrepair(cfg, data);
    end
    
    %%% Re-reference and high-pass filter
    cfg=[];
    cfg.reref      = 'yes';
    cfg.channel = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','Fp2'};
    cfg.refchannel = {'TP9','TP10'};
    data.label(1:32)=cfg.channel;
    
    data = ft_preprocessing(cfg,data);
    %     data.label(match_str(hdr.label,'Tp10'))={'TP10'};
    save([data_path filesep 're_ft_' File_Name(1:end-4)],'data','cfg');
end


