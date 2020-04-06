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

mkdir([data_path filesep 'fig_badCh']);

%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    
    %%% Filter
    cfg                = [];
    cfg.hpfilter       = 'yes';        % enable high-pass filtering
    cfg.hpfilttype     = 'but';
    cfg.hpfiltord         = 4;
    cfg.hpfreq         = 0.1;
    
    cfg.lpfilter       = 'yes';        % enable high-pass filtering
    cfg.lpfilttype     = 'but';
    cfg.lpfiltord         = 4;
    cfg.lpfreq         = 30;

    cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
    cfg.dftfreq        = [50]; % set up the frequencies for notch filtering
    data               = ft_preprocessing(cfg,data);

    
    %%%
    all_data=[];
    for nb=1:length(data.trial)
        all_data=[all_data data.trial{nb}];
    end
    figure; format_fig;
    imagesc(abs(all_data));
    caxis([0 1]*500);
    colorbar;
    this_title=File_Name;
    this_title(findstr(this_title,'_'))=' ';
    this_title(end-3:end)=[];
    title(this_title);
    set(gca,'YTick',1:length(data.label),'YTickLabel',data.label);
    print(gcf, [data_path filesep 'fig_badCh' filesep  File_Name '.png'], '-dpng');
end
