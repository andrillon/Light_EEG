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
% load('LightEEG_ICA_EyeMovements_v2.mat')
load LightEEG_ICA_V4.mat
List_Subj=dir([data_path filesep 're_*.mat']);

%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    load([data_path filesep 'If' File_Name],'comp');
    
    
    %%% Re-reference to average THOMAS ADDS DETREND
    cfg = [];
    eval(sprintf('cfg.component = %s;',BadComponents.Eyemovement{match_str(BadComponents.Participant,File_Name(8:13))})); % to be removed component(s)
    %     cfg.component=BadComponents{match_str(BadComponents(:,1),File_Name(7:12)),2};
    data = ft_rejectcomponent(cfg, comp, data);
    save([data_path filesep 'CIf' File_Name],'data');
    
end
