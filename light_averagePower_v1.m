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
List_Subj=dir([data_path filesep 'TFe_*.mat']);

%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    fprintf('... processing %s (%g/%g)\n',File_Name,nS,length(List_Subj))
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    if nS==1
        av_logPower=nan([length(List_Subj) size(TFRhann.powspctrm,1) size(TFRhann.powspctrm,2) size(TFRhann.powspctrm,3)]);
    end
    av_logPower(nS,:,:,:)=log(squeeze(mean(TFRhann.powspctrm,4)));
%     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
end
