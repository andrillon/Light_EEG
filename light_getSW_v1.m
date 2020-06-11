%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox
addpath(genpath(fooof_path))

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep 'CIfIfe_*.mat']);

f_range = [2, 40];
settings = struct();  % Use defaults
av_fooof_bg=[];
av_fooof_alpha=[];
%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    %     if strcmp(File_Name,'TFe_ft_DLT018.mat')
    %         continue;
    %     end
    fprintf('... processing %s (%g/%g)',File_Name,nS,length(List_Subj))
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
    all_Waves=[];
    for nB=1:5
        this_data=data.trial{nB};
        this_data=this_data-repmat(mean(this_data(match_str(data.label,{'TP9','TP10'}),:),1),[size(this_data,1), 1]);
        [twa_results]=twalldetectnew_TA_v2(this_data,data.fsample,0);
        
        for nE=1:size(this_data,1)
            all_Waves=[all_Waves ; [repmat([nS nB nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                cell2mat(twa_results.channels(nE).negzx)' ...
                cell2mat(twa_results.channels(nE).poszx)' ...
                cell2mat(twa_results.channels(nE).wvend)' ...
                cell2mat(twa_results.channels(nE).maxnegpk)' ...
                cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                cell2mat(twa_results.channels(nE).maxpospk)' ...
                cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                cell2mat(twa_results.channels(nE).mxdnslp)' ...
                cell2mat(twa_results.channels(nE).mxupslp)' ...
                cell2mat(twa_results.channels(nE).maxampwn)' ...
                cell2mat(twa_results.channels(nE).minampwn)' ...
                ]];
        end        
    end
    labels=data.label;
    Fs=data.fsample;
    save([data_path filesep 'SW_all_' File_Name],'all_Waves','labels','Fs');

end

