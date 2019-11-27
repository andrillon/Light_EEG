%% init
clear all
close all

% data_path='/Users/tand0009/Data/Cain_Light/';
data_path='/Volumes/tLab_BackUp1/Monash/Cain_Light/';

path_eeglab='/Users/tand0009/Work/local/eeglab';
path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';

addpath(genpath(path_LSCPtools));
addpath(path_fieldtrip);
ft_defaults;

%% Parameters Power Extraction
w_window  = 4; % in seconds
w_overlap = (w_window/2);

%%
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep '**/*.eeg']);
nc=0;

%%
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    if strcmp(File_Name,'DLT001.eeg') || strcmp(File_Name,'DLT004.eeg') || strcmp(File_Name,'DLT016.eeg') || strcmp(File_Name,'DLT024.eeg') || strcmp(File_Name,'DLT026.eeg')
        continue;
    end
    File_Path = List_Subj(nS).folder;
    nc=nc+1;
    fprintf('... processing %s (%g/%g)\n',File_Name,nS,length(List_Subj))
    
    hdr     = ft_read_header([File_Path filesep File_Name]);
    events  = ft_read_event([File_Path filesep File_Name]);
    data    = ft_read_data([File_Path filesep File_Name]);
    
    if SubInfo.Condition(match_str(SubInfo.PT_Code,File_Name(1:end-4)))=='D'
        condition(nc) = 0;
    elseif SubInfo.Condition(match_str(SubInfo.PT_Code,File_Name(1:end-4)))=='E'
        condition(nc) = 1;
    end
    
    %%% find events
    % clean events
    for ne=1:length(events)
        if isempty(events(ne).value)
            events(ne).value='null';
        end
    end
    events_times=[events.sample];
    
    % baseline
    baseline_beg = events_times(find(~(cellfun(@isempty,regexpi({events.value},'baseline')))));
    baseline_end = events_times(find(~(cellfun(@isempty,regexpi({events.value},'baseline'))))+1);
    
    % find boundaries of FG blocks
    for nb=1:4
        if strcmp(File_Name,'DLT024.eeg') && nb==1
            fix_beg(nb) = events_times(find(~(cellfun(@isempty,regexp({events.value},sprintf('FG %g real',nb))))));
        else
            fix_beg(nb) = events_times(find(~(cellfun(@isempty,regexp({events.value},sprintf('FG%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('FG %g',nb))))));
        end
        fix_end(nb) = events_times(find(~(cellfun(@isempty,regexp({events.value},sprintf('CT%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('CT %g',nb))))));
    end
    
    %%% extract power spectrum
    % %     % re-reference to av. mastoids
    % %     data=data-repmat(mean(data([10 21],:),1),[size(data,1), 1]);
    % re-reference to av. mastoids
           data=data-repmat(mean(data(match_str(hdr.label,{'TP9','TP10'}),:),1),[size(data,1), 1]);
 for nb=1:4
                this_data=squeeze(data(:,fix_end(nb)+((-4*60*hdr.Fs+1):0))); % select 4 minutes before the end of the FG period
                this_data=this_data-repmat(mean(this_data,2),1,size(this_data,2));
                    
    [twa_results]=twalldetectnew_TA(this_data,hdr.Fs,0);
                
    for nE=1:63
        all_Waves=[all_Waves ; [repmat([n npr nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
            cell2mat(twa_results.channels(nE).negzx)' ...
            cell2mat(twa_results.channels(nE).poszx)' ...
            cell2mat(twa_results.channels(nE).wvend)' ...
            cell2mat(twa_results.channels(nE).maxnegpk)' ...
            cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
            cell2mat(twa_results.channels(nE).maxpospk)' ...
            cell2mat(twa_results.channels(nE).maxpospkamp)' ...
            cell2mat(twa_results.channels(nE).mxdnslp)' ...
            cell2mat(twa_results.channels(nE).mxupslp)' ...
            ]];
    end
   end
    % loop across electrodes
    EEG_channels=hdr.label;
    EEG_channels(match_str(EEG_channels,{'x_dir','y_dir','z_dir','ECG1','ECG2','EOG R','EOG L'}))=[];
    pick_channels=hdr.label; %find(ismember(hdr.label,EEG_channels));
    for nelec=1:length(pick_channels)
    
    end
end

%% figures
