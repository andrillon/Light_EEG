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
List_Subj=dir([data_path filesep '**/*.eeg']);

%% Loop across participants to extract power
duration_epoch=[-4.5 0.5]; % in minutes
nc=0;
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    File_Path = List_Subj(nS).folder;
    nc=nc+1;
    fprintf('... processing %s (%g/%g)\n',File_Name,nS,length(List_Subj))
    
    hdr     = ft_read_header([File_Path filesep File_Name]);
    events  = ft_read_event([File_Path filesep File_Name]);
    data    = ft_read_data([File_Path filesep File_Name]);
    
    % clean events
    if strcmp(File_Name,'DLT001.eeg')
        events(4).value='Baseline';
        events(5).value='FG1'; events(5).type='Comment'; events(5).sample=583841;
        events(6).value='CT1'; events(6).sample=793500;
        events(7).value='FG2'; events(7).sample=883500;
        events(8).value='CT2'; events(8).sample=1033500;
        events(9).value='FG3'; events(9).sample=1123500;
        events(10).value='CT3'; events(10).sample=1273500;
        events(11).value='FG4'; events(11).sample=1363500;
        events(12).value='CT4'; events(12).sample=1513500;
    elseif strcmp(File_Name,'DLT004.eeg')
        events([5 6 9])=[];
    elseif strcmp(File_Name,'DLT016.eeg')
        events(9).value='CT3';
        events(10).value='FG4';
    elseif strcmp(File_Name,'DLT024.eeg')
        events(3).value='FG1'; events(2).value='Baseline';
    elseif strcmp(File_Name,'DLT026.eeg')
        events(1)=[];
        events(find(~(cellfun(@isempty,regexp({events.value},'^LostSamples')))))=[];
    elseif strcmp(File_Name,'DLT038.eeg')
        fprintf('... SKIPPING (baseline too short)\n')
       continue;
    end
    for ne=1:length(events)
        if isempty(events(ne).value)
            events(ne).value='null';
        end
    end
    %%% extract samples
    events_times=[events.sample];
    % baseline
    baseline_beg = events_times(find(~(cellfun(@isempty,regexpi({events.value},'baseline')))));
    baseline_end = events_times(find(~(cellfun(@isempty,regexpi({events.value},'baseline'))))+1);
    baseline_duration = (baseline_end-baseline_beg)/hdr.Fs/60;
    fprintf('... ... baseline: %g min\n',baseline_duration)
    % find boundaries of FG blocks
    for nb=1:4
        fix_beg(nb) = events_times(find(~(cellfun(@isempty,regexp({events.value},sprintf('FG%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('FG %g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('fg%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('fg %g',nb))))));
        fix_end(nb) = events_times(find(~(cellfun(@isempty,regexp({events.value},sprintf('CT%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('CT %g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('ct%g',nb)))) | ~(cellfun(@isempty,regexp({events.value},sprintf('ct %g',nb))))));
        fix_duration(nb) = (fix_end(nb)-fix_beg(nb))/hdr.Fs/60;
        fprintf('... ... FG block %g: %g min\n',nb,fix_duration(nb))
    end
    
    %%% Cut the data
    epoched_data=nan(size(data,1),(duration_epoch(2)-duration_epoch(1))*60*hdr.Fs+1,5);
    % baseline power
    epoched_data(:,:,1)=data(:,baseline_end+((duration_epoch(1)*60*hdr.Fs+1):(duration_epoch(2)*60*hdr.Fs+1))); % select 4 minutes before the end of the FG period
    %     pow_fix=[];
    for nb=1:4
        epoched_data(:,:,1+nb)=data(:,fix_end(nb)+((duration_epoch(1)*60*hdr.Fs+1):(duration_epoch(2)*60*hdr.Fs+1))); % select 4 minutes before the end of the FG period
    end
    
    %%% Save the data
    save([data_path filesep 'e_' File_Name(1:end-4)],'epoched_data','hdr','events','duration_epoch');
end
