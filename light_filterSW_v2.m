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
List_Subj=dir([data_path filesep 'SW_all_CIfIfre_*.mat']);

sw_thr=[];
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
    
    %%% parameters of SW detection
    paramSW.fixThr=[]; % if you want to use a fix threshold (eg 50) leave empty ([]) if you want to use the relative
    paramSW.prticle_Thr=90; % Choose percentile that you want to select: 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % Freq range you want to select: [1 4] or [4 10] in Hz
    paramSW.AmpCriterionIdx=4; % Criterion to select waves on: 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.art_ampl=150; % Rejection criterion
    paramSW.max_posampl=75; % Rejection criterion
    paramSW.max_Freq=7; % Rejection criterion
    
    %%% clean SW detection
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    %%% define SW threshold on baselone
    for nE=1:length(labels)
        thisE_Waves=all_Waves(all_Waves(:,2)==1 & all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave=paramSW.fixThr;
        else
            thr_Wave=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        sw_thr=[sw_thr ; [nS 1 nE CondSubj(nS)=='E' thr_Wave]];
        %         if thr_Wave>80
        %             labels{nE}
        %             pause;
        %         end
    end
    
    slow_Waves=[];
    for nE=1:length(labels)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        thr_Wave=sw_thr(sw_thr(:,1)==nS & sw_thr(:,3)==nE,5);
        if isempty(thr_Wave) || length(thr_Wave)>1
            continue;
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave,:)];
    end
    File_Name2=File_Name(bound{1}(2)+1:end);
    save([data_path filesep 'SW_90P2PbyE_depleted_' File_Name2(1:end-4)],'slow_Waves','labels','Fs','paramSW');
    
end

%%
cmap=colormap('parula'); %cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
figure; %set(gcf,'Position',[64          33        1097         952]);
Conds={'D','E','D-E'};
for nC=1:3
    subplot(1,3,nC);
    if nC==3
        ThresholdSW_All1=sw_thr(sw_thr(:, 4) == 1-1,5);
        ThresholdSW_Group1=sw_thr(sw_thr(:, 4) == 1-1,3);
        ThresholdSW1=grpstats(ThresholdSW_All1,ThresholdSW_Group1,'mean');
        
        ThresholdSW_All2=sw_thr(sw_thr(:, 4) == 2-1,5);
        ThresholdSW_Group2=sw_thr(sw_thr(:, 4) == 2-1,3);
        ThresholdSW2=grpstats(ThresholdSW_All2,ThresholdSW_Group2,'mean');
        
        ThresholdSW=ThresholdSW1-ThresholdSW2;
    else
        ThresholdSW_All=sw_thr(sw_thr(:, 4) == nC-1,5);
        ThresholdSW_Group=sw_thr(sw_thr(:, 4) == nC-1,3);
        ThresholdSW=grpstats(ThresholdSW_All,ThresholdSW_Group,'mean');
    end
    simpleTopoPlot_ft(ThresholdSW, layout,'on',[],0,1);
    title(sprintf('%s',Conds{nC}));
    colormap(cmap);
    colorbar;
    if nC==3
        caxis([-1 1]*12);
    else
        caxis([20 40]);
    end
    format_fig;
end

%%
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
figure; %set(gcf,'Position',[64          33        1097         952]);
Conds={'D','E'};
for nC=1:2
    subplot(1,2,nC);
    ThresholdSW_All=sw_thr(sw_thr(:, 4) == nC-1,5);
    ThresholdSW_Group=sw_thr(sw_thr(:, 4) == nC-1,3);
    ThresholdSW=grpstats(ThresholdSW_All,ThresholdSW_Group,'std');
    
    simpleTopoPlot_ft(ThresholdSW, layout,'labels',[],0,1);
    title(sprintf('%s',Conds{nC}));
    colormap(cmap);
    colorbar;
    %     caxis([20 40]);
    format_fig;
end