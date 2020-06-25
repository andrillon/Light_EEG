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
List_Subj=dir([data_path filesep 'TFCIfIfe_*.mat']);

%% Loop across participants to extract power
tlCondE=[];
tlCondD=[];
design=cell(1,2);
TFRhann_all=[];
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    %     if strcmp(File_Name,'TFe_ft_DLT018.mat')
    %         continue;
    %     end
    fprintf('... processing %s (%g/%g)',File_Name,nS,length(List_Subj))
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    if nS==1
        av_logPower=nan([length(List_Subj) size(TFRhann.powspctrm,1) size(TFRhann.powspctrm,2) size(TFRhann.powspctrm,3)]);
    end
    av_logPower(nS,:,:,:)=(squeeze(mean(log(TFRhann.powspctrm),4)));
    %     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
%     TFRhann_all{nS}=TFRhann;
%     TFRhann_all{nS}.powspctrm=squeeze(mean(TFRhann.powspctrm(2:5,:,:,:),1));
%     TFRhann_all{nS}.dimord='chan_freq_time';
%     
    if CondSubj(nS)=='E'
        design{1}=[design{1} [ones(1,4) ; 1:4]];
    elseif CondSubj(nS)=='D'
        design{2}=[design{2} [2*ones(1,4) ; 1:4]];
    end
    if CondSubj(nS)=='E'
        if isempty(tlCondE)
            tlCondE=TFRhann;
            tlCondE.dimord='subj_chan_freq';
%             temp=log(squeeze(mean(TFRhann.powspctrm,4)));
            tlCondE.powspctrm=[];
            tlCondE.powspctrm=squeeze(log(mean(TFRhann.powspctrm(2:5,:,:,:),4)));
            nc1=1;
        else
%             temp=log(squeeze(mean(TFRhann.powspctrm,4)));
            nc1=nc1+1;
            tlCondE.powspctrm=cat(1,tlCondE.powspctrm,squeeze(log(mean(TFRhann.powspctrm(2:5,:,:,:),4))));
        end
    elseif CondSubj(nS)=='D'
        if isempty(tlCondD)
            tlCondD=TFRhann;
            tlCondD.dimord='subj_chan_freq';
%             temp=log(squeeze(mean(TFRhann.powspctrm,4)));
            tlCondD.powspctrm=[];
            tlCondD.powspctrm=squeeze(log(mean(TFRhann.powspctrm(2:5,:,:,:),4)));
            nc2=1;
        else
%             temp=log(squeeze(mean(TFRhann.powspctrm,4)));
            nc2=nc2+1;
            tlCondD.powspctrm=cat(1,tlCondD.powspctrm,squeeze(log(mean(TFRhann.powspctrm(2:5,:,:,:),4))));
        end
    end
end

%%
Conds={'E','D'};
Colors={[100,100,100
    127,205,187
    65,182,196
    44,127,184
    37,52,148]/256,...
    [100,100,100
    254,178,76
    253,141,60
    240,59,32
    189,0,38]/256};
tlCondE=rmfield(tlCondE,'time');
tlCondD=rmfield(tlCondD,'time');


%%
cfg = [];
cfg.channel          = 'all';
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
% cfg.resampling='bootstrap';
cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 	500;
% prepare_neighbours determines what sensors may form clusters
cfglay=[];
cfglay.layout=layout;
layout = ft_prepare_layout(cfglay);
cfg_neighb=[];
cfg_neighb.method    = 'triangulation';
cfg_neighb.layout=cfglay.layout;
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TFRhann);

cfg.design           = [design{1}(1,:) design{2}(1,:)];%[ones(1,size(tlCondE.powspctrm,1)) 2*ones(1,size(tlCondD.powspctrm,1))];
cfg.ivar             = 1;

% [stat] = ft_freqstatistics(cfg, TFRhann_all{1}, TFRhann_all{2}, TFRhann_all{3}, TFRhann_all{4}, TFRhann_all{5}, TFRhann_all{6}, TFRhann_all{7}, TFRhann_all{8}, TFRhann_all{9}, TFRhann_all{10},...
%      TFRhann_all{11}, TFRhann_all{12}, TFRhann_all{13}, TFRhann_all{14}, TFRhann_all{15}, TFRhann_all{16}, TFRhann_all{17}, TFRhann_all{18}, TFRhann_all{19}, TFRhann_all{20},...
%      TFRhann_all{21}, TFRhann_all{22}, TFRhann_all{23}, TFRhann_all{24}, TFRhann_all{25}, TFRhann_all{26}, TFRhann_all{27}, TFRhann_all{28}, TFRhann_all{29}, TFRhann_all{30},...
%      TFRhann_all{31}, TFRhann_all{32}, TFRhann_all{33}, TFRhann_all{34}, TFRhann_all{35}, TFRhann_all{36}, TFRhann_all{37}, TFRhann_all{38}, TFRhann_all{39});
tlCondD2=tlCondD;
tlCondE2=tlCondE;
% tlCondD2.powspctrm(:,[2 3 30 24 7 29 1 32],50:55)=tlCondD2.powspctrm(:,[2 3 30 24 7 29 1 32],50:55)+50;
tlCondD2.powspctrm(:,[2 3],50:55)=tlCondD2.powspctrm(:,[2 3],50:55)+50;
% tlCondE2.powspctrm=permute(tlCondE2.powspctrm,[1 3 2]);
% tlCondD2.powspctrm=permute(tlCondD2.powspctrm,[1 3 2]);
% tlCondD2.dimord='subj_freq_chan';
% tlCondE2.dimord='subj_freq_chan';

[stat] = ft_freqstatistics(cfg, tlCondE,tlCondD);
%%
cfg = [];
% tlCondE=rmfield(tlCondE,'time');
% tlCondD=rmfield(tlCondD,'time');
freqE_cmb = ft_freqdescriptives(cfg, tlCondE);
freqD_cmb  = ft_freqdescriptives(cfg, tlCondD);


stat.raweffect = freqD_cmb.powspctrm - freqE_cmb.powspctrm;
cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'raweffect';
% cfg.zlim   = [-1e-27 1e-27];
cfg.layout = layout;
ft_clusterplot(cfg, stat);