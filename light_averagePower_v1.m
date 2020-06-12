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
    av_logPower(nS,:,:,:)=log(squeeze(mean(TFRhann.powspctrm,4)));
    %     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
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

thisChLabel='Cz';
freqs=TFRhann.freq;

figure;
for nCond=1:2
    subplot(1,2,nCond)
    hold on;
    hp=[];
    for nBl=1:5
        temp_pow=squeeze(av_logPower(CondSubj==Conds{nCond},nBl,match_str(TFRhann.label,thisChLabel),:));
        [~,hp(nBl)]=simpleTplot(freqs,temp_pow,0,Colors{nCond}(nBl,:),0,'-',0.5,1,0,[],2);
    end
    format_fig;
    xlabel('Freq (Hz)');
    ylabel('log(Power)');
    legend(hp,{'B0','B1','B2','B3','B4'});
    title(sprintf('%s - %s',Conds{nCond},thisChLabel));
    ylim([-4 6])
end


%%
thisChLabel='Pz';
figure;
for nCond=1:2
    hold on;
    hp=[];
    for nBl=1:5
        temp_pow=squeeze(av_logPower(CondSubj==Conds{nCond},nBl,match_str(TFRhann.label,thisChLabel),:));
        simpleTplot(freqs,temp_pow,0,Colors{nCond}(nBl,:),0,'-',0.5,1,0,0,2);
    end
    format_fig;
    xlabel('Freq (Hz)');
    ylabel('log(Power)');
%     legend(hp,{'B0','B1','B2','B3','B4'});
    title(sprintf('%s - %s',Conds{nCond},thisChLabel));
    ylim([-4 2])
   xlim([1 30])
end
%%
FOI=[6 11]; % Freq Band of Interest
figure; set(gcf,'Position',[64          33        1097         952]);
for nCond=1:2
    for nB=1:5
        subplot(3,5,5*(nCond-1)+nB);
        Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{nCond},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
        
        simpleTopoPlot_ft(Pow_AVG, layout,'on',[],0,1);
        title(sprintf('%s - %s',Conds{nCond},thisChLabel));
        colorbar;
        caxis([-2 0]);
        format_fig;
    end
end

for nB=1:5
    subplot(3,5,10+nB);
    Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1))-...
        squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
    
    simpleTopoPlot_ft(Pow_AVG, layout,'on',[],0,1);
    title(sprintf('%s - %s','E vs D',thisChLabel));
    colorbar;
    caxis([-1 1]*1);
    format_fig;
end