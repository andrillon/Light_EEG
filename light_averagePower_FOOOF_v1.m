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
List_Subj=dir([data_path filesep 'TFCIfIfe_*.mat']);

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
    
    temp_pow=squeeze(mean(TFRhann.powspctrm,4));
    for nB=1:5
        for nE=1:32
            try
            fooof_results = fooof(TFRhann.freq, squeeze(temp_pow(nB,nE,:))', f_range, settings,0);
%             if nE==17
%                 figure;
%                 fooof_plot(fooof_results, 0)
%             end
            av_fooof_bg=[av_fooof_bg ; [nS nB nE CondSubj(nS)=='E' fooof_results.background_params]];
            
            [closestvalue,index]=findclosest(fooof_results.peak_params(:,1),10);
            if closestvalue>7 && closestvalue<12
                av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' fooof_results.peak_params(index,:)]];
            else
                av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.peak_params,2))]];
            end
            catch
                       av_fooof_bg=[av_fooof_bg ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.background_params,2))]];
                   av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.peak_params,2))]];
 continue; 
            end
        end
    end
end

%%
Conds={'D','E'};
Colors={
    [100,100,100
    254,178,76
    253,141,60
    240,59,32
    189,0,38]/256,...
    [100,100,100
    127,205,187
    65,182,196
    44,127,184
    37,52,148]/256};


%% LME on background and alpha peak
table_bg=array2table(av_fooof_bg,'VariableNames',{'SubID','BlockN','ElecN','CondF','Offset','Slope'});
table_bg.SubID=categorical(table_bg.SubID);
table_bg.Cond=table_bg.CondF;
table_bg.Cond=categorical(table_bg.Cond);
table_bg.Cond(table_bg.Cond=='1')='E';
table_bg.Cond(table_bg.Cond=='0')='D';
table_bg.Cond=removecats(table_bg.Cond);

table_bg.Elec=table_bg.ElecN;
table_bg.Elec=categorical(table_bg.Elec);
for nEl=1:length(TFRhann.label)
    table_bg.Elec(table_bg.Elec==num2str(nEl))=TFRhann.label{nEl};
end
table_bg.Elec=removecats(table_bg.Elec);

mdl_0=fitlme(table_bg,'Slope~1+BlockN+(1|SubID)');
mdl_1=fitlme(table_bg,'Slope~1+BlockN+Elec+(1|SubID)');
mdl_2=fitlme(table_bg,'Slope~1+BlockN+Cond+(1|SubID)');
mdl_3=fitlme(table_bg,'Slope~1+BlockN*Cond+(1|SubID)');
% compare(mdl_0,mdl_2);


mdl_0b=fitlme(table_bg,'Offset~1+BlockN+(1|SubID)');
mdl_1b=fitlme(table_bg,'Offset~1+BlockN+Elec+(1|SubID)');
mdl_2b=fitlme(table_bg,'Offset~1+BlockN+Cond+(1|SubID)');
mdl_3b=fitlme(table_bg,'Offset~1+BlockN*Cond+(1|SubID)');
% compare(mdl_0b,mdl_2b)


table_alpha=array2table(av_fooof_alpha,'VariableNames',{'SubID','BlockN','ElecN','CondF','Freq','Amp','BW'});
table_alpha.SubID=categorical(table_alpha.SubID);
table_alpha.Cond=table_alpha.CondF;
table_alpha.Cond=categorical(table_alpha.Cond);
table_alpha.Cond(table_alpha.Cond=='1')='E';
table_alpha.Cond(table_alpha.Cond=='0')='D';
table_alpha.Cond=removecats(table_alpha.Cond);

table_alpha.Elec=table_alpha.ElecN;
table_alpha.Elec=categorical(table_alpha.Elec);
for nEl=1:length(TFRhann.label)
    table_alpha.Elec(table_alpha.Elec==num2str(nEl))=TFRhann.label{nEl};
end
table_alpha.Elec=removecats(table_alpha.Elec);

mdla_0=fitlme(table_alpha,'Freq~1+(1|SubID)');
mdla_1=fitlme(table_alpha,'Freq~1+Elec+(1|SubID)');
mdla_2=fitlme(table_alpha,'Freq~1+BlockN+Elec+(1|SubID)');
mdla_3=fitlme(table_alpha,'Freq~1+BlockN+Cond+(1|SubID)');
mdla_4=fitlme(table_alpha,'Freq~1+BlockN*Cond+(1|SubID)');


mdlb_0=fitlme(table_alpha,'Amp~1+(1|SubID)');
mdlb_1=fitlme(table_alpha,'Amp~1+Elec+(1|SubID)');
mdlb_2=fitlme(table_alpha,'Amp~1+BlockN+Elec+(1|SubID)');
mdlb_3=fitlme(table_alpha,'Amp~1+BlockN+Cond+(1|SubID)');
mdlb_4=fitlme(table_alpha,'Amp~1+BlockN*Cond+(1|SubID)');
%%
freqs=cfg.foi;
FOI=[9 11]; % Freq Band of Interest
figure; set(gcf,'Position',[64          33        1097         952]);
for nCond=1:2
    for nB=1:5
        subplot(2,5,5*(nCond-1)+nB);
        temp_topo=[];
        for nEl=1:32
        temp_topo(nEl)=nanmean(table_alpha.Amp(table_alpha.Cond==Conds{nCond} & table_alpha.BlockN==nB & table_alpha.ElecN==nEl));
        end
        simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%         title(sprintf('%s - %s',Conds{nCond},thisChLabel));
        colorbar;
        caxis([.2 1]);
        format_fig;
    end
end
