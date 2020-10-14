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
List_Subj=dir([data_path filesep 'TFCIfre_ft_*.mat']);

f_range = [2, 40];
settings = struct();  % Use defaults
av_fooof_bg=[];
av_fooof_alpha=[];
av_fooof_theta=[];
av_fooof_peaks=[];
av_fooof_maxalphatheta=[];
%% Loop across participants to extract power
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    if strcmp(File_Name,'TFCIfre_ft_DLT037.mat')
        continue;
    end
    fprintf('... processing %s (%g/%g)',File_Name,nS,length(List_Subj))
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    %     if size(TFRhann.powspctrm,1)==4
    %         continue;
    %     end
    %     if nS==1
    %         av_logPower=nan([length(List_Subj) size(TFRhann.powspctrm,1) size(TFRhann.powspctrm,2) size(TFRhann.powspctrm,3)]);
    %     end
    av_logPower(nS,:,:,:)=log(squeeze(mean(TFRhann.powspctrm,4)));
    
    %     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
    
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
    temp_pow=squeeze(mean(TFRhann.powspctrm(:,:,:,TFRhann.time<-30),4));
    for nB=1:5
        for nE=1:32
            try
                fooof_results = fooof(TFRhann.freq, squeeze(temp_pow(nB,nE,:))', f_range, settings,0);
                %             if nE==17
                %                 figure;
                %                 fooof_plot(fooof_results, 0)
                %             end
                av_fooof_bg=[av_fooof_bg ; [nS nB nE CondSubj(nS)=='E' fooof_results.background_params]];
                av_fooof_peaks=[av_fooof_peaks ; [repmat([nS nB nE CondSubj(nS)=='E'],size(fooof_results.peak_params,1),1) fooof_results.peak_params]];
                
                sub_freqpeaks=fooof_results.peak_params;
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,1)>9.5 & sub_freqpeaks(:,1)<14,:);
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,2)==max(sub_freqpeaks(:,2)),:);
                if ~isempty(sub_freqpeaks)
                    av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' sub_freqpeaks]];
                else
                    av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' nan(1,3)]];
                end
                %                 [closestvalue,index]=findclosest(fooof_results.peak_params(:,1),10);
                %                 if closestvalue>9 && closestvalue<12
                %                     av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' fooof_results.peak_params(index,:)]];
                %                 else
                %                     av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.peak_params,2))]];
                %                 end
                
                sub_freqpeaks=fooof_results.peak_params;
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,1)>5 & sub_freqpeaks(:,1)<9.5,:);
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,2)==max(sub_freqpeaks(:,2)),:);
                if ~isempty(sub_freqpeaks)
                    av_fooof_theta=[av_fooof_theta ; [nS nB nE CondSubj(nS)=='E' sub_freqpeaks]];
                else
                    av_fooof_theta=[av_fooof_theta ; [nS nB nE CondSubj(nS)=='E' nan(1,3)]];
                end
                
                sub_freqpeaks=fooof_results.peak_params;
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,1)>7 & sub_freqpeaks(:,1)<12,:);
                sub_freqpeaks=sub_freqpeaks(sub_freqpeaks(:,2)==max(sub_freqpeaks(:,2)),:);
                if ~isempty(sub_freqpeaks)
                    av_fooof_maxalphatheta=[av_fooof_maxalphatheta ; [nS nB nE CondSubj(nS)=='E' sub_freqpeaks]];
                else
                    av_fooof_maxalphatheta=[av_fooof_maxalphatheta ; [nS nB nE CondSubj(nS)=='E' nan(1,3)]];
                end
            catch
                av_fooof_bg=[av_fooof_bg ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.background_params,2))]];
                av_fooof_alpha=[av_fooof_alpha ; [nS nB nE CondSubj(nS)=='E' nan(1,size(fooof_results.peak_params,2))]];
                av_fooof_maxalphatheta=[av_fooof_maxalphatheta ; [nS nB nE CondSubj(nS)=='E' nan(1,3)]];
                figure;
                plot(TFRhann.freq, squeeze(temp_pow(nB,nE,:)))
                title(File_Name)
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
for nEl=1:length(layout.label)-5
    table_bg.Elec(table_bg.Elec==num2str(nEl))=layout.label{nEl};
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

figure; set(gcf,'Position',[62   198   548   787]);
subplot(2,1,1); hold on;
for nBl=1:5
    for nCond=1:2
        temp1=table_bg.Slope(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond});
        temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond}); temp1b=removecats(temp1b);
        mean1=grpstats(temp1,temp1b);
            simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
%         line([1 1]*nBl+(2*nCond-3)*0.2,[-1 1].*sem(mean1)+mean(mean1),'Color',Colors{nCond}(nBl,:),'LineWidth',2);
%         line([Pos-0.1*widthBar Pos+0.1*widthBar],[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError)
%         line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError)
%         
%         scatter(nBl+(2*nCond-3)*0.2,mean(mean1),'MarkerEdgeColor',Colors{nCond}(nBl,:),'MarkerFaceColor',Colors{nCond}(nBl,:),'LineWidth',2);
    end
end
title('Background')
xlabel('Block')
ylabel('Slope')
format_fig;
ylim([0.7 1.1])
xlim([0.5 5.5])

subplot(2,1,2);
for nBl=1:5
    for nCond=1:2
        temp1=table_bg.Offset(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond});
        temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond}); temp1b=removecats(temp1b);
        mean1=grpstats(temp1,temp1b);
%         simpleBarPlot(nBl+(2*nCond-3)*0.2,mean1,Colors{nCond}(nBl,:),0.35,'k',[],3);
            simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
    end
end
title('Background')
format_fig;
xlabel('Block')
ylabel('Offset')
ylim([-0.5 -0])
xlim([0.5 5.5])

%%
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nBl=1:5
        subplot(3,5,5*(nCond-1)+nBl);
        topo_toplot=[];
        for nE=1:length(layout.label)-5
            temp1=table_bg.Slope(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1=grpstats(temp1,temp1b);
            topo_toplot(nE)=nanmean(mean1);
        end
        simpleTopoPlot_ft(topo_toplot', layout,'on',[],0,1);
        title(sprintf('B%g: %s',nBl,Conds{nCond}));
        colorbar;
        caxis([.5 1.5]);
        format_fig;
    end
end

for nBl=1:5
    subplot(3,5,10+nBl);
    mean1=[];
    mean2=[];
        for nE=1:length(layout.label)-5
            temp1=table_bg.Slope(table_bg.BlockN==nBl & table_bg.Cond==Conds{1} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{1} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1(:,nE)=grpstats(temp1,temp1b);
            
            temp1=table_bg.Slope(table_bg.BlockN==nBl & table_bg.Cond==Conds{2} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{2} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean2(:,nE)=grpstats(temp1,temp1b);
        end
    [h,p,~,stats]=ttest2(mean1,...
        mean2);
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nBl,'D - E'));
    colorbar;
    caxis([-1 1]*3.5);
    format_fig;
end

%%
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nBl=1:5
        subplot(3,5,5*(nCond-1)+nBl);
        topo_toplot=[];
        for nE=1:length(layout.label)-5
            temp1=table_bg.Offset(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{nCond} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1=grpstats(temp1,temp1b);
            topo_toplot(nE)=nanmean(mean1);
        end
        simpleTopoPlot_ft(topo_toplot', layout,'on',[],0,1);
        title(sprintf('B%g: %s',nBl,Conds{nCond}));
        colorbar;
        caxis([-0.8 0]);
        format_fig;
    end
end

for nBl=1:5
    subplot(3,5,10+nBl);
    mean1=[];
    mean2=[];
        for nE=1:length(layout.label)-5
            temp1=table_bg.Offset(table_bg.BlockN==nBl & table_bg.Cond==Conds{1} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{1} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1(:,nE)=grpstats(temp1,temp1b);
            
            temp1=table_bg.Offset(table_bg.BlockN==nBl & table_bg.Cond==Conds{2} & table_bg.Elec==layout.label{nE});
            temp1b=table_bg.SubID(table_bg.BlockN==nBl & table_bg.Cond==Conds{2} & table_bg.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean2(:,nE)=grpstats(temp1,temp1b);
        end
    [h,p,~,stats]=ttest2(mean1,...
        mean2);
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nBl,'D - E'));
    colorbar;
    caxis([-1 1]*3.5);
    format_fig;
end

%%
table_alpha=array2table(av_fooof_alpha,'VariableNames',{'SubID','BlockN','ElecN','CondF','Freq','Amp','BW'});
table_alpha.SubID=categorical(table_alpha.SubID);
table_alpha.Cond=table_alpha.CondF;
table_alpha.Cond=categorical(table_alpha.Cond);
table_alpha.Cond(table_alpha.Cond=='1')='E';
table_alpha.Cond(table_alpha.Cond=='0')='D';
table_alpha.Cond=removecats(table_alpha.Cond);
table_alpha.Elec=table_alpha.ElecN;
table_alpha.Elec=categorical(table_alpha.Elec);
for nEl=1:length(layout.label)-5
    table_alpha.Elec(table_alpha.Elec==num2str(nEl))=layout.label{nEl};
end
table_alpha.Elec=removecats(table_alpha.Elec);

table_theta=array2table(av_fooof_theta,'VariableNames',{'SubID','BlockN','ElecN','CondF','Freq','Amp','BW'});
table_theta.SubID=categorical(table_theta.SubID);
table_theta.Cond=table_theta.CondF;
table_theta.Cond=categorical(table_theta.Cond);
table_theta.Cond(table_theta.Cond=='1')='E';
table_theta.Cond(table_theta.Cond=='0')='D';
table_theta.Cond=removecats(table_theta.Cond);
table_theta.Elec=table_theta.ElecN;
table_theta.Elec=categorical(table_theta.Elec);
for nEl=1:length(layout.label)-5
    table_theta.Elec(table_theta.Elec==num2str(nEl))=layout.label{nEl};
end
table_theta.Elec=removecats(table_theta.Elec);

mdla_0=fitlme(table_alpha,'Amp~1+(1|SubID)');
mdla_1=fitlme(table_alpha,'Amp~1+BlockN+(1|SubID)');
mdla_2=fitlme(table_alpha,'Amp~1+BlockN+Elec+(1|SubID)');
mdla_3=fitlme(table_alpha,'Amp~1+BlockN+Cond+(1|SubID)');
mdla_4=fitlme(table_alpha,'Amp~1+BlockN*Cond+(1|SubID)');


mdlb_0=fitlme(table_theta,'Amp~1+(1|SubID)');
mdlb_1=fitlme(table_theta,'Amp~1+BlockN+(1|SubID)');
mdlb_2=fitlme(table_theta,'Amp~1+BlockN+Elec+(1|SubID)');
mdlb_3=fitlme(table_theta,'Amp~1+BlockN+Cond+(1|SubID)');
mdlb_4=fitlme(table_theta,'Amp~1+BlockN*Cond+(1|SubID)');

%%

figure; set(gcf,'Position',[62   198   548   787]);
subplot(2,1,1); hold on;
this_table=table_alpha;
for nBl=1:5
    for nCond=1:2
        temp1=this_table.Amp(this_table.BlockN==nBl & this_table.Cond==Conds{nCond});
        temp1b=this_table.SubID(this_table.BlockN==nBl & this_table.Cond==Conds{nCond}); temp1b=removecats(temp1b);
        mean1=grpstats(temp1,temp1b);
            simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
    end
end
title('Peak Alpha Ampl')
format_fig;
xlabel('Block')
ylabel('Amp')

subplot(2,1,2);
this_table=table_theta;
for nBl=1:5
    for nCond=1:2
        temp1=this_table.Amp(this_table.BlockN==nBl & this_table.Cond==Conds{nCond});
        temp1b=this_table.SubID(this_table.BlockN==nBl & this_table.Cond==Conds{nCond}); temp1b=removecats(temp1b);
        mean1=grpstats(temp1,temp1b);
            simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
    end
end
title('Peak Theta Ampl')
format_fig;
xlabel('Block')
ylabel('Amp')
%%
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nBl=1:5
        subplot(3,5,5*(nCond-1)+nBl);
        topo_toplot=[];
        for nE=1:length(layout.label)-5
            temp1=table_theta.Amp(table_theta.BlockN==nBl & table_theta.Cond==Conds{nCond} & table_theta.Elec==layout.label{nE});
            temp1b=table_theta.SubID(table_theta.BlockN==nBl & table_theta.Cond==Conds{nCond} & table_theta.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1=grpstats(temp1,temp1b);
            topo_toplot(nE)=nanmean(mean1);
        end
        simpleTopoPlot_ft(topo_toplot', layout,'on',[],0,1);
        title(sprintf('B%g: %s',nBl,Conds{nCond}));
        colorbar;
        caxis([.4 .8]);
        format_fig;
    end
end

for nBl=1:5
    subplot(3,5,10+nBl);
    mean1=[];
    mean2=[];
        for nE=1:length(layout.label)-5
            temp1=table_theta.Amp(table_theta.BlockN==nBl & table_theta.Cond==Conds{1} & table_theta.Elec==layout.label{nE});
            temp1b=table_theta.SubID(table_theta.BlockN==nBl & table_theta.Cond==Conds{1} & table_theta.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1(:,nE)=grpstats(temp1,temp1b);
            
            temp1=table_theta.Amp(table_theta.BlockN==nBl & table_theta.Cond==Conds{2} & table_theta.Elec==layout.label{nE});
            temp1b=table_theta.SubID(table_theta.BlockN==nBl & table_theta.Cond==Conds{2} & table_theta.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean2(:,nE)=grpstats(temp1,temp1b);
        end
    [h,p,~,stats]=ttest2(mean1,...
        mean2);
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nBl,'D - E'));
    colorbar;
    caxis([-1 1]*2.5);
    format_fig;
end


%%
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nBl=1:5
        subplot(3,5,5*(nCond-1)+nBl);
        topo_toplot=[];
        for nE=1:length(layout.label)-5
            temp1=table_alpha.Amp(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{nCond} & table_alpha.Elec==layout.label{nE});
            temp1b=table_alpha.SubID(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{nCond} & table_alpha.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1=grpstats(temp1,temp1b);
            topo_toplot(nE)=nanmean(mean1);
        end
        simpleTopoPlot_ft(topo_toplot', layout,'on',[],0,1);
        title(sprintf('B%g: %s',nBl,Conds{nCond}));
        colorbar;
%         caxis([.4 .8]);
        format_fig;
    end
end

for nBl=1:5
    subplot(3,5,10+nBl);
    mean1=[];
    mean2=[];
        for nE=1:length(layout.label)-5
            temp1=table_alpha.Amp(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{1} & table_alpha.Elec==layout.label{nE});
            temp1b=table_alpha.SubID(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{1} & table_alpha.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean1(:,nE)=grpstats(temp1,temp1b);
            
            temp1=table_alpha.Amp(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{2} & table_alpha.Elec==layout.label{nE});
            temp1b=table_alpha.SubID(table_alpha.BlockN==nBl & table_alpha.Cond==Conds{2} & table_alpha.Elec==layout.label{nE}); temp1b=removecats(temp1b);
            mean2(:,nE)=grpstats(temp1,temp1b);
        end
    [h,p,~,stats]=ttest2(mean1,...
        mean2);
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nBl,'D - E'));
    colorbar;
    caxis([-1 1]*2.5);
    format_fig;
end


%%
table_peaks=array2table(av_fooof_peaks,'VariableNames',{'SubID','BlockN','ElecN','CondF','Freq','Amp','BW'});
table_thetaalpha=array2table(av_fooof_maxalphatheta,'VariableNames',{'SubID','BlockN','ElecN','CondF','Freq','Amp','BW'});

table_thetaalpha.SubID=categorical(table_thetaalpha.SubID);
table_thetaalpha.Cond=table_thetaalpha.CondF;
table_thetaalpha.Cond=categorical(table_thetaalpha.Cond);
table_thetaalpha.Cond(table_thetaalpha.Cond=='1')='E';
table_thetaalpha.Cond(table_thetaalpha.Cond=='0')='D';
table_thetaalpha.Cond=removecats(table_thetaalpha.Cond);
table_thetaalpha.Elec=table_thetaalpha.ElecN;
table_thetaalpha.Elec=categorical(table_thetaalpha.Elec);
for nEl=1:length(layout.label)-5
    table_thetaalpha.Elec(table_thetaalpha.Elec==num2str(nEl))=layout.label{nEl};
end
table_thetaalpha.Elec=removecats(table_thetaalpha.Elec);

figure;
for nBl=1:5
    subplot(1,5,nBl)
    histogram(table_thetaalpha.Freq(table_thetaalpha.BlockN==nBl & table_thetaalpha.CondF==0),1:0.5:30,'Normalization','probability','FaceColor',Colors{1}(5,:))
    hold on
    histogram(table_thetaalpha.Freq(table_thetaalpha.BlockN==nBl & table_thetaalpha.CondF==1),1:0.5:30,'Normalization','probability','FaceColor',Colors{2}(5,:))
    format_fig;
    xlabel('Peak Frequency (Hz)')
    ylabel('Probability')
    legend({'R','B'})
    xlim([5 14])
end

mdlc_0=fitlme(table_thetaalpha,'Freq~1+(1|SubID)');
mdlc_1=fitlme(table_thetaalpha,'Freq~1+BlockN+(1|SubID)');
mdlc_2=fitlme(table_thetaalpha,'Freq~1+BlockN+Elec+(1|SubID)');
mdlc_3=fitlme(table_thetaalpha,'Freq~1+BlockN+Cond+(1|SubID)');
mdlc_4=fitlme(table_thetaalpha,'Freq~1+BlockN*Cond+(1|SubID)');

% figure; set(gcf,'Position',[62   198   548   787]);
% subplot(2,1,1); hold on;
% for nBl=1:5
%     for nCond=1:2
%         temp1=table_thetaalpha.Freq(table_thetaalpha.BlockN==nBl & table_thetaalpha.Cond==Conds{nCond});
%         temp1b=table_thetaalpha.SubID(table_thetaalpha.BlockN==nBl & table_thetaalpha.Cond==Conds{nCond}); temp1b=removecats(temp1b);
%         mean1=grpstats(temp1,temp1b);
%             simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
% %         line([1 1]*nBl+(2*nCond-3)*0.2,[-1 1].*sem(mean1)+mean(mean1),'Color',Colors{nCond}(nBl,:),'LineWidth',2);
% %         line([Pos-0.1*widthBar Pos+0.1*widthBar],[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError)
% %         line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError)
% %         
% %         scatter(nBl+(2*nCond-3)*0.2,mean(mean1),'MarkerEdgeColor',Colors{nCond}(nBl,:),'MarkerFaceColor',Colors{nCond}(nBl,:),'LineWidth',2);
%     end
% end
% title('Theta/Alpha')
% xlabel('Block')
% ylabel('Freq')
% format_fig;
% % ylim([0.7 1.1])
% xlim([0.5 5.5])
% 
% subplot(2,1,2);
% for nBl=1:5
%     for nCond=1:2
%         temp1=table_thetaalpha.Amp(table_thetaalpha.BlockN==nBl & table_thetaalpha.Cond==Conds{nCond});
%         temp1b=table_thetaalpha.SubID(table_thetaalpha.BlockN==nBl & table_thetaalpha.Cond==Conds{nCond}); temp1b=removecats(temp1b);
%         mean1=grpstats(temp1,temp1b);
% %         simpleBarPlot(nBl+(2*nCond-3)*0.2,mean1,Colors{nCond}(nBl,:),0.35,'k',[],3);
%             simpleDotPlot(nBl+(2*nCond-3)*0.1,mean1,144,Colors{nCond}(nBl,:),0.65,Colors{nCond}(nBl,:),[],3);
%     end
% end
% title('Theta/Alpha')
% format_fig;
% xlabel('Block')
% ylabel('Amplitude')
% % ylim([-0.5 -0])
% xlim([0.5 5.5])


% look at proportion, amplitude and freq of theta and alpha band peaks by
% electrodes and condition

% combining alpha and theta, look at the highest peak and what is it's
% frequency