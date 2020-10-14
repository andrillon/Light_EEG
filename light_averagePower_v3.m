%% initiliase - clear all variables and scripts
clear all
% close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep 'laplTFCIfre_ft_*.mat']);

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
%     if length(TFRhann.label)~=32
%         warning(sprintf('not enough channels (n=%g)',length(TFRhann.label)))
%         continue;
%     end
%     if nS==1
%         av_logPower=nan([length(List_Subj) size(TFRhann.powspctrm,1) size(TFRhann.powspctrm,2) size(TFRhann.powspctrm,3)]);
%     end
    if size(TFRhann.powspctrm,1)==4
        continue;
    end
    %     av_logPower(nS,:,:,:)=(squeeze(mean(log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,TFRhann.freq>16,:),3),[1 1 size(TFRhann.powspctrm,3) 1])),4)));
    av_logPower(nS,:,:,:)=(squeeze(10*log10(mean(TFRhann.powspctrm(:,:,:,TFRhann.time<-30),4))));
    %     av_logPower_time(nS,:,:,:,:)=(squeeze(log(TFRhann.powspctrm(:,:,:,TFRhann.time<-30))));
    %     TFRhann.powspctrm_norm=10*log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(1,:,:,:),4),[size(TFRhann.powspctrm,1) 1 1 size(TFRhann.powspctrm,4)]));
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
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
%     ylim([-20 5])
end


%%
theseChLabels={'Fz','Cz','Oz'};
figure; set(gcf,'Position',[ 1           7        1482         798]);
for nBl=1:5
    for nChan=1:length(theseChLabels)
        subplot(3,5,5*(nChan-1)+nBl);
        thisChLabel=theseChLabels{nChan};
        for nCond=1:2
            hold on;
            hp=[];
            temp_pow=squeeze(av_logPower(CondSubj==Conds{nCond},nBl,match_str(TFRhann.label,thisChLabel),:));
            if nCond==1 && nBl==1
                %                 simpleTplot(freqs,(temp_pow),0,Colors{nCond}(nBl,:),0,'--',0.5,1,0,0,2);
                plot(freqs,mean(temp_pow),'Color',Colors{nCond}(nBl,:),'LineStyle','--','LineWidth',3)
            else
                %                 simpleTplot(freqs,(temp_pow),0,Colors{nCond}(nBl,:),0,'-',0.5,1,0,0,2);
                plot(freqs,mean(temp_pow),'Color',Colors{nCond}(nBl,:),'LineStyle','-','LineWidth',3)
            end
            format_fig;
            if nChan==length(theseChLabels)
                xlabel('Freq (Hz)');
            end
            if nBl==1
                ylabel(sprintf('log(Power) %s',theseChLabels{nChan}));
            end
            %             %     legend(hp,{'B0','B1','B2','B3','B4'});
            if nChan==1
                title(sprintf('Block %g',nBl));
            end
%             ylim([-12 0])
%             xlim([5 15])
        end
    end
end

%%
FOI=[20 30]; % Freq Band of Interest
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nB=1:5
        subplot(3,5,5*(nCond-1)+nB);
        Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{nCond},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
        
        simpleTopoPlot_ft(Pow_AVG, layout,'on',[],0,1);
        title(sprintf('B%g: %s',nB,Conds{nCond}));
        colorbar;
%         caxis([-20 -10]);
        %         caxis([-1 1]*max(abs(Pow_AVG)));
        format_fig;
    end
end

for nB=1:5
    subplot(3,5,10+nB);
    Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1))-...
        squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
    [h,p,~,stats]=ttest2(squeeze(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2)),...
        squeeze(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2)));
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nB,'D - E'));
    colorbar;
    %     caxis([-1 1]*.5);
        caxis([-1 1]*3);
    format_fig;
end

%%
FOI=[7 9]; % Freq Band of Interest
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nB=1:5
        subplot(3,5,5*(nCond-1)+nB);
        Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{nCond},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
        
        simpleTopoPlot_ft(Pow_AVG, layout,'on',[],0,1);
        title(sprintf('B%g: %s',nB,Conds{nCond}));
        colorbar;
        caxis([-10 -3]);
        %         caxis([-1 1]*max(abs(Pow_AVG)));
        format_fig;
    end
end

for nB=1:5
    subplot(3,5,10+nB);
    Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1))-...
        squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2),1));
    [h,p,~,stats]=ttest2(squeeze(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2)),...
        squeeze(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>FOI(1) & freqs<FOI(2)),4),2)));
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nB,'D - E'));
    colorbar;
    %     caxis([-1 1]*.5);
        caxis([-1 1]*3);
    format_fig;
end

%%
FOI=[7 10]; % Freq Band of Interest
figure; set(gcf,'Position',[64         311        1488         674]);
tight_subplot(3,5,[.01 .03],[.1 .01],[.01 .01])
for nCond=1:2
    for nB=1:5
        subplot(3,5,5*(nCond-1)+nB);
        Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{nCond},nB,:,freqs>=FOI(1) & freqs<FOI(2)),4),2),1));
        
        simpleTopoPlot_ft(Pow_AVG, layout,'on',[],0,1);
        title(sprintf('B%g: %s',nB,Conds{nCond}));
        colorbar;
                        caxis([-10 -3]);
        %         caxis([-1 1]*max(abs(Pow_AVG)));
        format_fig;
    end
end

for nB=1:5
    subplot(3,5,10+nB);
    Pow_AVG=squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>=FOI(1) & freqs<FOI(2)),4),2),1))-...
        squeeze(mean(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>=FOI(1) & freqs<FOI(2)),4),2),1));
    [h,p,~,stats]=ttest2(squeeze(mean(mean(av_logPower(CondSubj==Conds{1},nB,:,freqs>=FOI(1) & freqs<FOI(2)),4),2)),...
        squeeze(mean(mean(av_logPower(CondSubj==Conds{2},nB,:,freqs>=FOI(1) & freqs<FOI(2)),4),2)));
    simpleTopoPlot_ft(stats.tstat', layout,'on',[],0,1);
    title(sprintf('B%g: %s',nB,'D - E'));
    colorbar;
    %     caxis([-1 1]*.5);
        caxis([-1 1]*2);
    format_fig;
end
