%% init
clear all
close all

% data_path='/Users/tand0009/Data/Cain_Light/';
data_path='/Volumes/tLab_BackUp1/Monash/Cain_Light/';
sw_path='/Users/tand0009/Data/Cain_Light/SWdetection';

path_eeglab='/Users/tand0009/Work/local/eeglab14_1_2b/';
path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_export='/Users/tand0009/Work/local/export_fig/';

addpath(genpath(path_LSCPtools));
addpath(genpath(path_export));
addpath(path_fieldtrip);
ft_defaults;

addpath(genpath(path_eeglab))

%%
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep '**/*.eeg']);
nc=0;

%%
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
% fixThr=75;
fixThr=[];
art_ampl=150;
max_posampl=75;
max_Freq=7;

%%
elecs=readlocs('cain_elecloc_32ch.locs');
myLabels={elecs.labels};
for nS=1:length(List_Subj)
    
    %%% load waves
    File_Name = List_Subj(nS).name;
    if strcmp(File_Name,'DLT001.eeg') || strcmp(File_Name,'DLT004.eeg') || strcmp(File_Name,'DLT016.eeg') || strcmp(File_Name,'DLT024.eeg') || strcmp(File_Name,'DLT026.eeg') || exist([sw_path filesep File_Name(1:end-4) '_SW_all.mat'])==0 || nS==7
        continue;
    end
    File_Path = List_Subj(nS).folder;
    nc=nc+1;
    fprintf('... processing %s (%g/%g)\n',File_Name,nS,length(List_Subj))
    load([sw_path filesep File_Name(1:end-4) '_SW_all'],'all_Waves')
    hdr     = ft_read_header([File_Path filesep File_Name]);
    
    if SubInfo.Condition(match_str(SubInfo.PT_Code,File_Name(1:end-4)))=='D'
        condition(nc) = 0;
    elseif SubInfo.Condition(match_str(SubInfo.PT_Code,File_Name(1:end-4)))=='E'
        condition(nc) = 1;
    end
    
    %%% clean waves
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,AmpCriterionIdx)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,AmpCriterionIdx)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl,:)=[];
    
    %%% set threshold on first block (baseline)
    EEG_channels=myLabels;
    EEG_channels(match_str(EEG_channels,{'x_dir','y_dir','z_dir','ECG1','ECG2','EOG R','EOG L'}))=[];
    pick_EEGchannels=find(ismember(hdr.label,EEG_channels));
    
    
    %%% collect number of waves per block
    for nE=1:length(myLabels)
        if ismember(myLabels{nE},EEG_channels) && sum(ismember(hdr.label(all_Waves(:,3)),myLabels{nE}))~=0
            %             thr_Wave(nc,nE)=prctile(all_Waves(all_Waves(:,2)==1 & ismember(all_Waves(:,3),pick_EEGchannels),AmpCriterionIdx),prticle_Thr);
            thr_Wave(nc,nE)=prctile(all_Waves(all_Waves(:,2)==1 & all_Waves(:,3)==nE,AmpCriterionIdx),prticle_Thr);
            slow_Waves=all_Waves(all_Waves(:,AmpCriterionIdx)>thr_Wave(nc,nE),:);
            for nb=1:5
                dens_Waves(nc,nE,nb)=sum(slow_Waves(:,2)==nb & ismember(hdr.label(slow_Waves(:,3)),myLabels{nE}));
            end
        else
            thr_Wave(nc,nE)=nan;
            for nb=1:5
                dens_Waves(nc,nE,nb)=nan;
            end
        end
    end
end


%% figures
dens_Waves(:,sum(isnan(squeeze(nanmean(dens_Waves,3))))>0,:)=nan;

dens_Waves2=dens_Waves; %./repmat(dens_Waves(:,:,1),[1 1 size(dens_Waves,3)])*100;
figure; hp=[];
format_fig;
this_ch=match_str(hdr.label,'Fz');
hp(1)=errorbar(1:5,squeeze(mean(dens_Waves2(condition==0,this_ch,:),1)),squeeze(sem(dens_Waves2(condition==0,this_ch,:),1)),'LineWidth',2,'Color','r');
hold on;
scatter(1:5,squeeze(mean(dens_Waves2(condition==0,this_ch,:),1)),'LineWidth',2,'SizeData',72,'MarkerEdgeColor','r','MarkerFaceColor','w');


hp(2)=errorbar(1:5,squeeze(mean(dens_Waves2(condition==1,this_ch,:),1)),squeeze(sem(dens_Waves2(condition==1,this_ch,:),1)),'LineWidth',2,'Color','b');
hold on;
scatter(1:5,squeeze(mean(dens_Waves2(condition==1,this_ch,:),1)),'LineWidth',2,'SizeData',72,'MarkerEdgeColor','b','MarkerFaceColor','w');

set(gca,'XTick',1:5,'XTickLabel',{'bsl','b1','b2','b3','b4'})
format_fig;
ylabel('# of Waves on Fz')
legend(hp,{'D','E'})
xlim([0.5 5.5])
%%
fprintf('%2.0f/%2.0f\n',0,35)
clear topo_RManova topo_anova
for nE=1:35
    if ismember(hdr.label{nE},EEG_channels) && sum(isnan(squeeze(dens_Waves(:,nE,2))))~=length(squeeze(dens_Waves(:,nE,2)))
        fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nE,35)
        t = table(condition',squeeze(dens_Waves(:,nE,2)),squeeze(dens_Waves(:,nE,3)),squeeze(dens_Waves(:,nE,4)),squeeze(dens_Waves(:,nE,5)),...
            'VariableNames',{'Condition','meas1','meas2','meas3','meas4'});
        Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
        rm = fitrm(t,'meas1-meas4~Condition','WithinDesign',Meas);
        ranovatbl = ranova(rm);
        anovatbl = anova(rm);
        topo_RManova(:,nE)=[ranovatbl.F];
        topo_anova(:,nE)=[anovatbl.F];
        
        topo_RManova_pV(:,nE)=[ranovatbl.pValue];
        topo_anova_pV(:,nE)=[anovatbl.pValue];
    else
        topo_RManova(:,nE)=nan;
        topo_anova(:,nE)=nan;
        topo_RManova_pV(:,nE)=nan;
        topo_anova_pV(:,nE)=nan;
    end
end

%% Figure%%

%%
maxA=25/4;
figure; set(gcf,'Position',[14   141   264   657])
subplot(3,1,1); format_fig;
topo=squeeze(nanmean(nanmean(dens_Waves(condition==0,:,:),3),1))/4;
topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
% title('D')
colormap('parula')
colorbar; caxis([0 maxA])

subplot(3,1,2); format_fig;
topo=squeeze(nanmean(nanmean(dens_Waves(condition==1,:,:),3),1))/4;
topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
% title('E')
colormap('parula')
colorbar; caxis([0 maxA])

subplot(3,1,3); format_fig;
topo=squeeze(topo_anova(2,:));
topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
% title({'F-Value - Light Condition'})
colorbar; caxis([-1 1]*2.1)
colormap('parula')

export_fig(['/Users/tand0009/Work/Documents/Grants/ARC/2020/Light_DP2020/' filesep 'Fig2.fig'])
export_fig(['/Users/tand0009/Work/Documents/Grants/ARC/2020/Light_DP2020/' filesep 'Fig2.eps'],'-r 300')