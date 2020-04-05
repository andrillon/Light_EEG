%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox

%% Parameters Power Extraction
w_window  = 4; % in seconds
w_overlap = (w_window/2);

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep '**/*.eeg']);

%% Loop across participants to extract power
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
    
    if strcmp(File_Name,'DLT001.eeg')
        continue;
    elseif strcmp(File_Name,'DLT004.eeg')
        events([5 6 9])=[];
    elseif strcmp(File_Name,'DLT016.eeg')
        events(9).value='CT3';
        events(10).value='FG4';
    elseif strcmp(File_Name,'DLT024.eeg')
        events(3).value='FG1'; events(2)=[];
    elseif strcmp(File_Name,'DLT026.eeg')
        events(1)=[];
        events(find(~(cellfun(@isempty,regexp({events.value},'^LostSamples')))))=[];
    end
    
   
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
    % re-reference to av. mastoids
    %     data=data-repmat(mean(data([10 21],:),1),[size(data,1), 1]);
    % re-reference to the average
    data=data-repmat(mean(data(:,:),1),[size(data,1), 1]);
    
   
    % loop across electrodes
    EEG_channels=hdr.label;
    EEG_channels(match_str(EEG_channels,{'x_dir','y_dir','z_dir','ECG1','ECG2','EOG R','EOG L'}))=[];
    pick_channels=hdr.label; %find(ismember(hdr.label,EEG_channels));
    for nelec=1:length(pick_channels)
        if ismember(pick_channels{nelec},EEG_channels)
            fprintf('... ... power extraction elec %s (%g/%g)\n',hdr.label{nelec},nelec,length(pick_channels))
            % bandpass between 0.1 and 30Hz
            data(nelec,:)=bandpass(data(nelec,:), hdr.Fs, 0.1, 30, 4);
            
            % baseline power
            this_data=squeeze(data(nelec,baseline_end+((-4*60*hdr.Fs+1):0))); % select 4 minutes before the end of the FG period
            this_data=this_data-mean(this_data);
            [pow,faxis] = pwelch(this_data,w_window*hdr.Fs,round(w_overlap*hdr.Fs),[],hdr.Fs,'power');
            pow_baseline(nc,nelec,:)=pow';
            
            %     pow_fix=[];
            for nb=1:4
                this_data=squeeze(data(nelec,fix_end(nb)+((-4*60*hdr.Fs+1):0))); % select 4 minutes before the end of the FG period
                this_data=this_data-mean(this_data);
                [pow_fix(nc,nelec,nb,:),faxis] = pwelch(this_data,w_window*hdr.Fs,round(w_overlap*hdr.Fs),[],hdr.Fs,'power');
            end
        else
            pow_baseline(nc,nelec,:)=nan;
            for nb=1:4
                pow_fix(nc,nelec,nb,:)=nan;
            end
        end
    end
    
     % plot data
    figure; imagesc(data);
    
end

%%
pow_baseline(15,:)=nan;
for nb=1:4
    relpow_fix(:,:,nb,:)=squeeze(pow_fix(:,:,nb,:))./pow_baseline*100;
end

%% Figure: Time on Task Effect
my_channels={'Fz','Cz','Pz','Oz'};
my_colors=[1 0 0; 0.8 0.2 0; 0.1 0.1 0.4; 0.1 0.4 0.6 ];
figure;
for nch=1:length(my_channels)
    ch_idx=find(ismember(pick_channels,my_channels{nch}));
    
    subplot(2,2,nch); format_fig;
    hold on
%     plot(faxis,squeeze(nanmean(pow_baseline(:,ch_idx,:),1)),'LineWidth',2);
    for nb=1:4
        plot(faxis,squeeze(nanmean(pow_fix(:,ch_idx,nb,:),1))','LineWidth',2,'Color',my_colors(nb,:));
    end
    xlim([3 20])
    format_fig;
    if nch==1
        legend({'FG1','FG2','FG3','FG4'})
    end
    % legend({'FG1','FG2','FG3','FG4'})
    xlabel('Freq (Hz)')
    ylabel('Power')
    title(my_channels{nch})
    
end

%%
theta_band=[7.5 10];
alpha_band=[10 12];
addpath(genpath(path_eeglab))
for nplot=2
    figure;
    if nplot==1
        subplot(1,2,1);
        topo=squeeze(nanmedian(nanmean(pow_baseline(:,:,faxis>theta_band(1) & faxis<theta_band(2)),3),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title('Theta')
         caxis([0 7])
       
        subplot(1,2,2);
        topo=squeeze(nanmedian(nanmean(pow_baseline(:,:,faxis>alpha_band(1) & faxis<alpha_band(2)),3),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title('Alpha')
         caxis([0 7])
   else
        subplot(2,1,1);
        topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(:,:,:,faxis>theta_band(1) & faxis<theta_band(2)),3),4),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title('Theta'); format_fig
        caxis([-2 2])
       
        subplot(2,1,2);
        topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(:,:,:,faxis>alpha_band(1) & faxis<alpha_band(2)),3),4),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title('Alpha'); format_fig
        caxis([-2 2])
    end
end
% rmpath((path_eeglab))

%%
my_channels={'Fz','Cz','Pz','Oz'};
figure;
for nch=1:length(my_channels)
    subplot(1,2,1); hold on;
    tempT=squeeze(nanmean(pow_fix(:,match_str(pick_channels,my_channels{nch}),:,faxis>theta_band(1) & faxis<theta_band(2)),4));
    errorbar(1:4,mean(tempT),sem(tempT),'Color',my_colors(nch,:),'LineWidth',2);
    scatter(1:4,mean(tempT),'MarkerFaceColor','w','MarkerEdgeColor',my_colors(nch,:),'SizeData',72);
    format_fig;
    xlim([0.5 4.5])
    set(gca,'XTick',1:4,'XTickLabel',{'FG1','FG2','FG3','FG4'})
    ylabel('Power'); title('Theta')
    
    subplot(1,2,2); hold on;
    tempT=squeeze(nanmean(pow_fix(:,match_str(pick_channels,my_channels{nch}),:,faxis>alpha_band(1) & faxis<alpha_band(2)),4));
    hb(nch)=errorbar(1:4,mean(tempT),sem(tempT),'Color',my_colors(nch,:),'LineWidth',2);
    scatter(1:4,mean(tempT),'MarkerFaceColor','w','MarkerEdgeColor',my_colors(nch,:),'SizeData',72);
    format_fig;
    xlim([0.5 4.5])
    set(gca,'XTick',1:4,'XTickLabel',{'FG1','FG2','FG3','FG4'})
    ylabel('Power'); title('Alpha')
end
  subplot(1,2,2);   legend(hb,my_channels);


%% Figure: Effect of light
my_channels={'Fz','Cz','Pz','Oz'};
my_colors=[1 0 0; 0.8 0 0; 0.6 0 0.1; 0.4 0 0.3];
figure;
condition_names={'D','E','E-D'};
for nch=1:length(my_channels)
    ch_idx=find(ismember(pick_channels,my_channels{nch}));
    
    for ncond=1:3
        subplot(4,3,ncond+3*(nch-1)); format_fig;
        hold on
        %         plot(faxis,squeeze(nanmean(pow_baseline(condition==ncond-1,ch_idx,:),1)),'LineWidth',2);
        for nb=1:4
            if ncond==3
                tempplot=squeeze(nanmean(pow_fix(condition==1,ch_idx,nb,:),1))'-...
                    squeeze(nanmean(pow_fix(condition==0,ch_idx,nb,:),1))';
            else
                tempplot=squeeze(nanmean(pow_fix(condition==ncond-1,ch_idx,nb,:),1))';
            end
            plot(faxis,tempplot,'LineWidth',2,'Color',my_colors(nb,:));
        end
        xlim([3 20])
        format_fig;
        if nch==1
%             legend({'FG1','FG2','FG3','FG4'})
        end
        % legend({'FG1','FG2','FG3','FG4'})
        if nch==4
        xlabel('Freq (Hz)')
        end
        if ncond==1
        ylabel('Power')
        end
        title(sprintf('%s - %s',condition_names{ncond},my_channels{nch}))
        if ncond==3
            ylim([-3 2])
        else
            ylim([0 5])
        end
    end
end

%%
my_channels={'Fz','Cz','Pz','Oz'};
my_colors=[1 0 0; 0.8 0 0; 0.6 0 0.1; 0.4 0 0.3];
figure;
condition_names={'D','E','E-D'};
for nch=1:length(my_channels)
    ch_idx=find(ismember(pick_channels,my_channels{nch}));
    
    for ncond=3
        subplot(2,2,nch); format_fig;
        hold on
        for nb=1:4
            if ncond==3
                tempplot=squeeze(nanmean(nanmean(pow_fix(condition==1,ch_idx,:,:),1),3))'-...
                    squeeze(nanmean(nanmean(pow_fix(condition==0,ch_idx,:,:),1),3))';
            else
                tempplot=squeeze(nanmean(pow_fix(condition==ncond-1,ch_idx,nb,:),1))';
            end
            plot(faxis,tempplot,'LineWidth',2,'Color',my_colors(nb,:));
        end
        xlim([3 20])
        format_fig;

        if nch>2
        xlabel('Freq (Hz)')
        end
        if ncond==1 || ncond==3
        ylabel('Power')
        end
        title(sprintf('%s - %s',condition_names{ncond},my_channels{nch}))
        if ncond==3
            ylim([-2 1])
        else
            ylim([0 5])
        end
        line(xlim,[0 0],'Color','k','LineStyle','--')
    end
end

% %%
% theta_band=[7.5 11];
% sigma_band=[11 15];
% for nplot=1:4
%     figure;
%     for ncond=1:2
%         subplot(2,3,3-(ncond));
%         topo=squeeze(nanmedian(nanmean(pow_fix(condition==ncond-1,:,nplot,faxis>theta_band(1) & faxis<theta_band(2)),4),1));
%         topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
%         title(sprintf('%s - %s',condition_names{ncond},'Theta'))
%         
%         subplot(2,3,6-(ncond));
%         topo=squeeze(nanmedian(nanmean(pow_fix(condition==ncond-1,:,nplot,faxis>sigma_band(1) & faxis<sigma_band(2)),4),1));
%         topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
%         title(sprintf('%s - %s',condition_names{ncond},'Alpha'))
%     end
%     subplot(2,3,3);
%     topo=squeeze(nanmedian(nanmean(pow_fix(condition==1,:,nplot,faxis>theta_band(1) & faxis<theta_band(2)),4),1))-...
%         squeeze(nanmedian(nanmean(pow_fix(condition==0,:,nplot,faxis>theta_band(1) & faxis<theta_band(2)),4),1));
%     topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
%     title(sprintf('%s - %s',condition_names{3},'Theta'))
%     
%     subplot(2,3,6);
%     topo=squeeze(nanmedian(nanmean(pow_fix(condition==1,:,nplot,faxis>sigma_band(1) & faxis<sigma_band(2)),4),1))-...
%         squeeze(nanmedian(nanmean(pow_fix(condition==0,:,nplot,faxis>sigma_band(1) & faxis<sigma_band(2)),4),1));
%     topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
%     title(sprintf('%s - %s',condition_names{3},'Alpha'))
% end
%%
theta_band=[7.5 10];
sigma_band=[11 15];
 figure;
  set(gcf,'Position',[ 318    49   922   756])
  for ncond=1:2
        subplot(2,3,3-(ncond));
        topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==ncond-1,:,:,faxis>theta_band(1) & faxis<theta_band(2)),3),4),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title(sprintf('%s - %s',condition_names{ncond},'Theta'))
          caxis([-2 2]); format_fig; colorbar;

        subplot(2,3,6-(ncond));
        topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==ncond-1,:,:,faxis>sigma_band(1) & faxis<sigma_band(2)),3),4),1));
        topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
        title(sprintf('%s - %s',condition_names{ncond},'Sigma'))
          caxis([-2 2]); format_fig; colorbar;
  end
    subplot(2,3,3);
    topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==1,:,:,faxis>theta_band(1) & faxis<theta_band(2)),3),4),1))-...
        squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==0,:,:,faxis>theta_band(1) & faxis<theta_band(2)),3),4),1));
    topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
    title(sprintf('%s - %s',condition_names{3},'Theta'))
    caxis([-1 1]*0.5); format_fig; colorbar;

    subplot(2,3,6);
    topo=squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==1,:,:,faxis>sigma_band(1) & faxis<sigma_band(2)),3),4),1))-...
        squeeze(nanmedian(nanmean(nanmean(pow_fix(condition==0,:,:,faxis>sigma_band(1) & faxis<sigma_band(2)),3),4),1));
    topoplot(topo', 'cain_elecloc_32ch.locs','style','both','whitebk','on','electrodes','off');
    title(sprintf('%s - %s',condition_names{3},'Sigma'))
caxis([-1 1]); format_fig; colorbar;

%%
% figure;
% totwin=2:1:20;
% for nwin=totwin
% [pow_baseline,faxis] = pwelch(this_data,nwin*hdr.Fs,[],[],hdr.Fs,'power');
%     plot(faxis,pow_baseline); xlim([3 20])
%     hold on;
%     pause;
% end

%%
% elec = ft_read_sens('cain_elecloc_32ch.sfp');
% cfg=[];
% cfg.elec=elec;
% [layout] = ft_prepare_layout(cfg);
% save(['cain_elecloc_32ch_layout'],'layout');