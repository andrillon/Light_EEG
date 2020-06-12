%% initiliase - clear all variables and scripts
clear all
close all

%% set up path for project using local file
run localdef_ligthEEG.m

addpath(genpath(path_LSCPtools)); % Thomas' general toolkit
addpath(path_fieldtrip); % Filedtrip toolbox (EEG)
ft_defaults; % Set up fieldtrip toolbox
addpath(genpath(fooof_path))
addpath(genpath(path_raincloud))

%% List files and retrieve layout
load('light_subinfo.mat');
load('cain_elecloc_32ch_layout.mat');
List_Subj=dir([data_path filesep 'SW_CIfIfe_*.mat']);

%% Loop across participants to extract power
SW_properties=[];
SW_density=[];
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
    
    
    for nBl=1:5
        temp_slow_Waves=slow_Waves(slow_Waves(:,2)==nBl,:);
        nout=histc(temp_slow_Waves(:,3),1:length(labels));
        if size(nout,1)>size(nout,2)
            nout=nout';
        end
        nout=nout/5; % SW density in minutes
        
        SW_density=[SW_density ; [nS nBl CondSubj(nS)=='E' nout]];
    end
    SW_properties=[SW_properties ; [slow_Waves repmat(CondSubj(nS)=='E',size(slow_Waves,1),1)]];
    
end

%%
figure;

data=[];
for nBl = 1:5
    for nC = 1:2
        data{nBl, nC} = mean(SW_density(SW_density(:, 3) == nC-1 & SW_density(:, 2) == nBl,4:end),2);
    end
end

% make figure
figure;
cl=[Colors{2}(end,:) ; Colors{1}(end,:)];
h   = rm_raincloud(data, cl);
% set(gca, 'YLim', [-0.3 1.6]);

data2=data;
for nBl = 2:5
    for nC = 1:2
        data2{nBl, nC} = data2{nBl, nC}./data2{1, nC}*100;
    end
end

% make figure
figure;
cl=[Colors{2}(end,:) ; Colors{1}(end,:)];
h   = rm_raincloud(data2(2:5,:), cl);
% set(gca, 'YLim', [-0.3 1.6]);

%%
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)

figure; set(gcf,'Position',[64          33        1097         952]);
Conds={'D','E'};
for nC=1:2
    for nBl=1:5
        subplot(2,5,5*(nC-1)+nBl);
        Dens_AVG=squeeze(mean(SW_density(SW_density(:, 3) == nC-1 & SW_density(:, 2) == nBl,4:end),1));
        
        simpleTopoPlot_ft(Dens_AVG', layout,'on',[],0,1);
        title(sprintf('%s - %g',Conds{nC},nBl));
        colormap(cmap);
        colorbar;
        caxis([0 1]*6);
        format_fig;
    end
end

%%
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)

figure; set(gcf,'Position',[64          33        1097         450]);
Conds={'D','E'};
for nBl=1:5
    subplot(1,5,nBl);
    Dens_AVG=squeeze(mean(SW_density(SW_density(:, 3) == 1 & SW_density(:, 2) == nBl,4:end),1))-...
        squeeze(mean(SW_density(SW_density(:, 3) == 0 & SW_density(:, 2) == nBl,4:end),1));
    
    simpleTopoPlot_ft(Dens_AVG', layout,'on',[],0,1);
    title(sprintf('E-D - %g',nBl));
    colormap(cmap);
    colorbar;
    caxis([-1 1]*3);
    format_fig;
end
