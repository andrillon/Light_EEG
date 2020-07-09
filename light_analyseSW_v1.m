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
List_Subj=dir([data_path filesep 'SW_fix37uV_CIfIfe_*.mat']);

%% Loop across participants to extract power
SW_properties=[];
SW_properties_all=[];
SW_density=[];
all_slow_Waves=[];
for nS=1:length(List_Subj)
    
    %%% load data
    File_Name = List_Subj(nS).name;
    %     if strcmp(File_Name,'TFe_ft_DLT018.mat')
    %         continue;
    %     end
    fprintf('... processing %s (%g/%g)',File_Name,nS,length(List_Subj))
    File_Path = List_Subj(nS).folder;
    load([data_path filesep File_Name]);
    slow_Waves(slow_Waves(:,5)>120000,:)=[];
    
    %%% extract info
    bound{1}=findstr(File_Name,'_');
    bound{2}=findstr(File_Name,'.');
    CodeSubj=File_Name(bound{1}(end)+1:bound{2}(1)-1);
    CondSubj(nS)=SubInfo.Condition(find(~cellfun(@isempty,regexpi(SubInfo.PT_Code,CodeSubj))));
    fprintf('... condition %s\n',CondSubj(nS))
    
    all_slow_Waves=[all_slow_Waves ; slow_Waves];
    for nBl=1:5
        temp_slow_Waves=slow_Waves(slow_Waves(:,2)==nBl,:);
        nout=histc(temp_slow_Waves(:,3),1:length(labels));
        if size(nout,1)>size(nout,2)
            nout=nout';
        end
        nout=nout/5; % SW density in minutes
        
        SW_density=[SW_density ; [nS nBl CondSubj(nS)=='E' nout]];
        for nE=1:32
            SW_properties=[SW_properties ; [nS nBl CondSubj(nS)=='E' nE sum(temp_slow_Waves(:,3)==nE) mean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,[4 9 11 12 13]),1)]];
        end
        SW_properties_all=[SW_properties_all ; [nS nBl CondSubj(nS)=='E' nE size(temp_slow_Waves,1) mean(temp_slow_Waves(:,[4 9 11 12 13]),1)]];
        
    end
    
end

%%
table_SW=array2table(SW_properties,'VariableNames',{'SubID','BlockN','Cond','Elec','NumW','P2P','NegP','PosP','negS','posS'});
table_SWall=array2table(SW_properties_all,'VariableNames',{'SubID','BlockN','Cond','Elec','NumW','P2P','NegP','PosP','negS','posS'});

table_SW.SubID=categorical(table_SW.SubID);
table_SW.Cond=categorical(table_SW.Cond);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.Cond(table_SW.Cond=='1')='E';
table_SW.Cond(table_SW.Cond=='0')='D';
table_SW.Cond=removecats(table_SW.Cond);

table_SWall.SubID=categorical(table_SWall.SubID);
table_SWall.Cond=categorical(table_SWall.Cond);
table_SWall.Elec=categorical(table_SWall.Elec);
table_SWall.Cond(table_SWall.Cond=='1')='E';
table_SWall.Cond(table_SWall.Cond=='0')='D';
table_SWall.Cond=removecats(table_SWall.Cond);

NumW_mdl0=fitlme(table_SW,'NumW~1+BlockN+(1|SubID)');
NumW_mdl1=fitlme(table_SW,'NumW~1+BlockN*Cond+(1|SubID)');

P2P_mdl0=fitlme(table_SW,'P2P~1+BlockN+(1|SubID)');
P2P_mdl1=fitlme(table_SW,'P2P~1+BlockN*Cond+(1|SubID)');

% NegP_mdl0=fitlme(table_SW,'NegP~1+BlockN+(1|SubID)');
% NegP_mdl1=fitlme(table_SW,'NegP~1+BlockN*Cond+(1|SubID)');
%
% PosP_mdl0=fitlme(table_SW,'PosP~1+BlockN+(1|SubID)');
% PosP_mdl1=fitlme(table_SW,'PosP~1+BlockN*Cond+(1|SubID)');
%
% negS_mdl0=fitlme(table_SW,'negS~1+BlockN+(1|SubID)');
% negS_mdl1=fitlme(table_SW,'negS~1+BlockN*Cond+(1|SubID)');
%
% posS_mdl0=fitlme(table_SW,'posS~1+BlockN+(1|SubID)');
% posS_mdl1=fitlme(table_SW,'posS~1+BlockN*Cond+(1|SubID)');
%%
figure;
Cond={'D','E'};
data=[];
for nBl = 1:5
    for nC = 1:2
        data{nBl, nC} = (table_SW.NumW(table_SW.Elec == '13' & table_SW.Cond == Cond(nC) & table_SW.BlockN == (nBl)));
        %          data{nBl, nC} = (table_SWall.NumW(table_SWall.Cond == Cond(nC) & table_SWall.BlockN == (nBl)));
        meandata(nBl, nC) =nanmean(data{nBl, nC}./data{1, nC});
        semdata(nBl, nC) =sem(data{nBl, nC}./data{1, nC});
        
    end
end
hold on;
hb(1)=errorbar((1:5)-0.05,meandata(:,1),semdata(:,1),'Color',Colors{2}(end,:),'LineWidth',3);
hb(2)=errorbar((1:5)+0.05,meandata(:,2),semdata(:,2),'Color',Colors{1}(end,:),'LineWidth',3);
scatter((1:5)-0.05,meandata(:,1),'Marker','o','SizeData',144,'MarkerEdgeColor',Colors{2}(end,:),'MarkerFaceColor',Colors{2}(end,:),'MarkerFaceAlpha',0.7);
scatter((1:5)+0.05,meandata(:,2),'Marker','o','SizeData',144,'MarkerEdgeColor',Colors{1}(end,:),'MarkerFaceColor',Colors{1}(end,:),'MarkerFaceAlpha',0.7);

legend(hb,{'D','E'})
format_fig;
xlabel('Block')
ylabel('Num SW (norm)')
% % make figure
% figure;
% cl=[Colors{2}(end,:) ; Colors{1}(end,:)];
% h   = rm_raincloud(data, cl);
% % set(gca, 'YLim', [-0.3 1.6]);
%
% data2=data;
% for nBl = 2:5
%     for nC = 1:2
%         data2{nBl, nC} = data2{nBl, nC}./data2{1, nC}*100;
%     end
% end
%
% % make figure
% figure;
% cl=[Colors{2}(end,:) ; Colors{1}(end,:)];
% h   = rm_raincloud(data2(2:5,:), cl);
% % set(gca, 'YLim', [-0.3 1.6]);

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


%%
%%
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
Conds={'D','E'};
Names={'A','D','E'};
figure; set(gcf,'Position',[64          33        1097         450]);
for np=1:3
    subplot(1,3,np);
    toplot=[];
    for nE=1:32
        if np==1
            toplot(nE)=squeeze(mean(table_SW.NumW(table_SW.BlockN>1 & table_SW.Elec==num2str(nE))));
        elseif np==2
            toplot(nE)=squeeze(mean(table_SW.NumW(table_SW.BlockN>1 & table_SW.Elec==num2str(nE) & table_SW.Cond=='D')));
        elseif np==3
            toplot(nE)=squeeze(mean(table_SW.NumW(table_SW.BlockN>1 & table_SW.Elec==num2str(nE) & table_SW.Cond=='E')));
        end
    end
    simpleTopoPlot_ft(toplot', layout,'on',[],0,1);
    title(Names{np});
    colormap(cmap);
    colorbar;
    caxis([0 1]*30);
    format_fig;
end

%%
figure; set(gcf,'Position',[64          33        1097         450]);
for np=1:3
    subplot(1,3,np);
    toplot=[];
    for nE=1:32
        if np==1
            toplot(nE)=squeeze(nanmean(table_SW.P2P(table_SW.BlockN>1 & table_SW.Elec==num2str(nE))));
        elseif np==2
            toplot(nE)=squeeze(nanmean(table_SW.P2P(table_SW.BlockN>1 & table_SW.Elec==num2str(nE) & table_SW.Cond=='D')));
        elseif np==3
            toplot(nE)=squeeze(nanmean(table_SW.P2P(table_SW.BlockN>1 & table_SW.Elec==num2str(nE) & table_SW.Cond=='E')));
        end
    end
    simpleTopoPlot_ft(toplot', layout,'on',[],0,1);
    title(Names{np});
    colormap(cmap);
    colorbar;
%     caxis([350 750]);
    format_fig;
end