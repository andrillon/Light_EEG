%%
run localdef_ligthEEG.m
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath((path_fieldtrip));
ft_defaults;
load('cain_elecloc_32ch_layout.mat');

%% choose and load subject
List_Subj=dir([data_path filesep 'Ifre_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
load([data_path filesep ListNames{pick}])
oridata=data;

%%
figure;
cfg = [];
cfg.component = 1:32;       % specify the component(s) that should be plotted
cfg.layout    = layout; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp)
set(gcf,'Position',[1           1        1871         984]);


%% reject trials
pickComponents=input('Select component you want to plot: ');

cfg          = [];
cfg.channel  = pickComponents; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = layout; % specify the layout file that should be used for plotting
cfg.allowoverlap='true';
ft_databrowser(cfg, comp)

