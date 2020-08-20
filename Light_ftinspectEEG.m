%%
run localdef_ligthEEG.m
% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath((path_fieldtrip));
ft_defaults;

%% choose and load subject
List_Subj=dir([data_path filesep 'CIfIfre_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
load([data_path filesep ListNames{pick}])
oridata=data;

%% display data
cfg=[];
cfg.continuous='no';
cfg.allowoverlap='true';
cfg.viewmode='vertical';
cfg = ft_databrowser(cfg, oridata);

% %% reject trials
% cfg          = [];
% cfg.method   = 'summary';
% cfg.alim     = 5e-5;
% data        = ft_rejectvisual(cfg,oridata);

