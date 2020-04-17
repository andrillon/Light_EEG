function trl = lightEEG_trialfun(cfg);

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

trl = [];

% clean events
findsep=findstr(cfg.dataset,'/');
findext=findstr(cfg.dataset,'.');
File_Name=cfg.dataset(findsep(end)+1:findext-1);

if strcmpi(File_Name,'DLT001')
    event(4).value='Baseline';
    event(5).value='FG1'; event(5).type='Comment'; event(5).sample=583841;
    event(6).value='CT1'; event(6).sample=793500;
    event(7).value='FG2'; event(7).sample=883500;
    event(8).value='CT2'; event(8).sample=1033500;
    event(9).value='FG3'; event(9).sample=1123500;
    event(10).value='CT3'; event(10).sample=1273500;
    event(11).value='FG4'; event(11).sample=1363500;
    event(12).value='CT4'; event(12).sample=1513500;
elseif strcmpi(File_Name,'DLT004')
    event([5 6 9])=[];
elseif strcmpi(File_Name,'DLT016')
    event(9).value='CT3';
    event(10).value='FG4';
elseif strcmpi(File_Name,'DLT018')
    event(2).value='Baseline';
    event(3).value='FG1';
    event(6).value='CT2';
elseif strcmpi(File_Name,'DLT024')
    event(3).value='FG1'; event(2).value='Baseline';
elseif strcmpi(File_Name,'DLT026')
    event(1)=[];
    event(find(~(cellfun(@isempty,regexp({event.value},'^LostSamples')))))=[];
elseif strcmpi(File_Name,'DLT038')
    event(3)=[];
end
for ne=1:length(event)
    if isempty(event(ne).value)
        event(ne).value='null';
    end
end

for i=1:length(event)
% if strcmp(event(i).type, cfg.trialdef.eventtype)
  % it is a trigger, see whether it has the right value
  if ismember(event(i).value, cfg.trialdef.eventvalue)
    % add this to the trl definition
    begsample     = event(i).sample - cfg.trialdef.prestim*hdr.Fs;
    endsample     = event(i).sample + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;   
    trl(end+1, :) = [round([begsample endsample offset])];
  end
% end
end
