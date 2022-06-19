function [recording, filename] = getrawandwfs(unit, spknum, map)

% first argument is which unit you're looking for 
% ('[unit]_wf' almost always has more info, so i add '_wf' later)
% second argument is which spike you want to surround for the 1 min. raw clip
% (input is requested if there are more than a hundred timestamps)
% third input is "map", i.e. input the 4-vector of raw channels that correspond to the remapped data you're going for, 
% get raw data and waveforms/metrics out in struct; 
% option of filename to be 2nd output

[fname, pathname] = uigetfile('*.plx', 'Select a PLX file');
filename = fullfile(pathname, fname);
nexname = strcat(filename(1:end-3), 'nex');
sheetname = fname(1:end-3);
if length(sheetname) > 22
    sheetname = sheetname(1:22);
end

recording = struct('raw', [], 'neuron', []);
recording.neuron = struct('wfs', [], 'ts', [], 'meanwf', [], 'wfandstd', []);

% go get wfs ... 
[~, ~, ts, ~, wfs] = nex_wf(nexname, strcat(unit, '_wf'));
        if any(ts) == 1 || any(wfs) == 1
            ts = ts';
            wfs = wfs';
            meanwf = mean(wfs);
            stdwf = std(wfs);
            recording(4).neuron.wfs = wfs;
            recording(4).neuron.ts = ts;
            wfandstd = [meanwf+stdwf; meanwf; meanwf-stdwf]';
            almostinfo = recording(4).neuron.wfs(:, 8:256);
            unitinfo = [recording(4).neuron.ts almostinfo];
        elseif any(ts) == 0 && any(wfs) == 0
            disp('neither timestamps nor waveforms...')
        end
% then go get raw ...                    
for i = 1:numel(map)
    if i == 1
        if size(ts, 1) > 100
            disp (ts(end/2:end/2+50))
            disp ('last ts / 2 through last ts / 2 + 50 are shown; index for end / 2 + 50 follows...')
            disp (size(ts,1)/2+50)
            spknum = input('which INTEGER spike (index) should raw data center on? last idx displayed is end/2 + 50...');
        end
    end
    start = (recording(4).neuron.ts(spknum)-.5)*40e3; % must be < 40kHz * seconds of raw data to fit into a *.xls
    finish = (recording(4).neuron.ts(spknum)+.5)*40e3 - 1;
    [~, ~, ad] = plx_ad_span_v(filename, map(i), start, finish);
    recording(i).raw = ad;
end

%write everything to excel (after making rate hist as well...)
ratehist = binspikes(recording(4).neuron.ts, 1);
xlswrite(strcat(sheetname,'_', unit(4:6),'.xls'), ratehist, strcat(sheetname, 'ratehist'), 'B2');
xlswrite(strcat(sheetname,'_', unit(4:6),'.xls'), wfandstd, strcat(sheetname, 'wfandstd'), 'B2');
if any(ad) && numel(ad) > 1
    xlswrite(strcat(sheetname,'_', unit(4:6),'.xls'), [recording(1).raw recording(2).raw ...
    recording(3).raw recording(4).raw], strcat(sheetname, 'raw'), 'B2');
end
xlswrite(strcat(sheetname,'_', unit(4:6),'.xls'), unitinfo, strcat(sheetname, unit), 'B2');
plx_close(filename);
