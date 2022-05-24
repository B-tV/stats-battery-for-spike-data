clear recording
[fname, pathname] = uigetfile('*.nex', 'Select a NEX file');
filename = fullfile(pathname, fname);
nex_info(filename);  %just to get the length right quick...
h = ['01'; '05'; '09'; '13'];
j = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'; 'j'; 'k'];
recording=struct('times',0);
for i = 1:length(h)
    for k = 1:length(j)
        %unit=strcat('SPK', h(i,:), j(k)); % for troubleshooting what
        %happened to some unit or another... COULD REQUIRE "TETSPK..." 
        [t ts] = nex_ts(filename,strcat('SPK', h(i,:), j(k)));
    if ts==0;
       disp('nothing recognized using SPK... trying TETSPK')
       [t ts] = nex_ts(filename,strcat('TETSPK', h(i,:), j(k)));
    end
    if ts~=0
        recording(numel(recording)+1) = struct('times', ts);
    end
    end
end

recording(1)=[];  % repeat (with other indices) as necessary to remove low-quality units
clear h i j k t ts


% save the first structure array as 'recording1', then run this script 
% on the next file, and make sure the newer data has the same # of units as recording1 
% before writing more data to recording1 i.e. only keep units you captured for all x# of file splits
% 
% for m = 1:numel(recording)
%     recording1(m).times(end+1:end+numel(recording(m).times))=recording(m).times;
% end

%%% CAUTION HERE IF YOU DON'T WANT AN XLS OF THE OVERALL RATE HIST's OUTPUT! (i.e. comment out if you don't want it) 
overallratehist = binspikes(recording, 1);
xlswrite(strcat(fname,'.xls'), overallratehist);
