function [allresults, FILE, params_updated] = statsbatteryIIbatch(pathname,fnames,params,includeinj,D1ag,check)

% [allresults, FILE, params_updated] = statsbatteryIIbatch(pathname,fnames,params,includeinj,D1ag)
% ALSO outputs "newdata.xls" to the matlab path where yo shiznit be stowed
% pathname = a single character array, the initial part of the filename ( = pathname + fname)
% fnames = filenames specifically listed in a cell array (latter part of the filename)
% params = column 1-4 specifying the inputs for statsbatteryII; overloaded by KBD inputs through "check"
%   1st: delay onset (i.e. duration after beginning of recording to start pre period)
%   2nd: time of drug injection 
%   3rd: pre offset (i.e. duration before drug injection to end pre period, i.e. subtracted from 2nd)
%   4th: drug onset delay (i.e. duration to wait before starting drug period
%   should be as many rows as files; can be 1 row of "default parameters" row if unknown
% includeinj = input 1 to include all data immediately before/after injections; 
%   0 to ignore a default of about half a minute on either side of injection (or inform this using params)
% the following only matter WHEN CHECKING NEXFile FOR EVENT INFO:
% D1ag (ONLY matters when including inj) = input 1 to allow pre offset to end pre period at 2370; 
%   input 0 to allow pre-offset to default the pre period end to 4770, i.e. for l-DOPA recordings
% check = input 1 to check the event info from the file and to specify injection times as desired

if size(params, 1) == 1 % i.e. default parameters for any file without KBDevent information when params are unspecified
    params=repmat(params,size(fnames,1),1);  
end

allresults = zeros(size(fnames,1),56); % or a table?? ...just happens to be 56 columns (from ? metrics) 
FILE(size(fnames,1)) = struct('name',[],'KBD1',[],'KBD2',[],'KBD3',[],'KBD4',[],'KBD5',[],'KBD6',[],'KBD7',[],'KBD8',[]); 
% append as many KBD# events as you'd like and include them below ...
% others may be appended from within the script and indexed after these (e.g. 'saline', etc.)

for L=1:size(fnames,1) % # of filenames you have to check (should be same as number of files in that path!)
    clear recording
    filename = strcat(pathname,fnames{L},'.nex');
    %filename = fullfile(pathname, fnames(L));

    % for use in labeling results in various places (excel, structs with certain info, etc.)
    sheetname = filename; 
    sheetname(end-3:end) = [];
    sheetname(1:length(pathname)) = [];
    if length(sheetname) > 31 % because excel can't handle more than 31 character sheet names...
        sheetname(end-(end-31):end) = [];
    end
    
    nex_info(filename);  %just to get the length right quick...
% for labeling units to retrieve
    h = ['01'; '05'; '09'; '13'];
    j = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'; 'j'; 'k'];
    recording=struct('times',0);
    for i = 1:length(h)
        for k = 1:length(j)
        % unit=strcat('SPK', h(i,:), j(k)); 
        % for troubleshooting what happened to some unit or another... COULD REQUIRE "TETSPK..." 
        [~, ts] = nex_ts(filename,strcat('SPK', h(i,:), j(k)));
        if ts==0;
            disp('nothing recognized using SPK... trying TETSPK...')
            [~, ts] = nex_ts(filename,strcat('TETSPK', h(i,:), j(k)));
        end
        if ts==0;
            [~, ts] = nex_ts(filename,strcat('SPK', h(i,:), j(k),'_wfU'));
        end
        if ts~=0;
            recording(numel(recording)+1) = struct('times', ts);
        end               
        end
    end
    
% IF YOU WANT EVENT INFO FROM NEX FILE ...
    if check == 1
        NexFile = readNexFile(filename); 
        
% ... reorganize KBD events into something useful!
        if isfield(NexFile,'events')
            disp(sheetname)
            for i=1:length(NexFile.events)
                a = strcmpi(NexFile.events{i,1}.name,'KBD1');
                b = strcmpi(NexFile.events{i,1}.name,'KBD2');
                c = strcmpi(NexFile.events{i,1}.name,'KBD3');
                d = strcmpi(NexFile.events{i,1}.name,'KBD4');
                e = strcmpi(NexFile.events{i,1}.name,'KBD5');
                f = strcmpi(NexFile.events{i,1}.name,'KBD6');
                g = strcmpi(NexFile.events{i,1}.name,'KBD7');
                h = strcmpi(NexFile.events{i,1}.name,'KBD8');
                if a==1
                    FILE(L).KBD1 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD1')
                    disp(FILE(L).KBD1)
                elseif b==1
                    FILE(L).KBD2 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname; % redundant, but necessary in case there was only one event (not knowing which)
                    disp('KBD2')
                    disp(FILE(L).KBD2)
                elseif c==1
                    FILE(L).KBD3 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD3')
                    disp(FILE(L).KBD3)                                        
                elseif d==1
                    FILE(L).KBD4 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD4')
                    disp(FILE(L).KBD4)
                elseif e==1
                    FILE(L).KBD5 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD5')
                    disp(FILE(L).KBD5)
                elseif f==1
                    FILE(L).KBD6 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD6')
                    disp(FILE(L).KBD6)
                elseif g==1
                    FILE(L).KBD7 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD7')
                    disp(FILE(L).KBD7)
                elseif h==1
                    FILE(L).KBD8 = NexFile.events{i,1}.timestamps;
                    FILE(L).name = sheetname;
                    disp('KBD8')
                    disp(FILE(L).KBD8)
                end                   
            end
        else
            disp('NO NexFile EVENT info! ...using default params')
            FILE(L).NoEvtInfo = sheetname;
        end
    % in case there was no KBD3
        if isempty(FILE(L).KBD3)
            disp(filename)
            warning('no KBD3...inj?')
        end
   
% NOW go through the events and pick which ones are to be used as injections    
    % ignore saline injection when params ignore injections + have fudge factor for ~a-minute-later injection artifacts 
        if includeinj == 0
            if isempty(FILE(L).KBD1) && isempty(FILE(L).KBD2) && isempty(FILE(L).KBD3) && isempty (FILE(L).KBD4) ...
                    && isempty(FILE(L).KBD5) && isempty(FILE(L).KBD6) && isempty(FILE(L).KBD7) && isempty(FILE(L).KBD8)
                disp ('NO EVT info; using default params')
              % FILE(L).saline = params(L,1); 
            end
            x = input('what time should be used as SALINE injection? ("FILE(?).KBD?(?)" or timestamp) = ');
            y = input('ADD how many seconds to get to PREdrug ONset? (in seconds, NOT timestamp NOR field(x)) = ');
            % should be ~65 sec after injection in order to keep consistent with drug injection...
            params(L,1) = x+y;
            FILE(L).saline = params(L,1);
        end
    % drug injection obviously necessary regardless of how much around it is ignored 
        x = input('what time should be used as DRUG injection? ("FILE(?).KBD?(?)" or timestamp) = ');  
        params(L,2) = x;
        FILE(L).drug = x;
        if D1ag == 0
            if includeinj == 0
                y = input('SUBTRACT how many seconds to get to PREdrug OFFset? (in seconds, NOT timestamp NOR field(x)) = ');
                z = input('ADD how many seconds to get to POSTdrug ONset? (in seconds, NOT timestamp NOR field(x)) = ');
                params(L,3) = params(L,3)+y;
                params(L,4) = params(L,4)+z;
            else
                params(L,3) = params(L,3)+x-4770;   
            end
        else
            if includeinj == 0
                y = input('SUBTRACT how many seconds to get to PREdrug OFFset? (in seconds, NOT timestamp NOR field(x)) = ');
                z = input('ADD how many seconds to get to POSTdrug ONset? (in seconds, NOT timestamp NOR field(x)) = ');
                params(L,3) = params(L,3)+y;
                %if the duration after subtracting your pre-offset is more
                %than the effect period, then fix that!
                if params(L,2)-params(L,1) > 2400 
                    params(L,3) = (params(L,2)-params(L,1))-2400;
                end
                params(L,4) = params(L,4)+z;
            else
                params(L,3) = params(L,3)+x-2370;
            end
        end
    end


    recording(1)=[];  % repeat (with other indices) as necessary to remove low-quality units
    [results]= statsbatteryII(recording, params(L,1),params(L,2),params(L,3),params(L,4));
    %params(L,:) = [params(L,1) params(L,2) params(L,3) params(L,4)];  % should be pointless i think...
    
    % write to excel
    writetable(results,'newdata.xls', 'Sheet', sheetname);
    
    %and finally output for further whatever...
    if L==1
        allresults = table2array(results);
    else
        allresults = vertcat(allresults,table2array(results)); % unavoidable array growth without pre-specifying # of units
    end  
end
clear h i j k l t ts
params_updated = params; 

%%BREEEAHG!
disp('double fin!');
