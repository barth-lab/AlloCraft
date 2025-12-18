function [fid] = add2log(mydir, yourMsg, fileName)
%logFile Adds message to log file 
%
%% Usage
%   add2log(mydir,yourMsg)
%   add2log(mydir,yourMsg, fileName)
%
% * mydir is the directory where you want your log file to be saved
%
% * yourMsg is the message you want added to the log file. The message
% could be a string or it could be a cell array of strings. The log will
% print the current date and time beside the message
% 
% * fileName (Optional) is the name of the log file, defaults to log.txt

% Set default file name
if nargin<3
    fileName =  'log.txt';  
end
    
fid = fopen(fullfile(mydir, fileName), 'a');
if fid == -1
  error('Cannot open log file.');
end

if  iscell(yourMsg)
    for thisLine = 1:length(yourMsg)
        if thisLine==1
            fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg{1});
        else
            fprintf(fid, '%s\n', yourMsg{thisLine});
        end
    end
else   
    fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
end

fclose(fid);

end

