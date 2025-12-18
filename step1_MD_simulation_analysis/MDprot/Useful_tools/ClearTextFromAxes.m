function ClearTextFromAxes(handleToAxes)
%% Clears text on a figure 
try
    handlesToChildObjectsInAxes = findobj(handleToAxes, 'Type', 'text');
    if ~isempty(handlesToChildObjectsInAxes)
        delete(handlesToChildObjectsInAxes);
    end
catch ME
    errorMessage = sprintf('Error in program %s, function %s(), at line %d.\n\nError Message:\n%s', ...
        mfilename, ME.stack(1).name, ME.stack(1).line, ME.message);
    uiwait(warndlg(errorMessage));
end
return; % from ClearTextFromAxes
end