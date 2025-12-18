function standardizedName = standardizeProtonatedStateName(name)
    % Define a function to standardize amino acid names from protonated
    % states (ex.: HSD, HSE) to standard names (ex.:HIS)
    % 
    %% Example usage:
    % inputNames = {'HSD', 'HSE', 'HSP', 'GLH', 'ASH', 'LYN', 'ALA', 'CYS'};
    % standardizedNames = cell(size(inputNames));
    % 
    % for i = 1:length(inputNames)
    %     standardizedNames{i} = standardizeAminoAcid(inputNames{i});
    % end
    % 
    % % Display the results
    % disp('Original Names:');
    % disp(inputNames);
    % disp('Standardized Names:');
    % disp(standardizedNames);

%     name = deblank(name);
    % Check for trailing blank character
    hasTrailingBlank = false;
    if length(name) > 3 && name(end) == ' '
        hasTrailingBlank = true;
        name = strtrim(name); % Remove the trailing blank character for processing
    end
    
    % Create a map of protonated forms to their main amino acid names
    protonationMap = containers.Map(...
        {'HSD', 'HSE', 'HSP', 'GLH', 'ASH', 'LYN'},... % List of protonated forms
        {'HIS', 'HIS', 'HIS', 'GLU', 'ASP', 'LYS'});  % Corresponding standard names
    
    % Check if the given name is in the map
    if isKey(protonationMap, name)
        standardizedName = protonationMap(name);
    else
        % If not found in the map, return the name as is
        standardizedName = name;
    end
    
    % Add the trailing blank character back if it was present in the input
    if hasTrailingBlank
        standardizedName = [standardizedName, ' '];
    end
end