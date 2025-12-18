%%find_helices:
% This script finds the 7 transmembrane helices in the input PDB structure
% required input: ssProt

% Find TM helices, GPCRs have 7 TM helices

function helices = findHelices(ssProt)
    isHelixBefore = false; % was the stretch before this a helix?
    helixStretch = false; % How long is this helix?
    helixCounter = true; % Which TM helix are we at?

    helixLength = [30 30 35 30 40 40 35]; % Helices have a typical length, please double check these values

    helices = zeros(7,2);  % 7 TM helices
    for letter = 1:(length(ssProt) - 2)
        n = ssProt(letter);
        np1 = ssProt(letter + 1);
        np2 = ssProt(letter + 2);

        isHelix = ((n == 'H') && (np1 == 'H') && (np2 == 'H')) ...
            && ~((n ~= 'H') && (np1 ~= 'H') && (np2 ~= 'H'));

        if isHelix % If current residue is a helix, add to helixStretch
            helixStretch = helixStretch + 1;
        end

        if helixStretch > helixLength(helixCounter) % If helix is stretching too much, loosen the constraints
            if (n ~= 'H') && (np1 ~= 'H')
                isHelix = false; % Helix seems to end!
            end
        end

        if isHelix && ~isHelixBefore % Helix started
            helices(helixCounter, 1) = letter;
        elseif ~isHelix && isHelixBefore % Helix ended
            helices(helixCounter, 2) = letter;
            helixCounter = helixCounter + 1; % Next helix pls!
            helixStretch = 0;
        end
        % update old helix index:
        isHelixBefore = isHelix;

        if helixCounter == 8 % we need a total of 7 helices
            break;
        end
    end

    assert(isempty(find(helices == 0, 1)), "Code had difficulty finding helices automatically! Please input helices variable manually and try again!");
end
