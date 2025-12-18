function [BW, BW_cell] = pdb2BW(residues,helices)
%pdb2BW Changes PDB numbering of a GPCR to Ballesteros-Weinstein numbering
%   Detailed explanation goes here
% helices takes the form:
%   StartPDB  ENDPDB     StartBW   End BW
%   31.0000   62.0000    1.2900    1.6000
%
% helices for dopamine D2: 
% helices = [31 62 1.29 1.60;  67 98 2.37 2.68; 103 138 3.21 3.56; 144 173 4.34 4.63; ...
%    186 226 5.35 5.75; 362 399 6.24 6.61; 404 429 7.31 7.56];

BW = zeros(length(residues),1);
BW_cell = cell(length(residues),1);
loops = ['Nter';'ICL1';'ECL1';'ICL2';'ECL2';'ICL3';'ECL3';'Cter'];
formatSpec = '%.2f'; % Specify precision to get BW numbering

for resi = 1:length(residues)
    isTM = 0; % Is it TM? Or perhaps loops
    
    for TMH = 1:7 % 7 TM helices
       if  residues(resi) >= helices(TMH,1) && residues(resi) <= helices(TMH,2) % Is it in this helix?
           BW(resi) = helices(TMH,3) +(residues(resi)-helices(TMH,1))/100;
           BW_cell{resi} = num2str(BW(resi),formatSpec);
           isTM = 1;
       end
    end
    
    if isTM == 0 % It's not TM! prolly loops then
         for TMH = 1:7 % 7 TM helices
             
            if  residues(resi) < helices(TMH,1) 
                BW_cell{resi} = loops(TMH,:); % Give proper loop position
                break % we found it! stop the loop
            elseif TMH == 7 && residues(resi) > helices(TMH,2) % For the last one, check if it's Cterm
                BW_cell{resi} = loops(TMH+1,:);
            end 

         end
    end
end
end

