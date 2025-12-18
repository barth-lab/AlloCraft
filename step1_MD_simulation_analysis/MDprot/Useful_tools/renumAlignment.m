function [arrayRenum] = renumAlignment(array2Transform,pdb1,pdb2)
%renumAlignment Renumbers array2Transform from pdb1 to pdb2 numbering after
%sequence alignment
%
%% Usage:
% renumAlignment(array2Transform,pdb1,pdb2)

seqIn = pdb2seq(pdb2);
seqAc = pdb2seq(pdb1);

[score,alignment] = nwalign(seqIn{1},seqAc{1});

inaNdx = zeros(length(alignment),1);
actNdx = zeros(length(alignment),1);
counter = 1;
counterIna = 1;
for i = 1:length(alignment)
    if ~strcmp(alignment(3,i),'-')
    actNdx(i) = counter;
    counter = counter + 1;
    end

    if ~strcmp(alignment(1,i),'-')
    inaNdx(i) = counterIna;
    counterIna = counterIna + 1;
    end
end

 arrayRenum = arrayfun(@(x) inaNdx(actNdx == x),array2Transform);
end