function alignment = alignSequences(sequences)
    seqCount = length(sequences);
    realSeq = zeros(seqCount, 1);

    for seqIndex = 1:seqCount
        seq = sequences{seqIndex};
        realSeq(seqIndex) = (~isempty(seq) && ~strcmp(seq, repmat('-', 1, length(seq))));
    end

    realSeqCount = sum(realSeq);
    realSeqIndices = find(realSeq);

    if realSeqCount < 1
        alignment = repmat('-', seqCount, 0);
    elseif realSeqCount <= 1
        seq = sequences{realSeqIndices(1)};
        alignment = repmat('-', seqCount, length(seq));
        alignment(realSeqIndices, :) = seq;
    elseif realSeqCount <= 2
        [~, rawAlignment] = nwalign(string(sequences{realSeqIndices(1)}), string(sequences{realSeqIndices(2)}));
        alignment = repmat('-', seqCount, size(rawAlignment, 2));
        alignment(realSeqIndices, :) = rawAlignment([1 3], :);
    else
        rawAlignment = multialign(arrayfun(@(x) sequences{x}, realSeqIndices, 'UniformOutput', false));
        alignment = repmat('-', seqCount, size(rawAlignment, 2));
        alignment(realSeqIndices, :) = rawAlignment;
    end
end
