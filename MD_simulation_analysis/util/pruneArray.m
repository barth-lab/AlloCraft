function [newdata, row, col] = pruneArray(data)
    [r, c] = find(data);
    row = unique(r);
    col = unique(c);
    newdata = data(row, col);
end
