function col = myhistc(x, t)
% col = myhistc(x, t)
% as [~, col] = histc(x,t) but with specific col value for x == max(t)
% that does not take the last index, meaning the last bin contains also the
% right bound

[~, col] = histc(x,t); %#ok
n = length(t);
b = col == n;
if any(b)
    tmax = t(end);
    for colmax = n:-1:1
        if t(colmax) ~= tmax
            break
        end
    end
    col(b) = colmax;
end

end % myhistc
