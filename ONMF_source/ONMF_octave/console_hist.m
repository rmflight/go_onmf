function [table] = console_hist(X)

ind = unique(X);
freq = [];
for i=ind
    freq = [freq sum(X == i)];
end

table = [ind; freq];