function X = getSparse(mat)

X = sparse(mat(:,1)+1, mat(:,3)+1, mat(:,2));
