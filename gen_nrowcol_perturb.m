function A = gen_nrowcol_perturb(n, t)
    %% generate symmetric perturbation matrix
    % randomly select t nodes and wire them to 100 other random nodes
    % Author: Liangyue Li (Arizona State University)
    A = zeros(n,n);
    for i = 1:t
        row_pos = randi([1,n],1,1);
        nonzero = 100;
        nzpos = randi([1,n],1,nonzero);
        w = zeros(1, n);
        w(nzpos) = 1;
        A(:,row_pos) = w;
        A(row_pos,:) = w';
    end
    A = triu(A,1) + tril(A,-1);
    A = (A+A')/2;
    A = sparse(A);
end