function [U, Lam] = updateEigen(U0, Lam0, delta_A, dr)
%{
 Goal : Update eigen decomposition
 Author: Liangyue Li (Arizona State University)

 Input:
   U0: eigen vectors of A (old matrix)
   Lam0: eigen values of A
   delta_A: perturbation
   dr: top dr eigen decomposition of delta_A

 Output:
   U: eigen vectors of new A
   Lam: eigen values of new A
%}
    % low-rank structure of delta_A
    
    [X, Y] = eigs(delta_A, dr);
    
    % perform qr on [U0, X], if don't care about speed, use matlab's qr
    [Q, R] = myQR(U0, X);
    %[Q, R] = qr([U0, X], 0);
    
    M = diag([diag(Lam0);diag(Y)]);
    Z = R*M*R';
    Z = (Z+Z')/2;
    
    r = rank(Z);
    if r < size(Z, 1)
        [V, Lam] = eigs(Z, r);
    else
        [V, Lam] = eig(Z);
    end
    
    U = Q*V;
    [~,ind] = sort(abs(diag(Lam)),'descend');
    top_r = min(size(U0,2),length(ind));
    ind = ind(1:top_r);
    Lam = Lam(ind,ind);
    U = U(:,ind);
end