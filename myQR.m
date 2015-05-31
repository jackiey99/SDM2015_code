function [Q,R] = myQR(U, X)
    % QR-decomposition like algorithm
    % Author: Liangyue Li (Arizona State University)
    [m, l] = size(U);
    t = size(X,2);
    
    R = zeros(t);
    Q = X - U*(U'*X);
    for j = 1:t
        R(j,j) = norm(Q(:,j));
        Q(:,j) = Q(:,j)/R(j,j);
        for k = j+1:t
            R(j,k) = Q(:,j)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(j,k)*Q(:,j);
        end
    end
    Q = [U,Q];
    I = eye(l);
    A = U'*X;
    B = zeros(t,l);
    R = [I,A;B,R];
end