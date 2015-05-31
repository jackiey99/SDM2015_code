function score = original_kernel(c, A1, A2)
% 
% Calculate graph kernel using original equation
% Author: Liangyue Li (Arizona State University)

% Input: 
%   c: decay factor
%   A1: Adjacency matrix of the first graph
%   A2: Adjacency matrix of the second graph
%   score: graph kernel of A1 and A2
    n = size(A1,1);
    q1 = ones(n,1)/n;
    q = kron(q1,q1);
    p = q;
    
    score = q'* ((speye(n*n) - c*kron(A1,A2))\p);
end