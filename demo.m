clear all;
load as50days

A1 = graphs{1}.A;
A2 = graphs{2}.A;

n1 = size(A1, 1);
n2 = size(A2, 1);

% initial top r eig decomposition of A
[U1, Lam1] = eigs(A1, 20);
[U2, Lam2] = eigs(A2, 20);

% top dr eig decomposition on delta A (i.e., the update matrix)
dr = 20;

% decay parameter
c = 0.00001;

% start and end probability
q1 = ones(n1,1)/n1;
p1 = q1;
q2 = ones(n2,1)/n2;
p2 = q2;

% number of timestamps
num_seq = 10;

for t = 1:num_seq
    fprintf('Update step: %d \n', t); 
    
    % generate perturbation matrix
    delta_A1 = gen_nrowcol_perturb(n1, 20);
    delta_A2 = gen_nrowcol_perturb(n2, 20);
    
    A1 = A1 + delta_A1;
    A1 = (A1 + A1')/2;
    A2 = A2 + delta_A2;
    A2 = (A2 + A2')/2;
    
    % update the low rank approximation (ref. Alg 1 and Alg 2 in SDM paper)
    [U1, Lam1] = updateEigen(U1, Lam1, delta_A1, dr);
    [U2, Lam2] = updateEigen(U2, Lam2, delta_A2, dr);

    Lam12 = kron(diag(Lam1)', diag(Lam2)'); 
    Lam = Lam12./(1- c*Lam12);
    
    L = kron(q1'*U1, q2'*U2);
    R = kron(U1'*p1, U2'*p2);
    score = (q1'*p1)*(q2'*p2) + c*(L.*Lam)*R;
    format shortE
    fprintf('Graph kernel score: %e \n', score);
end