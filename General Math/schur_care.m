function X = schur_care(A_mat,B_mat,Q_mat,R_mat)
% function X = schur_care(A_mat,B_mat,Q_mat,R_mat)
%
% Compute the solution to
%   A'*X + X*A - X*B*R^(-1)*B'*X + Q = 0
% using the Schur Algorithm.
%
% This code is based on
% Bini, D. A., Iannazzo, B., & Meini, B. (2012). Numerical Solution of 
%   Algebraic Riccati Equations. SIAM. Retrieved from 
%   http://books.google.com/books?id=Gi1jZAz3pncC&pg=PA102&dq=Numerical+Solution+of+Algebraic+Riccati+Equations+initialize+newton&hl=en&sa=X&ei=JV9aU_DaCZHgsAT17oLQAg&ved=0CC0Q6AEwAA#v=onepage&q=Numerical Solution of Algebraic Riccati Equations initialize newton&f=false
%
% Inputs
%   A_mat: system dynamics (n by n)
%   B_mat: force influence matrix (n by m)
%   Q_mat: state penalty (n by n, symmetric, positive semidefinite)
%   R_mat: input penalty (m by m, symmetric, positive definite)
% Outputs
%   X: solution to A'*X + X*A - X*B*R^(-1)*B'*X + Q = 0
%
% Written by:
%   David Allen
%   davidallen@vt.edu
%   25 April 2014
%
% |=======================================================================|
% | Version    Author      Date         Comments                          |
% | 1.0.0      D Allen     25 Apr 14    Original                          |
% |                                                                       |
% |=======================================================================|

% Get sizes
[nA1,nA2]=size(A_mat);
[nB,mB]=size(B_mat);
[nQ1,nQ2]=size(Q_mat);
[mR1,mR2]=size(R_mat);


% Check sizes
if nA1~=nA2
    error('A_mat must be a square matrix')
else
    n=nA1;
end

if nB~=n
    error('B_mat must be a n by m matrix')
else
    m=mB;
end

if nQ1~=nQ2
    error('Q_mat must be a square matrix')
elseif nQ1~=n
    error('Q_mat must be a n by n matrix')
end

if mR1~=mR2
    error('R_mat must be a square matrix')
elseif mR1~=m
    error('R_mat must be a m by m matrix')
elseif mR1 == 1
    disp('WARNING: schur_care may not provide good results for m=1.')
end

% Check for controllability
rankC=rank(ctrb(A_mat,B_mat));
if rankC~=n
    disp('WARNING: System not controllable, solution may not be valid.')
end


% Define matricies to agree with Bini et al
A=A_mat;
B=B_mat*(R_mat\(B_mat'));
C=Q_mat;

% This is the algorithm from Bini et al
H=[A,-B;-C,-A'];
[U,T]=schur(H,'real');% compute the schur form of H
e=ordeig(T);
[es,is]=sort(real(e),'ascend');
sel = zeros(2*n,1);
sel(is(1:n)) = 1;
Q= ordschur(U,T,sel); % Sorting the Schur form of H
X=Q(n+1:2*n,1:n)/Q(1:n,1:n);

end