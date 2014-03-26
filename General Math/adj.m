function adjA=adj(A)
% function adjA=adj(A)
%
% This function finds the adjoint of matrix A
%
% Inputs:
%   A: The matrix whose adjoint is to be calculated
% Outputs:
%   adjA: the adjoint of A
%
% Written by:
%   David Allen
%   davidallen@vt.edu
%   25 March 2014

SizeA=size(A);% Get the size of A
C=zeros(SizeA);% Initialize C

for i=1:SizeA(1)
    for j=1:SizeA(2)
        tempA=A;
        tempA(i,:)=[];% Delete ith row
        tempA(:,j)=[];% Delete jth column
               
        C(i,j)=(-1)^(i+j)*det(tempA);% (i,j) element of C
    end
end

adjA=C';
end