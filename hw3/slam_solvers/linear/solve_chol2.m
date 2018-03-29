% SOLVE_CHOL2
% 16-831 Fall 2016 - *Stub* Provided
% Solves linear system using second Cholesky method
%
% Arguments: 
%     A     - A matrix from linear system that you generate yourself
%     b     - b vector from linear system that you generate yourself
%
% Returns:
%     x     - solution to the linear system, computed using the specified
%             version of the Cholesky decomposition
%     R     - R factor from the Cholesky decomposition
%
function [x, R] = solve_chol2(A, b)

S = A'*A;
r = symamd(S);

T = chol(S(r,r),'lower');
y = forward_sub(T, A(:,r)'*b);
x = back_sub(T', y);

R=T';
x(r)=x;


end
