function H = tl_hilbert(N)
%
% H = tl_hilbert(N)
%
% Generate NxN matrix H representing the Hilbert transform for an Nx1
% vector.  Because hilbert() is linear, the H matrix is also the
% tangent-linear operator: tl_uH = H*tl_u. Furthermore, can use H' for
% adjoint coding
%
% The algorithm is brute force: compute N hilbert transforms on delta
% functions.
%
% TEST: compare hilbert() vs. tl_hilbert() for a random Nx1 vector:
% >> N=100;
% >> r=rand(N,1);
% >> plot(hilbert(r))
% >> hold on
% >> plot(tl_hilbert(N)*r)
%

H = hilbert(eye(N));
