function [u_dmd, Phi, omega] = dmd_func(X, t, dt, r, u0)

if nargin < 5
    u0 = X(:, 1);
end

% STEP 0: data split
X1 = X(:,1:end-1);
X2 = X(:,2:end);

% STEP 1: singular value decomposition (SVD)
% r = 10;
[U, S, V] = svd(X1, 'econ');

Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

% STEP 2: low-rank subspace matrix
%         (similarity transform, least-square fit matrix, low-rank subspace matrix)
Atilde = Ur'*X2*Vr*Sr^(-1);

% STEP 3: eigen decomposition
% W: eigen vectors
% D: eigen values
[W, D] = eig(Atilde);

% STEP 4: original space DMD mode
Phi = X2*Vr*Sr^(-1)*W;	% DMD modes

lambda = diag(D);       % eigen value
omega = log(lambda)/dt; % log of eigen value

% STEP 5: reconstruct the signal
b = pinv(Phi)*u0;  % pseudo-inverse initial conditions

u_modes = zeros(r,length(t));  % DMD reconstruction for every time point

for i = 1:length(t)
    u_modes(:,i) =(b.*exp(omega*(t(i))));
end

u_dmd = Phi*u_modes;   % DMD resconstruction with all modes

end