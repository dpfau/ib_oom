function Ixz = mi(oom,p_z)
% Mutual information between a latent state Z_t-1 and the observation X_t.
% If no input is provided for p(z_t-1), the stationary distribution of the
% OOM is assumed.

p_xz = squeeze(sum(oom))'; % p(x|z)
if nargin == 1
%     [v ~] = eig(sum(oom,3));
%     p_z = v(:,1)/sum(v(:,1)); % p(z), stationary distribution
%   Taking eigenvectors of large matrices is cumbersomely slow, not to
%   mention inefficient since we only want the first one, so just calculate
%   the stationary distribution by brute force instead
    T = sum(oom,3); % transition matrix for latent states
    dz = 1;
    p_z = rand(size(oom,1));
    p_z = p_z/sum(p_z);
    while dz > 1e-12
        p_z1 = T*p_z;
        dz = sum(abs(p_z1-p_z))/length(p_z);
        p_z = p_z1;
    end
end
p_x = p_xz*p_z;

Ixz = entropy(p_x) - entropy(p_xz)*p_z;