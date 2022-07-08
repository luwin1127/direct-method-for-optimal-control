%--------------------------------------------------------------------------
% LGR_nodes.m
% determines Lagrange-Gauss-Lobatto (LGL) nodes
%--------------------------------------------------------------------------
% tau = LGR_nodes(N)
%   N: number of nodes minus 1, should be an integer greater than 0
% tau: LGR nodes
%--------------------------------------------------------------------------
% Author: Lingwei Li, Ph.D Student, Northwestern Polytechnical University
% Date: 06/24/2022
%--------------------------------------------------------------------------
function tau = LGR_nodes(N)
    syms t
    % calculate node locations
    for k = 0:N
        tau(k+1,1) = -cos(pi*k/t); % symbolically to maintain precision
    end
    tau = double(subs(tau,'t',N));
end