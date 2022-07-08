%--------------------------------------------------------------------------
% LGR_Dmatrix.m
% determines approximate differentiation matrix for Legendre-based method
% with LGR nodes
%--------------------------------------------------------------------------
% D = LGR_Dmatrix(tau)
% tau: LGR nodes
%   D: differentiation matrix
%--------------------------------------------------------------------------
% Author: Lingwei Li, Ph.D Student, Northwestern Polytechnical University
% Date: 06/24/2022
%--------------------------------------------------------------------------
function D = LGR_Dmatrix(tau)
    % number of nodes
    N = length(tau)-1;

    % See Page 104 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    % Algorithms, Analysis and Applications, Springer Series in Compuational
    % Mathematics, 41, Springer, 2011. 
    % Uses the function: lepoly()
    % Original function: D = legslbdiff(n,x) located at
    % http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslbdiff.m
    % in fact, this code is found in GPOPS-I, gpopsCollocD.m
    n = N + 1;
    if n==0, D = [];                % null differentiation matrix
        return; 
    end   
    xxPlusEnd = tau; 
    M = length(xxPlusEnd);
    M1 = M+1;
    M2 = M*M;
    
    % Compute the barycentric weights
    Y = repmat(xxPlusEnd,1,M);
    Ydiff = Y - Y'+eye(M);
    
    WW = repmat(1./prod(Ydiff,2),1,M);
    D = WW./(WW'.*Ydiff);

    D(1:M1:M2) = 1-sum(D);
    D = -D';                %Diff matrix
%     D = D(1:end-1,:);       %Augment for LGR points
end