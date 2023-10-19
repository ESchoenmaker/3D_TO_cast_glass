function varargout = eig3(A, Options)
% D = eig3(A)
% [V,D] = eig3(A)
% eig(..., Options)
%
% Compute in one shot the eigen-values of multiples (3 x 3) matrices
%
% INPUTS:
%   A: (3 x 3 x n) array
%   Options, structure with optional fields
%       orthoflag: scalar boolean, if TRUE ortgonalization eigenvectors in
%           the subspace of the eigenvector with multiplicity.
%           default [FALSE]
%       reltol: double scalar, relative tolerance to detect multiplicity of
%           eogen value, default [1e-7]
%       warnflag: scalar boolean, if TRUE verbose warning if there is ambiguity
%           of eigen valurs with multiplicity is detected, default [TRUE]     
% OUTPUTS:
%   D: (3 x n). EIG3 returns in D(:,k) three eigen-values of A(:,:,k)
%   V: (2 x 2 x n) array, eigenvectors of A
%       A(:,:,k)*V(:,:,k) = V(:,:,k)*diag(D(:,k)) for k=1, 2, ...
%
% See also: CardanRoots, eig2, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%       Original 27-May-2010
%       25-March-2021: implement eigen vectors
if size(A,1) ~= 3 || size(A,2) ~= 3
    error('A must be [3x3xn] array');
end
if nargin < 2
    Options = struct();
end
Aorg = A;
A = reshape(A, 9, []).';
P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,5)+A(:,9));
% Principal minors
M11 = A(:,5).*A(:,9) - A(:,8).*A(:,6);
M22 = A(:,9).*A(:,1) - A(:,3).*A(:,7);
M33 = A(:,1).*A(:,5) - A(:,4).*A(:,2);
P1 = (M11 + M22 + M33);
% Determinant
P0 = - A(:,1).*M11 ...
     + A(:,4).*(A(:,2).*A(:,9)-A(:,8).*A(:,3)) ...
     - A(:,7).*(A(:,2).*A(:,6)-A(:,5).*A(:,3));
D = CardanRoots(P3, P2, P1, P0).';
D = sort(D,1);
if nargout <= 1
    varargout = {D};
else
    orthoflag   = getoptions(Options, 'orthoflag',  false);
    reltol      = getoptions(Options, 'reltol',     1e-7);
    warnflag    = getoptions(Options, 'warnflag',   true);
    
    V = zeros(size(Aorg));
    % Different cases for groupping eigenvalues
    D2 = D(3,:);
    D2 = real(D2).^2+imag(D2).^2;
    tol2 = reltol^2*D2;
    dL = diff(D,1,1);
    dL2 = real(dL).^2+imag(dL).^2;
    caseid = [1 2]*(dL2>tol2);
    for i=0:3
        pkeep = caseid == i;
        if ~any(pkeep)
            continue
        end
        switch i
            case 0
                L = D(1,pkeep);
                m = 3;
            case 1
                L = D([1 2],pkeep);
                m = [1 2];
            case 2
                L = D([1 3],pkeep);
                m = [2 1];
            case 3
                L = D(:,pkeep);
                m = [1 1 1];
        end
        nL = size(L,1); % length(m)
        Akeep = reshape(Aorg(:,:,pkeep),9,[]);
        ni = size(Akeep,2);
        col0 = cumsum([0 m]);
        for k=1:nL
            Ak = Akeep;
            Ak([1 5 9],:) = Ak([1 5 9],:)-L(k,:);
            Ak3x3 = reshape(Ak,[3,3,ni]);
            % QR decomposition
            try
                [~, R, p] = MultipleQR(Ak3x3,0);
            catch
                error('EIGS3 requires https://www.mathworks.com/matlabcentral/fileexchange/68976-multipleqr package to be installed');
            end
            p = double(p);
            mk = m(k);
            rk = (3-mk);          
            
            if warnflag
                % Get the diagonal elements of R
                Rii = reshape(R,[9,ni]);
                Rii = Rii([1 5 9],:);
                if ~isreal(Rii)
                    % working with sqr(abs(R)) is faster than abs(R) for complex
                    abs2Rii = real(Rii).^2+imag(Rii).^2;
                    R11 = abs2Rii(1,:);
                    % Find the numerical ranks
                    b = [abs2Rii >= R11*(reltol^2);
                        false(1,ni)];
                else
                    absRii = abs(Rii);
                    R11 = absRii(1,:);
                    % Find the numerical ranks
                    b = [absRii >= R11*reltol;
                        false(1,ni)];
                end
                [~,r] = min(b,[],1);
                r = reshape(r-1,1,[]);
                r = min(r,3);
     
                if any(r+mk ~= 3)
                    warning('eig3:Degenerate', 'Some matrices has too close eigenvalues');
                    r(:) = rk;
                end
            else
                r = rk + zeros(1,ni);
            end
                       
            Z = reshape(R(rk:-1:1,:,:),[rk, 3*ni]);
            Z = reshape(Z(:,(rk+1:3)'+3*(0:ni-1)),[rk, mk, ni]);           
            R = R(rk:-1:1,rk:-1:1,:);
            nop = repmat(uint32((rk-1:-1:0)'),[1 ni]);
            % Back substitution
            X = backsubs_mex(R, Z, r, nop);     % (rk x mk x ni) 
            I = eye(mk);
            X = [-X; repmat(I,[1,1,ni])];       % (3 x mk x ni)
            X(p + 3*(0:mk*ni-1)) = X;
            if orthoflag
                X = MultipleQR(X,0);
            end
            col = col0(k)+(1:mk);
            V(:,col,pkeep) = X;
        end
    end
    V = V ./ sqrt(sum(V.*conj(V),1)); 
    
    varargout = {V, D};
end    
end % eig3
%%
function value = getoptions(options, name, defaultvalue)
% function value = getoptions(options, name, defaultvalue)
fields = fieldnames(options);
found = strcmpi(name,fields);
if any(found)
    value = options.(fields{found});
    if isempty(value)
        value = defaultvalue;
    end
else
    value = defaultvalue;
end
end % getoptions
