function [Xout,Fout] = qrlinsolve(A,B)
%qrlinsolve Solve linear equations with (A*A')X=B using QR decomposition

%   Reference
%   [1] Guerci, Space-Time Adaptive Processing for Radar, 2003, pp173

[~,AR,PIVOTPERM] = qr(A',0);         % AR is upper triangular, see [1]
% Forward substitution after backward substitution
F = AR'\B(PIVOTPERM,:);
X = AR\F;
%X(PIVOTPERM,:) = X;
[~,PIVOTPERMSRT] = sort(PIVOTPERM);
Xout = X(PIVOTPERMSRT,:);
if nargout > 1
    Fout = F(PIVOTPERMSRT,:);
else
    Fout = F;
end

end