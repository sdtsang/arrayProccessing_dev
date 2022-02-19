%---------------------------------
function [pos,ang,N_elem,D_src] = parseInputsSteerVec(pos_in,ang_in)
eml_assert_no_varsize(1:nargin,pos_in,ang_in);
if size(pos_in,1) == 1
    pos = [zeros(1,size(pos_in,2));pos_in;zeros(1,size(pos_in,2))];
elseif size(pos_in,1) == 2
    pos = [zeros(1,size(pos_in,2));pos_in];
else
    pos = pos_in;
end
sigdatatypes.validate3DCartCoord(pos,'sensorcov','POS');

N_elem = size(pos,2);

if size(ang_in,1) == 1
    ang = [ang_in;zeros(1,size(ang_in,2))];
else
    ang = ang_in;
end
sigdatatypes.validateAzElAngle(ang,'sensorcov','ANG');

D_src = size(ang,2);

% if isscalar(ncov_in)
%     validateattributes(ncov_in,{'double'},{'scalar','real',...
%         'finite','nonnegative'},'sensorcov','NCOV');
%     ncov = ncov_in*eye(N_elem);
% elseif isrow(ncov_in)
%     validateattributes(ncov_in,{'double'},{'real',...
%         'finite','positive','size',[1 N_elem]},'sensorcov','NCOV');
%     ncov = diag(ncov_in);
% else
%     validateattributes(ncov_in,{'double'},{'finite','size',[N_elem N_elem]},...
%         'sensorcov','NCOV');
%     tol = 10*eps(max(abs(diag(ncov_in))));  % based on stats cholcov
%     cond = any(any(abs(ncov_in - ncov_in') > tol));
%     if cond
%         coder.internal.errorIf(cond,...
%              'phased:phased:notHermitian','NCOV');
%     end
%     % Check for positive definite
%     [~,nFlagPD] = chol(ncov_in);
%     cond = logical(nFlagPD);
%     if cond
%         coder.internal.errorIf(cond,...
%              'phased:phased:notPositiveDefinite','NCOV');
%     end
%     ncov = ncov_in;
% end
% 
% 
% 
% if isscalar(scov_in)
%     validateattributes(scov_in,{'double'},{'scalar','real',...
%         'finite','positive'},'sensorcov','SCOV');
%     scov = scov_in*eye(N_ang);
% elseif isrow(scov_in)
%     validateattributes(scov_in,{'double'},{'real',...
%         'finite','positive','size',[1 N_ang]},'sensorcov','SCOV');
%     scov = diag(scov_in);
% else
%     validateattributes(scov_in,{'double'},{'finite','size',[N_ang N_ang]},...
%         'sensorcov','SCOV');
%     tol = 10*eps(max(abs(diag(scov_in))));   % based on stats cholcov
%     cond = any(any(abs(scov_in - scov_in') > tol));
%     if cond
%         coder.internal.errorIf(cond,...
%              'phased:phased:notHermitian','SCOV');
%     end
% 
%     % Check for positive semi definite
%     [~,sED] = eig((scov_in+scov_in')/2);  % ensure Hermitian
%     sED = diag(sED);
%     tol = eps(max(abs(sED))); % based on stats cholcov
%     sED(abs(sED)<=tol)=0;
%     cond = any(sED<0);
%     if cond
%         coder.internal.errorIf(cond,...
%              'phased:phased:notPositiveSemiDefinite','SCOV');
%     end
% 
%     scov = scov_in;
end