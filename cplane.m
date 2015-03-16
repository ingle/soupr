function [finterp,xgrid,ygrid,mse] = cplane(xvec,yvec,xd,yd,fd,varargin)
% FINTERP = INTERPCART3D(XVEC, YVEC, XD, YD, FD,...)
% Inputs:
% XVEC is a vector of equispaced x-cooridinate values generated using linspace
% XD is a vector with x-coordinates of known data
% FD is vector of known data values for each (XD,YD)
% Note: NUMEL(FD), NUMEL(XD), NUMEL(YD)  must all be equal
% 
% Optional arguments:
% PENALTY 'laplacian' default, other options are 'gradient' and 'gradlap'
% SOLVER 'lsqr' default, other option is 'backslash'
% VERBOSITY (keep off if using crossval!)
%
% Date: Oct 1, 2013
% Author: Atul Ingle
%
% Acknowledgement: John D'Errico for his "gridfit" submission on Matlab FEX
% which inspired this version to use histc() for calculating interpolation 
% ratios quickly

if isempty(varargin)
    penalty = 'laplacian'; solver='lsqr'; lambda=0.01; verbosity = 0;
elseif numel(varargin)==1
    penalty = varargin{1}; solver='lsqr'; lambda=0.01; verbosity = 0;
elseif numel(varargin)==2
    penalty = varargin{1}; solver = varargin{2}; lambda=0.01; verbosity = 0;
elseif numel(varargin)==3
    penalty = varargin{1}; solver = varargin{2}; lambda=varargin{3}; verbosity = 0;
elseif numel(varargin)==4
    penalty = varargin{1}; solver = varargin{2}; lambda=varargin{3}; verbosity = varargin{4};
else
    error('varargin takes at most 4 arguments: cplane(..., penalty, solver, lambda, verbosity)');
end

N = numel(xd);
assert( numel(yd)==N && numel(fd)==N );
xd=xd(:); yd=yd(:); fd=fd(:); % recall: we want to fit f = f(x,y).
xvec = xvec(:); yvec = yvec(:);
Nx = numel(xvec); Ny = numel(yvec);

if verbosity
	disp(['fraction sampled=',num2str(N/(Nx*Ny)*100),'%']);
end
% determine which cell in the array each point lies in
[~,indx] = histc(xd,xvec);
[~,indy] = histc(yd,yvec);
indov = (indx==Nx); indx(indov) = indx(indov)-1;
indov = (indy==Ny); indy(indov) = indy(indov)-1;
ind = indy + Ny*(indx-1);

delx = diff(xvec); dely = diff(yvec);
% interpolation ratios
tx = (xd - xvec(indx))./delx(indx);
ty = (yd - yvec(indy))./dely(indy);

      
% bilinear interpolation
if verbosity
	disp('forming interpolation matrix...');
end
A = sparse(repmat((1:N)',1,4),[ind,ind+1,ind+Ny,ind+Ny+1], ...
[(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], ...
N,Nx*Ny);



if verbosity
	disp('forming penalty term matrix...');
end
mindel = min([ delx(1), dely(1)]); % assuming uniform spaced points
% mindel is for scaling purposes, otherwise 1/delx^2 and 1/dely^2 
% may become too large (because delx and dely are too small).
alljj = 1:Nx*Ny; [tmpy, tmpx] = ind2sub([Ny, Nx], alljj);

switch penalty
	case 'gradient'
		error('not implemented yet');

	case 'laplacian'
        B = sparse([1:Nx*Ny,alljj(tmpy-1>=1),alljj(tmpy+1<=Ny),alljj(tmpx-1>=1),alljj(tmpx+1<=Nx)],...
                [1:Nx*Ny,sub2ind([Ny, Nx], tmpy(tmpy-1>=1)-1, tmpx(tmpy-1>=1)),...
                sub2ind([Ny, Nx], tmpy(tmpy+1<=Ny)+1, tmpx(tmpy+1<=Ny)),...
                sub2ind([Ny, Nx], tmpy(tmpx-1>=1), tmpx(tmpx-1>=1)-1),...
                sub2ind([Ny, Nx], tmpy(tmpx+1<=Nx), tmpx(tmpx+1<=Nx)+1)],...
            [ones(1,Nx*Ny)*(-2)*((mindel/delx(1))^2 + (mindel/dely(1))^2),...
             ones(1,sum(tmpy-1>=1))*(mindel/dely(1))^2,...
             ones(1,sum(tmpy+1<=Ny))*(mindel/dely(1))^2,...
             ones(1,sum(tmpx-1>=1))*(mindel/delx(1))^2,...
             ones(1,sum(tmpx+1<=Nx))*(mindel/delx(1))^2,...
             ],...
            Nx*Ny, Nx*Ny);
end

if verbosity
	disp('solving system of equations...');
end
switch penalty
    case 'gradient'
		error('gradient penalty not implemented yet');
        strct = whos('A');
        disp(['size of A=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dx');
        disp(['size of Dx=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dy');
        disp(['size of Dy=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dy');
        disp(['size of Dy=', num2str(strct.bytes/1024/1024),' MB']);
        switch solver
            case 'lsqr'
                disp('using lsqr solver');
                finterp = lsqr((A'*A) + lambda*(Dx'*Dx + Dy'*Dy),(A'*fd), 1e-6, 10000);
            otherwise
                disp('using backslash solver');
                finterp = ((A'*A) + lambda*(Dx'*Dx + Dy'*Dy))\(A'*fd);
        end
    case 'gradlap'
		error('gradlap penalty not implemented yet.');
        strct = whos('A');
        disp(['size of A=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dx');
        disp(['size of Dx=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dy');
        disp(['size of Dy=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('Dy');
        disp(['size of Dy=', num2str(strct.bytes/1024/1024),' MB']);
        strct = whos('B');
        disp(['size of B=', num2str(strct.bytes/1024/1024),' MB']);
        switch solver
            case 'lsqr'
                disp('using lsqr solver');
                finterp = lsqr((A'*A) + lambda*(Dx'*Dx + Dy'*Dy)+ lambda*(B'*B),(A'*fd), 1e-6, 10000);
            otherwise
                disp('using backslash solver');
                finterp = ((A'*A) + lambda*(Dx'*Dx + Dy'*Dy)+lambda*(B'*B))\(A'*fd);
        end
    otherwise %laplacian
		if verbosity
			strct = whos('A');
			disp(['size of A=', num2str(strct.bytes/1024/1024),' MB']);
			strct = whos('B');
			disp(['size of B=', num2str(strct.bytes/1024/1024),' MB']);
		end
        switch solver
            case 'lsqr'
				if verbosity
					disp('using lsqr solver');
				end
                finterp = lsqr((A'*A)+lambda*(B'*B), (A'*fd), 1e-6, 1000);
            otherwise
				if verbosity
					disp('using backslash solver');
				end
                finterp = ((A'*A)+lambda*(B'*B))\(A'*fd);
        end

end

if nargout==3
    [xgrid,ygrid]=meshgrid(xvec,yvec);
elseif nargout==4
	[xgrid,ygrid]=meshgrid(xvec,yvec);
    mse = norm( A*finterp - fd, 2 );
end

finterp = reshape(finterp, [Ny, Nx]);



