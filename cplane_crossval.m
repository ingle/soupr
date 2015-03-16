% FUNCTION [zgrid, xgrid, ygrid, ocvlambda, lambdaVec, valVec] = ...
%    CPLANE_CROSSVAL( x, y, f, xvec, yvec )
% Inputs: (x,y) known data locations
%           f   known data point at location (x,y)
%        (xvec,yvec) query location(s), both monotonic increasing vectors
% Outputs: zgrid imputed values at grid locations
%        (xgrid,ygrid) use with surf()
%          ocvlambda   lambda to minimize OCV leave-out-one score
%          lambdaVec   vector of lambdas used for gridsearch
%          valVec      vector of OCV scores on the lambdaVec grid
% Generates the "best" cplane fit by repeatedly fitting with one data
% point left out and using the lambda that makes the smallest ocv score
% lambda is searched using a gridsearch.
% Uses laplacian+backslash solver options in cplane()

function [zgrid, xgrid, ygrid, ocvlambda, lambdaVec, valVec] = ...
    cplane_crossval( x, y, f, xvec, yvec )

lambdaVec = [logspace( -5, 0, 10 ), 10]; % change this to your favorite range

valVec = [];
penalty = 'laplacian';
solver = 'backslash';
    
lno = 0;
for lambda = lambdaVec
    lno = lno+1;
    fprintf('[%2d/%2d]: %0.2e\n', lno, numel(lambdaVec), lambda);
    val = 0;
    
    Nd = numel(x);
    for ii=1:Nd % we will skip the ii^th data point
        [zgrid, ~,~,~] = cplane(xvec, yvec, x([1:ii-1, ii+1:end]), ...
            y([1:ii-1, ii+1:end]), f([1:ii-1, ii+1:end]),...
            penalty, solver, lambda, 0);

        % now predict the value at fd(Nd) using bilinear interpolation
        % for notation see Wikipedia article on bilinear interpolation
        xidx = find( xvec < x(ii), 1, 'last');
        yidx = find( yvec < y(ii), 1, 'last');
        x1 = xvec(xidx); x2 = xvec(xidx+1);
        y1 = yvec(yidx); y2 = yvec(yidx+1);
        fq11 = zgrid(yidx, xidx  ); fq12 = zgrid(yidx+1, xidx  );
        fq21 = zgrid(yidx, xidx+1); fq22 = zgrid(yidx+1, xidx+1);

        fr1 = (x2 - x(ii))/(x2-x1)*fq11 + (x(ii)-x1)/(x2-x1)*fq21;
        fr2 = (x2 - x(ii))/(x2-x1)*fq12 + (x(ii)-x1)/(x2-x1)*fq22;

        fdpred = (y2-y(ii))/(y2-y1)*fr1 + (y(ii)-y1)/(y2-y1)*fr2;
        val = val + (fdpred - f(ii))^2;

    end

    val = val/Nd;

    valVec = [valVec; val];
end

ocvlambda = lambdaVec( valVec == min(valVec) );
if ocvlambda==lambdaVec(end) || ocvlambda==lambdaVec(1)
    warning('optimal lambda was at the endpoint of gridsearch');
end

[fgrid,xgrid,ygrid,~] = cplane(xvec, yvec, x, y, f, penalty, solver, ocvlambda);

end
