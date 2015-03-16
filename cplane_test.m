xvec = linspace(-2,2,50);
yvec = linspace(-2,2,50);
sliz = linspace(-3,3,100);
% consider an ellipse: f(x,y) = 1 for x^2/1^2 + y^2/1^2 + z^2/2^2 <= 1
% and f(x,y) = 0 outside the ellipse.
P = 12;
sig = 0.2;

% set up radial sampling lines on each C-plane
R = linspace(0,2,40);
THETA = (0:2*P-1)*pi/P;
[Rm, Tm] = meshgrid(R, THETA);
Rm = Rm(:); Tm = Tm(:);

%------------------- estimate ocvlambda --------------------------------
zslice=0; % consider the largest circular crosssection of the ellipsoid
xd = zeros(size(Rm)); yd = xd; fd = xd;
flag=0;
for ii = 1:numel(Rm)
    xc = Rm(ii)*cos(Tm(ii));
    yc = Rm(ii)*sin(Tm(ii));
    if xc^2/1 + yc^2/1 + zslice^2/2^2 <=1
        fc=1;
    else
        fc=0;
    end
    if ( xc<max(xvec) && xc>min(xvec) && yc<max(yvec) && yc>min(yvec) )
        flag=flag+1;
        xd(flag)=xc;
        yd(flag)=yc;
        fd(flag)=fc;
    end
end
xd = xd(1:flag); yd = yd(1:flag); fd = fd(1:flag);

fdno = fd + sig*randn(size(fd));
[fgrid, xgrid, ygrid, ocvlambda, lambdaVec, valVec] = cplane_crossval( xd, yd, fdno, xvec, yvec );


%------------------------ reconstruct all cplanes with ocvlambda smoothing parameter -----------------

cplaneStack = zeros( size(fgrid,1), size(fgrid,2), numel(sliz) ); % holds the 3D stack of C-planes
for zslice = 1:numel(sliz)

    xd = zeros(size(Rm)); yd = xd; fd = xd;
    flag=0;
    for ii = 1:numel(Rm)
        xc = Rm(ii)*cos(Tm(ii));
        yc = Rm(ii)*sin(Tm(ii));
        if xc^2/1 + yc^2/1 + sliz(zslice)^2/2^2 <=1
            fc=1;
        else
            fc=0;
        end
        if ( xc<max(xvec) && xc>min(xvec) && yc<max(yvec) && yc>min(yvec) )
            flag=flag+1;
            xd(flag)=xc;
            yd(flag)=yc;
            fd(flag)=fc;
        end
    end
    xd = xd(1:flag); yd = yd(1:flag); fd = fd(1:flag);

    fdn = fd + sig*randn(size(fd));
     
    [finterp, xgrid, ygrid,mse] = cplane( xvec, yvec, xd, yd, fdn, 'laplacian', 'backslash', ocvlambda );
    cplaneStack(:,:,zslice) = finterp;

end

figure(1); stem3( xd, yd, fdno ); title('noisy data points on a c-plane');
figure(2); surf(fgrid, 'edgecolor', 'none'); title('surface fit to the c-plane of figure1');
% Uncomment the next line if you have slidingviewer.m from Matlab FileExchange:
% figure(3); slidingviewer(cplaneStack); colormap hot
