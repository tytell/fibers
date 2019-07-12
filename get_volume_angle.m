function [C,a,b, R] = get_volume_angle(V, wind, varargin)

opt.method = 'hough';   % or 'hessian' or 'svd'
opt.threshold = 0.8;
opt = parsevarargin(opt,varargin, 3);

if (nargin == 1) || isempty(wind)
    %set up the windowing transform (kaiser)
    w1 = kaiser(size(V,1),0.5);
    w2 = kaiser(size(V,2),0.5);
    w3 = kaiser(size(V,3),0.5);
    wind = repmat(w1,[1 length(w2) length(w3)]) .* ...
        repmat(w2',[length(w1) 1 length(w3)]) .* ...
        repmat(shiftdim(w3,-2),[length(w1) length(w2) 1]);
end

V = V.*wind;

%do the 3D fft, multiply by the complex conjugate, and
%inverse fft to get the autocorrelation
V1f1 = fftn(flip(flip(flip(V,1),2),3));
V1f2 = fftn(V);
C = ifftn(V1f1.*V1f2);

%this puts the peak in the center
C = fftshift(C);
C = real(C);            % just take the real part

ind1 = 1:size(C,1);
ind1 = ind1 - (ind1(end)+ind1(1))/2;
ind2 = 1:size(C,2);
ind2 = ind2 - (ind2(end)+ind2(1))/2;
ind3 = 1:size(C,3);
ind3 = ind3 - (ind3(end)+ind3(1))/2;

if strcmp(opt.method,'hough')
    [a1,a2,a3] = ndgrid(ind1,ind2,ind3);
    
    dang = atan(1/ind1(end));
    dang = pi/2/ceil(pi/2/dang);
    
    %generate evenly spaced angles on the sphere
    ang1 = 0:dang:pi;
    ang2 = -pi/2:dang:pi/2;
    
    %with radii spaced out evenly
    r = ind1;
    
    %and generate the x,y,z coordinates that correspond to those angles
    [ang1,ang2,r] = ndgrid(ang1,ang2,r);
    i1 = r.*sin(ang1).*cos(ang2);
    i2 = r.*sin(ang1).*sin(ang2);
    i3 = r.*cos(ang1);
    
    %interpolate on to the spherical coordinate system
    C2 = interpn(a1,a2,a3, C, i1,i2,i3, 'linear',NaN);
    
    %and sum over the radii
    angvote = sum(C2,3);

    %look for the angles that got the most votes
    [maxvote,k] = max(angvote(:));
    [i,j] = ind2sub(size(ang1(:,:,1)),k);
    
    ang1vote = ang1(i,j,1);
    ang2vote = ang2(i,j,1);
    
    a(1,1) = sin(ang1vote)*cos(ang2vote);
    a(2,1) = sin(ang1vote)*sin(ang2vote);
    a(3,1) = cos(ang1vote);

    R(1) = maxvote / sum(C2(:));
    
    %now we'll sum over the cylinder oriented along the current primary
    %axis
    
    %first make the angles
    ang1 = -pi/2:dang:pi/2;     % angle in the cylinder
    r = ind1;                   % radius from the central axis
    d = ind1;                   % distance along the central axis
    
    [ang1,r,d] = ndgrid(ang1,r,d);
    
    %get the vectors perpendicular to a
    [~,k] = max(abs(a));
    k = mod(k,3) + 1;
    m = zeros(3,1);
    m(k) = 1;
    
    %remove the component parallel to a
    m = m - dot(a,m)*a;
    m = m/norm(m);
    %and get the last perpedicular component
    n = cross(a,m);

    %now we'll construct coordinates in an a, m, n coordinate system, then 
    %project back on to the a1,a2,a3 coordinates
    mdist = r.*cos(ang1);
    ndist = r.*sin(ang1);
    
    i1 = a(1)*d + m(1)*mdist + n(1)*ndist;
    i2 = a(2)*d + m(2)*mdist + n(2)*ndist;
    i3 = a(3)*d + m(3)*mdist + n(3)*ndist;
    
    %interpolate on to the spherical coordinate system
    C3 = interpn(a1,a2,a3, C, i1,i2,i3, 'linear',NaN);
    
    %and sum over the radii and over the distance along the central axis
    angvote = nansum(nansum(C3,3),2);
    angn = sum(sum(isfinite(C3),3),2);
    angn(angn < 5) = NaN;
    
    [maxvote,k] = max(angvote(:)./angn(:));
    
    angvote = ang1(k,1,1);
    b(1,1) = m(1)*cos(angvote) + n(1)*sin(angvote);
    b(2,1) = m(2)*cos(angvote) + n(2)*sin(angvote);
    b(3,1) = m(3)*cos(angvote) + n(3)*sin(angvote);
    
    R(2) = maxvote / nansum(C3(angn >= 5));
elseif strcmp(opt.method,'svd')
    ctr1 = (size(C,1)+1)/2;
    ctr2 = (size(C,2)+1)/2;
    ctr3 = (size(C,3)+1)/2;
    
    if (ceil(ctr1) ~= ctr1)
        ctr1 = floor(ctr1) + [0; 1];
    end
    if (ceil(ctr2) ~= ctr2)
        ctr2 = floor(ctr2) + [0; 1];
    end
    if (ceil(ctr3) ~= ctr3)
        ctr3 = floor(ctr3) + [0; 1];
    end
    
    mid = mean(flatten(C(ctr1,ctr2,ctr3)));
    
    C = C ./ mid;
    
    [a1,a2,a3] = ndgrid(ind1,ind2,ind3);

    isthresh = C > opt.threshold;
    
    a1 = a1(isthresh);
    a2 = a2(isthresh);
    a3 = a3(isthresh);
    
    [U,S,V] = svd([a1(:) a2(:) a3(:)]);
    R = diag(S);
    
    a = V(:,1);
    b = V(:,2);
elseif strcmp(opt.method,'hessian')
    %get the hessian at the center of the correlation matrix
    H = hessian3(ind2,ind1,ind3, C, 'difference'); %,difford,diff1,diff2);
    hctr = ceil(length(ind1)/2);
    
    %run through them all and display the x,y components of the
    %smallest eigenvectors
    %                 hu = zeros(length(rgnx),length(rgny));
    %                 hv = zeros(length(rgnx),length(rgny));
    %                 ed = zeros(length(rgnx),length(rgny));
    %                 for h1 = 1:length(rgny),
    %                     for h2 = 1:length(rgnx),
    %                         h3 = ceil(size(H,5)/2);
    %                         if (all(isfinite(flatten(H(:,:,h1,h2,h3))))),
    %                             %and take the eigenvalues
    %                             [ev,ed1] = eig(H(:,:,h1,h2,h3));
    %                             %minimum eigenvalue corresponds to the axis of least
    %                             %curvature
    %                             [ed(h1,h2),ind] = min(abs(diag(ed1)));
    %
    %                             hu(h1,h2) = ev(1,ind);
    %                             hv(h1,h2) = ev(2,ind);
    %                         end;
    %                     end;
    %                 end;
    
    
    H = H(:,:,hctr,hctr,hctr);
    
    %and take the eigenvalues
    [ev,ed] = eig(H);
    ed = diag(ed);
    %minimum eigenvalue corresponds to the axis of least
    %curvature
    [~,ord] = sort(abs(ed));
    
    a = ev([2 1 3],ord(1));
    b = ev([2 1 3],ord(2));
end

