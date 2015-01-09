function varargout = fiberAngle3(x,y,z, V, grid, corrfile, dbg)
%function [x,y,z,u,v,w, eigval] = fiberAngle3(x,y,z, V, grid, corrfile, dbg)
%x,y are by pixel coordinates (i.e., not just the endpoints)
%z is by slice coordinates
%V is the image volume, of size [length(y) length(x) length(z)]
%grid is a two element vector, with the size of the correlation volume as
%the first element and the step between correlation volumes as the second
%element
%(optional) corrfile is the name of the file to write the actual correlation
%matrices to
%(optional) dbg is true or false for debugging.  With debugging on, does not
%actually perform the correlations

%sort out parameters
if (nargin < 7),
  dbg = 0;
  if (nargin < 6),
    corrfile = [];
  end;
end;

%check whether we need to save the correlations
if (isempty(corrfile)),
  saveCorr = false;
else
  saveCorr = true;
end;

%same for the eigenvalues
if (nargout == 7),
    saveEig = true;
else
    saveEig = false;
end;

%check the data type of the image
if (~strcmp(class(V),'double')),
    convertToDouble = true;
    maxIntensity = double(max(V(:)));
else
    convertToDouble = false;
end;

%set up an inverse scale (world units to pixels)
if (length(z) > 1),
  scale = 1./abs([x(2)-x(1) y(2)-y(1) z(2)-z(1)]);
else
  scale = 1./abs([x(2)-x(1) y(2)-y(1)]);
  scale(3) = 0;
end;

%get the limits of the volume
lim = [min(x) min(y) min(z); max(x) max(y) max(z)];
%and the number of pixels
sz = size(V);

%can directly specify the position of the correlation volumes, if you
%want, using a cell array for grid
if (iscell(grid)),
    ctrx = grid{1};
    ctry = grid{2};
    ctrz = grid{3};

    ictr1 = round((ctrx-lim(1,1))*scale(1)) + 1;
    ictr2 = round((ctry-lim(1,2))*scale(2)) + 1;
    ictr3 = round((ctrz-lim(1,3))*scale(3)) + 1;

    igrid = round(grid{4}.*scale);
else
    %otherwise, do it the normal way
    %grid can also specify non-square regions if it's six elements long,
    %with the first 3 elements being the size in x y and z and the second
    %3 elements being the step in x y and z
    if (length(grid) == 2),
        grid = grid([1 1 1 2 2 2]);
    end;

    %get the grid size in pixels
    igrid = round(grid(1:3).*scale);
    %can't have a grid larger than the volume
    if (any(igrid(1:3) > sz(1:3))),
        warning('Grid sizes are larger than the volume.  Adjusting down.');
        k = find(igrid(1:3) > sz(1:3));
        igrid(k) = sz(k);
    end;
end;

%set up the coordinates for each correlation volume
%centers of the volumes in x y and z
%the "i" prefix is for an index/pixel coordinate
irgnctr = floor(igrid/2);

%pixel coordinates (irgn1) and real coordinates (rgnx)
irgn1 = -(irgnctr(1)-1):irgnctr(1);
rgnx = irgn1/scale(1);

irgn2 = -(irgnctr(2)-1):irgnctr(2);
rgny = irgn2/scale(2);

irgn3 = -(irgnctr(3)-1):irgnctr(3);
rgnz0 = irgn3/scale(3);

%generate the centers of the correlation volumes so that they're centered
%within the entire volume
if (~iscell(grid)),
    igridstep = round(grid(4:6).*scale);

    ictr1 = round(irgnctr(1)+1:igridstep(1):sz(2)-irgnctr(1));
    ictr1 = ictr1 + round((sz(2)-irgnctr(1)-ictr1(end))/2);
    ictr2 = round(irgnctr(2)+1:igridstep(2):sz(1)-irgnctr(2));
    ictr2 = ictr2 + round((sz(1)-irgnctr(2)-ictr2(end))/2);
    ictr3 = round(irgnctr(3)+1:igridstep(3):sz(3)-irgnctr(3));
    %make sure we have at least *some* points in the z direction
    if (isempty(ictr3)),
        ictr3 = ceil(sz(3)/2);
        k = find((ictr3 + irgn3 >= 1) & (ictr3 + irgn3 <= sz(3)));
        irgn3 = irgn3(k);
        rgnz0 = rgnz0(k);
    end;
    ictr3 = ictr3 + round((sz(3)-irgnctr(3)-ictr3(end))/2);

    %get the real coordinates
    ctrx = (ictr1-1)/scale(1) + lim(1,1);
    ctry = (ictr2-1)/scale(2) + lim(1,2);
    ctrz = (ictr3-1)/scale(3) + lim(1,3);
end;

%scale up the z dimension so that its scale matches the x and y
%dimensions.  Normally z is much coarser than x and y, so we'll have to
%interpolate, but it's important to make the FFT work reasonably
if (scale(3) < scale(1)),
    sc = scale(1);
    
    irgnctr(3) = floor(-rgnz0(1).*sc)+1;
    rgnz = (-(irgnctr(3)-1):irgnctr(3))/sc;
else
    rgnz = rgnz0;
end;

% h = (-3:3)';
% ih = repmat(irgnctr,[length(h) 1]) + repmat(h,[1 3]);
difford = 5;
hctr = ceil(length(rgnx)/2);

%set up the windowing transform (kaiser)
w1 = kaiser(length(rgny),0.5);
w2 = kaiser(length(rgnx),0.5);
w3 = kaiser(length(rgnz),0.5);
wind = repmat(w1,[1 length(w2) length(w3)]) .* ...
    repmat(w2',[length(w1) 1 length(w3)]) .* ...
    repmat(shiftdim(w3,-2),[length(w1) length(w2) 1]);

diff1 = [];
diff2 = [];

if (saveCorr),
    corrfid = fopen(corrfile,'wb');
    fwrite(corrfid,[length(ictr1) length(ictr2) length(ictr3)],'int8');
    fwrite(corrfid,[length(rgnx) length(rgny) length(rgnz)],'int8');
end;
    
timedWaitBar(0, 'Fiber angles...');
a = 1;
N = length(ictr1)*length(ictr2)*length(ictr3);
for i = 1:length(ictr1),
    for j = 1:length(ictr2),
        for k = 1:length(ictr3),
            V1 = V(ictr2(j)+irgn2, ictr1(i)+irgn1, ictr3(k)+irgn3);
            if (convertToDouble),
                V1 = double(V1)/maxIntensity;
            end;

            %check to make sure that there's a signal in this region
            if (~dbg & (std(V1(:)) ~= 0))
                %interpolate in z if necessary
                if (length(rgnz) > length(rgnz0)),
                    V1 = interp3(rgnx,rgny,rgnz0,V1, rgnx,rgny',rgnz);
                end;

                %apply the window
                V1 = V1.*wind;

                %do the 3D fft, multiply by the complex conjugate, and
                %inverse fft to get the autocorrelation
                V1f1 = fftn(flipdim(flipdim(flipdim(V1,1),2),3));
                V1f2 = fftn(V1);
                C = ifftn(V1f1.*V1f2);

                %this puts the peak in the center
                C = fftshift(C);
                C = real(C);            % just take the real part

                %get the hessian at the center of the correlation matrix
                H = hessian3(rgnx,rgny, rgnz, C, 'difference'); %,difford,diff1,diff2);
                
                %run through them all and display the x,y components of the
                %smallest eigenvectors
                hu = zeros(length(rgnx),length(rgny)); 
                hv = zeros(length(rgnx),length(rgny));
                ed = zeros(length(rgnx),length(rgny));
                for h1 = 1:length(rgny),
                    for h2 = 1:length(rgnx),
                        h3 = ceil(size(H,5)/2);
                        if (all(isfinite(flatten(H(:,:,h1,h2,h3))))),
                            %and take the eigenvalues
                            [ev,ed1] = eig(H(:,:,h1,h2,h3));
                            %minimum eigenvalue corresponds to the axis of least
                            %curvature
                            [ed(h1,h2),ind] = min(abs(diag(ed1)));

                            hu(h1,h2) = ev(1,ind);
                            hv(h1,h2) = ev(2,ind);
                        end;
                    end;
                end;
                        
                    
                H = H(:,:,hctr,hctr,hctr);

                %and take the eigenvalues
                [ev,ed] = eig(H);
                ed = diag(ed);
                %minimum eigenvalue corresponds to the axis of least
                %curvature
                [m,ind] = min(abs(ed));

                %save correlation matrix and/or eigenvalues
                if (saveCorr),
                    fwrite(corrfid,C,'double');
                    fwrite(corrfid,H,'double');
                end;
                if (saveEig),
                    eigval(:,j,i,k) = ed([ind 1:ind-1 ind+1:end]);
                end;

                %save the eigenvector
                u(j,i,k) = ev(1,ind);
                v(j,i,k) = ev(2,ind);
                w(j,i,k) = ev(3,ind);
            else
                u(j,i,k) = 0;
                v(j,i,k) = 0;
                w(j,i,k) = 0;
            end;
            if (~timedWaitBar(a/N)),
                if (saveCorr),
                    fclose(fid);
                end;
                varargout = {ctrx,ctry,ctrz, u,v,w};
                return;
            end;
            a = a+1;
        end;
    end;
end;

if (saveCorr),
    fclose(corrfid);
end;
if (dbg),
    varargout = cell(1,nargout);
else
    varargout = {ctrx,ctry,ctrz, u,v,w};
    if (saveEig),
        varargout(7) = {eigval};
    end;
end;


