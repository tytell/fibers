function fiberAngle3(h5imagefile,h5outfile, igrid, varargin)
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

opt.savecorr = false;
opt.saveeig = false;
opt.debug = false;
opt.method = 'hough';
opt = parsevarargin(opt,varargin,3);

debug = opt.debug;
issavecorr = opt.savecorr;
issaveeig = opt.saveeig;

info = h5info(h5imagefile,'/Image');
sz = info.Dataspace.Size;

attrnames = {info.Attributes.Name};
if ismember('Scale',attrnames)
    scale = h5readatt(h5imagefile,'/Image','Scale');
else
    scale = 1;
end
if ismember('Maximum',attrnames)
    imagemax = h5readatt(h5imagefile,'/Image','Maximum');
    imagemax = double(imagemax);
else
    imagemax = 1;
end

if numel(scale) == 1
    scale = scale*[1 1 1];
end

%note that x axis is the first dimension in the data set, not the second
x = ((1:sz(1)) - (sz(1)+1)/2) * scale(1);
y = ((1:sz(2)) - (sz(2)+1)/2) * scale(2);
z = ((1:sz(3)) - (sz(3)+1)/2) * scale(3);

%get the limits of the volume
lim = [min(x) min(y) min(z); max(x) max(y) max(z)];

%can directly specify the position of the correlation volumes, if you
%want, using a cell array for grid
if (iscell(igrid)),
    ictr1 = igrid{1};
    ictr2 = igrid{2};
    ictr3 = igrid{3};

    ctrx = round((ictr1-1)*scale(1)) + lim(1);
    ctry = round((ictr2-1)*scale(2)) + lim(2);
    ctrz = round((ictr3-1)*scale(3)) + lim(3);

    igrid = igrid{4};
    iscenter = true;
else
    %otherwise, do it the normal way
    %grid can also specify non-square regions if it's six elements long,
    %with the first 3 elements being the size in x y and z and the second
    %3 elements being the step in x y and z
    if (length(igrid) == 2),
        igrid = igrid([1 1 1 2 2 2]);
    end;

    %get the grid size in real units
    vgrid = igrid.*[scale scale];
    %can't have a grid larger than the volume
    if (any(igrid(1:3) > sz(1:3))),
        warning('Grid sizes are larger than the volume.  Adjusting down.');
        k = find(igrid(1:3) > sz(1:3));
        igrid(k) = sz(k);
    end;
    iscenter = false;
end;

%set up the coordinates for each correlation volume
%centers of the volumes in x y and z
%the "i" prefix is for an index/pixel coordinate
irgnctr = floor(igrid/2);
nrgn = 2*irgnctr(1:3);

%pixel coordinates (irgn1) and real coordinates (rgnx)
irgn1 = -(irgnctr(1)-1):irgnctr(1);
rgnx = irgn1/scale(1);

irgn2 = -(irgnctr(2)-1):irgnctr(2);
rgny = irgn2/scale(2);

irgn3 = -(irgnctr(3)-1):irgnctr(3);
rgnz = irgn3/scale(3);

%generate the centers of the correlation volumes so that they're centered
%within the entire volume
if (~iscenter),
    igridstep = igrid(4:6);

    ictr1 = round(irgnctr(1)+1:igridstep(1):sz(2)-irgnctr(1));
    ictr1 = ictr1 + round((sz(2)-irgnctr(1)-ictr1(end))/2);
    ictr2 = round(irgnctr(2)+1:igridstep(2):sz(1)-irgnctr(2));
    ictr2 = ictr2 + round((sz(1)-irgnctr(2)-ictr2(end))/2);
    ictr3 = round(irgnctr(3)+1:igridstep(3):sz(3)-irgnctr(3));
    ictr3 = ictr3 + round((sz(3)-irgnctr(3)-ictr3(end))/2);

    %get the real coordinates
    ctrx = (ictr1-1)*scale(1) + lim(1,1);
    ctry = (ictr2-1)*scale(2) + lim(1,2);
    ctrz = (ictr3-1)*scale(3) + lim(1,3);
end;

% h = (-3:3)';
% ih = repmat(irgnctr,[length(h) 1]) + repmat(h,[1 3]);
difford = 5;
hctr = ceil(length(rgnx)/2);

%set up the windowing transform (kaiser)
w1 = kaiser(length(rgnx),0.5);
w2 = kaiser(length(rgny),0.5);
w3 = kaiser(length(rgnz),0.5);
wind = repmat(w1,[1 length(w2) length(w3)]) .* ...
    repmat(w2',[length(w1) 1 length(w3)]) .* ...
    repmat(shiftdim(w3,-2),[length(w1) length(w2) 1]);

diff1 = [];
diff2 = [];

if (opt.savecorr)
    h5opencreate(h5outfile,'/Correlations',[length(rgnx) length(rgny) length(rgnz) ...
        length(ictr1) length(ictr2) length(ictr3)], ...
        'ChunkSize',[length(rgnx) length(rgny) length(rgnz) 1 1 1]);
    h5opencreate(h5outfile,'/Hessians',[3 3 ...
        length(ictr1) length(ictr2) length(ictr3)], ...
        'ChunkSize',[3 3 1 1 1]);
end;
if (opt.saveeig)
    h5opencreate(h5outfile,'/Eigenvalues',[3 ...
        length(ictr1) length(ictr2) length(ictr3)], ...
        'ChunkSize',[3 1 1 1]);
end    
h5opencreate(h5outfile,'/Eigenvectors',[3 3 ...
    length(ictr1) length(ictr2) length(ictr3)], ...
    'ChunkSize',[3 1 1 1 1]);
h5opencreate(h5outfile,'/Grid/X', size(ctrx));
h5write(h5outfile,'/Grid/X', ctrx);
h5opencreate(h5outfile,'/Grid/Y', size(ctry));
h5write(h5outfile,'/Grid/Y', ctry);
h5opencreate(h5outfile,'/Grid/Z', size(ctrz));
h5write(h5outfile,'/Grid/Z', ctrz);
h5writeatt(h5outfile,'/Grid','CorrelationVolume',vgrid(1:3));

if strcmp(opt.method,'hough')
    [xvote,yvote,zvote] = ndgrid(irgn1,irgn2,irgn3);
    mag = sqrt(xvote.^2 + yvote.^2 + zvote.^2);
    xvote = xvote ./ mag;
    yvote = yvote ./ mag;
    zvote = zvote ./ mag;
end

a = 1;
N = length(ictr1)*length(ictr2)*length(ictr3);
progress(a,N, 'Computing fiber angles');
for k = 1:length(ictr3),
    for i = 1:length(ictr1),
        for j = 1:length(ictr2),
            V1 = h5read(h5imagefile,'/Image',[ictr1(i)-irgnctr(1) ictr2(j)-irgnctr(2) ictr3(k)-irgnctr(3)], ...
                nrgn);
            V1 = double(V1) / imagemax;
            
            %check to make sure that there's a signal in this region
            if (~debug && (rms(V1(:)) > 0.05))
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

                if strcmp(opt.method,'hough')
                    C = C / sum(C(:));
                    x1 = sum(xvote(:).*C(:));
                    y1 = sum(yvote(:).*C(:));
                    z1 = sum(zvote(:).*C(:));
                    
                    % CONTINUE - check significance of (x1,y1,z1) vector
                    % then project onto plane perpendicular to vector and
                    % get mean vector
                else
                    
                    %get the hessian at the center of the correlation matrix
                    H = hessian3(rgnx,rgny, rgnz, C, 'difference'); %,difford,diff1,diff2);
                    
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
                    [~,ind] = min(abs(ed));
                end
                
                %save correlation matrix and/or eigenvalues
                if (savecorr),
                    h5write(h5outfile,'/Correlations',C,[1 1 1 j i k],...
                        [size(C) 1 1 1]);
                    h5write(h5outfile,'/Hessians',H,[1 1 j i k],...
                        [size(H) 1 1 1]);
                end;
                if (saveeig),
                    h5write(h5outfile,'/Eigenvalues',ed([ind 1:ind-1 ind+1:end]),[1 j i k],...
                        [3 1 1 1]);
                end
                h5write(h5outfile,'/Eigenvectors',ev(:,[ind 1:ind-1 ind+1:end]),[1 1 j i k],...
                    [3 3 1 1 1]);
            end;
            a = a+1;
            progress(a,N);
        end;
    end;
end;

