filename = 'Drerio4.h5';
datasetname = '/image';
outfile = ''; %'Drerio4-saggital.mp4';
fps = 10;

imageoutname = 'Drerio4-saggital';

info = h5info(filename,datasetname);
sz = info.Dataspace.Size;
fprintf('Total size: %d z, %d y, %d x\n', sz);

x = round((sz(3)+1)/2);

Iyz = h5read(filename,datasetname,[1 1 x],[sz(1) sz(2) 1]);
Iyz = flipud(squeeze(Iyz)');

fig = figureseries('saggital');
set(fig,'Units','pixels','WindowStyle','normal');
pos = get(fig,'Position');

w = 1024;
h = round(sz(2)/sz(1) * 1024);
pos(3:4) = [w h];

set(fig,'Position',pos);

clf;
axes('Position',[0 0 1 1]);
imshow(Iyz,'InitialMagnification','fit');

if ~isempty(imageoutname)
    imwrite(im2uint8(Iyz), [imageoutname '1.jpg']);
end

if ~isempty(outfile)
    vid = VideoWriter(outfile, 'MPEG-4');
    set(vid, 'FrameRate',fps);
    open(vid);
end

hold on;

dx = 1;
for i = 1:sz(3)*2
    Iyz = h5read(filename,datasetname,[1 1 x],[sz(1) sz(2) 1]);
    Iyz = flipud(squeeze(Iyz)');

    cla;
    imshow(Iyz,'InitialMagnification','fit');
    drawnow;

    if ~isempty(imageoutname) && (i == round(0.12*sz(3)))
        imwrite(im2uint8(Iyz), [imageoutname '2.jpg']);
    end
    
    if ~isempty(outfile)
        frame = getframe;
        writeVideo(vid,frame);
    end
    
    if (x == sz(3))
        dx = -1;
    elseif (x == 1)
        dx = 1;
    end
    x = x+dx;
end
if (~isempty(outfile))
    close(vid);
end

hold off;

