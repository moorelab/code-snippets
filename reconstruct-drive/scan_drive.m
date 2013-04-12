%% scan electrode array
% see www.mit.edu/~jvoigts for documentation
%
% this code is mostly untested and intended only as
% example for how to apply this method. 
%
% apr 2013 jakob voigts (jvoigts@mit.edu)

%% connecto to arduino, using the MATLAB arduino target code
delete(instrfind({'Port'},{'COM4'}))
a=arduino('COM4')
a.servoAttach(9);


%%
imaqreset;
%% set up avt pike
vid = videoinput('avtmatlabadaptor64_r2009b', 1, 'F0M5_Mono8_640x480');
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;
vid.LoggingMode = 'memory';

src.Shutter = 800; % edit this for different exposures

%% start preview
preview(vid);

%% rotate in full arc to check centering and focus

deg_from=10; % define arc here
deg_to=180;

for i=1:5
    
    a.servoWrite(9,deg_from);
    pause(2);
    a.servoWrite(9,deg_to);
    pause(2);
end;

%% stop preview
stoppreview(vid);

%% run scan
c=0;

outpath='C:\Users\jvoigts\Dropbox\drive_reconstruction\drive2\scan1\'
a.servoWrite(9,deg_to);
pause(1);
mkdir(outpath);

for deg=deg_from:deg_to
    c=c+1;
    a.servoWrite(9,max(0,deg-20)); % run back and forth to increase precision (servos get stuck on small adjustments sometimes)
    pause(.5);
    a.servoWrite(9,deg);
    pause(3);
    
    start(vid);
    pause(1);
    im = squeeze(getdata(vid,vid.FramesAvailable));
    
    imagesc(im);
    colormap(gray)
    daspect([1 1 1]);
    title(['deg',num2str(deg)]);
    drawnow;
    
    imwrite(im,fullfile(outpath,['scan_',num2str(c),'.png']));
    
end;

