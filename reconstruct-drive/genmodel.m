function f = genmodel(X)
% see www.mit.edu/~jvoigts for documentation
%
% this code is mostly untested and intended only as
% example for how to apply this method. 
%
% apr 2013 jakob voigts (jvoigts@mit.edu)
%
% takes as parameters
% offsets: ox oy
% camera axial angle: phi
% camera elevation angle
% N* x,y,z coordinates
%
% known: camera azimuth angles

% out: N_angles X 2*N matrix, x first then y in image plane

ox=X(1);
oy=X(2);

phi=deg2rad(X(3));

el=X(4);

% estimate fov manually once for lens/spacer
%fov=1;  % 0.1674
fov=3; %  0.167
%fov=6; %  0.169481
%fov=8; % 0.17126 
%fov=10.2; %0.172
%fov=15; % 0.176
%fov=10;

N=40;
x=X([1:N]+4);
y=X([1:N]+4+N);
z=X([1:N]+4+(N*2));




rot=[cos(phi) -sin(phi) ;sin(phi) cos(phi)];
c=0;
for d=1:8:171
    c=c+1;
    A = viewmtx(d,el,fov);
    [m,n] = size(x);
    x4d = [x(:),y(:),z(:),ones(m*n,1)]';
    x2d = A*x4d;
    x2 = zeros(m,n); y2 = zeros(m,n);
    x2(:) = x2d(1,:)./x2d(4,:);
    y2(:) = x2d(2,:)./x2d(4,:);
    
    %rotate around camera view axis
    tmp=[x2;y2]'*rot;
    x2=tmp(:,1);y2=tmp(:,2);
    
    
    pts(c,:)=[x2'+ox,y2'+oy];
    
end;

f=pts;
