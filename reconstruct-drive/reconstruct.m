%% manually reconstruct electrode array
% see www.mit.edu/~jvoigts for documentation
%
% this code is mostly untested and intended only as
% example for how to apply this method. 
%
% apr 2013 jakob voigts (jvoigts@mit.edu)

%% count images

readInDirectory=['/home/jvoigts/Dropbox/drive_reconstruction/nt_a9/scan2/'];

files = dir([readInDirectory '*.png']);
numImages=numel(files);

stack= [];
c=0;
reconstruct.coord=[];

%% LOAD
load(fullfile(readInDirectory,'tracked_points.mat'))

%% cache images
imgs=[];
for i=1:numImages
    i
    I=imread([readInDirectory 'scan_' num2str(i) '.png']);
    imgs(:,:,i)=I;
end;

%% init /settings
run=1;
c=1;
i=1;
lastxy=[ 0 0];

reconstruct.Npoints=40; 
reconstruct.coord=zeros(reconstruct.Npoints,numImages,3);

for i=1:reconstruct.Npoints
reconstruct.labels{i}='N';
end;

%% do manual reconstruction

while run
    %disp(i);
    
    reconstruct.coord(reconstruct.coord==0)=NaN;
    
    i=min(numImages,max(i,1))
    %  fnum=ceil(rand*(numImages-2));
    %fnum=i;
    
    I=imgs(:,:,i);
    clf;
    imagesc(I); hold on;
    title(['frame ' num2str(i) ]);
    
    % plot selection
    for j=1:reconstruct.Npoints
        
        py=(j/reconstruct.Npoints)*size(imgs,1);
        py_t=((j+1)/reconstruct.Npoints)*size(imgs,1);
        
        col=hsv2rgb((j/reconstruct.Npoints),1,.7);
        if c==j % if this one is selected
            
            fill([-20 -0 -0 -20],[py py py_t py_t]-((1/reconstruct.Npoints)*size(imgs,1)),'b','FaceColor',[1 1 1].*.7);
            
        end;
        text(-15,py-5,['P_{' num2str(j) '}'],'color',col);
        plot([-50,0],[1 1].*py,'k');
        plot([-50,0],[1 1].*py_t,'k');
        plot([-20 -20],[0 size(imgs,1)],'k');
        plot([-40 -40],[0 size(imgs,1)],'k');
        
        text(-30,py-5,[reconstruct.labels{j}]);
        
        text(-47,py-5,'x');
        
        scaleby=[min(reconstruct.coord) max(reconstruct.coord)];
        
        tp=(i/numImages)*size(imgs,2);
        plot([tp tp],[0 -100],'k')
        
        plot(linspace(0,size(imgs,2),numImages),-50*(([reconstruct.coord(j,:,1)]-min(reconstruct.coord(:)))./range(reconstruct.coord(:))),'-','color',col);
        plot(linspace(0,size(imgs,2),numImages),-50*(([reconstruct.coord(j,:,2)]-min(reconstruct.coord(:)))./range(reconstruct.coord(:)))-50,'-','color',col);
        
    end;
    
    
    set(gca, 'position', [0 0 1 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    
    %drawnow;
    daspect([1 1 1]);
    colormap(gray)
    
    for cc=1:reconstruct.Npoints
      %  if 1%cc~=c
             col=hsv2rgb((cc/reconstruct.Npoints),1,.7);
            plot(reconstruct.coord(cc,:,1),reconstruct.coord(cc,:,2),'r.','MarkerSize',4,'color',col);
            plot(reconstruct.coord(cc,i,1),reconstruct.coord(cc,i,2),'bo','MarkerSize',6,'color',col);
       % end
    end;
    
    
    plot(reconstruct.coord(c,:,1),reconstruct.coord(c,:,2),'r-');
    plot(reconstruct.coord(c,i,1),reconstruct.coord(c,i,2),'bo','MarkerSize',16,'color',col);
    
    
    
    
    % plot closeup of last position for precision
    
    wsize=100;
    close_x=lastxy(1)+[-wsize:wsize];
    close_y=lastxy(2)+[-wsize:wsize];
    close_x=round(max(min(close_x,size(I,2)),1));
    close_y=round(max(min(close_y,size(I,1)),1));
    
    Iclose=I(close_y,close_x);
    
    imagesc([0 wsize*4]+size(I,2),[0 wsize*4],Iclose);
    plot( [wsize*2]+size(I,2),[0 wsize*2],'r.','MarkerSize',2)
    
    xlim([-50 size(imgs,2)+wsize*4]);
    ylim([-100 size(imgs,1)]);
    
    [x,y,m]=ginput(1)
    
    switch m
        case 120
            run=0;
        case 3 % right mouse
          %  reconstruct.coord(c,i,:)=[x,y,3];
            i=i+1;
              lastxy=[x y];
        case 1 % left mouse
            if (y<0) % top timeline
                i=round((x/size(imgs,2))*numImages);
            end
            
            if (x>size(imgs,2))&&(y>0) % closeup image
                xs=lastxy(1)+( ((x-size(imgs,2))-(wsize*2))/2 );
                ys=lastxy(2)+( (y-(wsize*2))/2 );
                  reconstruct.coord(c,i,:)=[xs,ys,1];
                  i=i+1;
                  m=0.5;
                  lastxy=((m.*[xs ys])+((1-m).*lastxy));
                
            end;
            if (x>0)&&(y>0)&& (x<size(imgs,2)) % main image
                 reconstruct.coord(c,i,:)=[x,y,1];
                 i=i+1;
                  lastxy=[x y];
            else % menu on left
                
                
                for j=1:reconstruct.Npoints
                    py=((j-1)/reconstruct.Npoints)*size(imgs,1);
                    py_t=((j)/reconstruct.Npoints)*size(imgs,1);
                    
                    if  (x<0)&&(x>-20) % direct click on point selector
                        if (y<py_t)&&(y>=py)
                            c=j;
                            
                        end
                    end
                    
                    if    (x<-20)&&(x>-40) % click on label
                        if (y<py_t)&&(y>=py) % slect which point was clocked
                            
                            prompt={'Enter label'};
                            name='Input';
                            numlines=1;
                            defaultanswer={'x'};
                            
                            answer=inputdlg(prompt,name,numlines,defaultanswer);
                            
                            reconstruct.labels{j}=answer;
                        end;
                    end
                    
                    if (x<-40)&&(x>-50) % delete button
                        
                        if (y<py_t)&&(y>=py) % slect which point was clocked
                            button = questdlg('really delete?','delete?','delete','rather not','nope');
                            if strcmp(button,'delete')
                                
                                reconstruct.coord(j,:,:)=NaN;
                            end;
                        end;
                        
                        
                    end;
                end;
                
                
            end;
        case 3
            reconstruct.coord(c,i,:)=[x,y,3];
        case 102 %  Fadd and fwd
            reconstruct.coord(c,i,:)=[x,y,1];
            i=i+1;
        case 100 % D add and bkw
            reconstruct.coord(c,i,:)=[x,y,1];
            i=i-1;
        case 101 %  E bkw
            i=i-1;
        case 114 % R  bkw
            i=i+1;
            
        otherwise
    end;
    
  
    
    %ts(i)=mean(mean(I(370:376,272:287)));
end

%% save
save(fullfile(readInDirectory,'tracked_points.mat'),'reconstruct');

%% preprocess reconstruct.coordinates

p_c=[];
p_c(:,:,1)=reconstruct.coord(:,:,1);
p_c(:,:,2)=reconstruct.coord(:,:,2);

p_c(p_c==0)=NaN;

plot(p_c(:,:,1)',p_c(:,:,2)');

for i=1:size(reconstruct.coord,1)
    for c=1:2
        x=p_c(i,:,c);
        
        x(1:3)=NaN; % remove first few where we can't fully trust servo angles
        x(end-2:end)=NaN;
        
        in=find(isnan(x));
        nn=find(~isnan(x));
        if numel(nn)>3
            sp = spaps([1:numel(x)],x,80) ;
            xi = fnval(sp,[1:numel(x)]);
            xi(in)=NaN;
            p_c(i,:,c)=xi;
        end;
    end;

end;

ptss=[];
c=0;
for i=1:reconstruct.Npoints
    c=c+1;
    %ptss=reshape(p_c,[reconstruct.Npoints*2 numImages 1]);
    ptss(:,c)=p_c(i,:,1);
end;
for i=1:reconstruct.Npoints
    c=c+1;
    ptss(:,c)=p_c(i,:,2);
    
end;
clf;
plot(ptss);
plot([p_c(:,:,1) p_c(:,:,2)]')

ylabel('x/y position');
xlabel('rotation(deg)');

pts=ptss; % copy from data



%%  scale coordinates 
figure(1);
clf; hold on;


phi=deg2rad(15);
ox=0;
oy=0;


X=[];
X(1)=ox;X(2)=oy;X(3)=phi;

N=40;

X(1)=.4; % just make up some random ones fgor testing/vis
X(2)=.7;
X(3)=-16;%phi
X(4)=-20; %el

X([1:N]+4)=0; % initialize empty at 0
X([1:N]+4+N)=0;
X([1:N]+4+(N*2))=0;

X([1:4]+4)=  [.6  .8  0 0];
X([1:4]+4+N)=[ 0.1  .6  0 .6];
X([1:4]+4+(N*2))=[ 0.1 0.1 .1 .1];

%X=Xest;


pts_est=genmodel(X);


pts=pts./(  max(abs(pts(:)))  * 2); %scale coordinates


figure(3); clf; hold on;

plot(pts_est(:,1:(end/2)),pts_est(:,(end/2)+1:end));
plot(pts(:,1:(end/2)),pts(:,(end/2)+1:end),'k');

plot(pts_est(1,1:(end/2)),pts_est(1,(end/2)+1:end),'ko');
plot(pts(13,1:(end/2)),pts(13,(end/2)+1:end),'ro');

daspect([1 1 1]);

%% fit data

x0=ones(size(X));
x0(1)=230;
x0(2)=400;
x0(4)=30;
%x0(5)=0; % fov

x0=X;

costfun= @(x) sqrt(sum(sum(nan2zero(genmodel(x)- pts(1:8:end,:) ).^2))); % ignore missing data here

options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','iter');
[Xest,fval] = fminunc(costfun,x0,options);

figure(2); clf; hold on;
plot(Xest);
plot(X,'r');

figure(3); clf; hold on;
pts_est=genmodel(Xest);
plot(pts_est(:,1:(end/2)),pts_est(:,(end/2)+1:end));
plot(pts(:,1:(end/2)),pts(:,(end/2)+1:end),'k');


daspect([1 1 1]);
xe=Xest([1:N]+4);
ye=Xest([1:N]+4+N);
ze=Xest([1:N]+4+(N*2));

%% calibrate EDIT IDs FOR CALIBRATION POINTS MANUALLY HERE
% picking 2 that are 2mm apart (or multiple sets of 2)
c=0;
for i=1:reconstruct.Npoints
    if strcmp(reconstruct.labels{i},'c')
        c=c+1;
       cp(c)=i;
    end
end;

c=0;


c=c+1; % first calibration measurement
from=17; to=20; % EDIT HERE
dist_mm=2;  %mm
px2mm(c)=dist_mm/norm([xe(from),ye(from),ze(from)]-[xe(to),ye(to),ze(to)],2);


c=c+1; % first calibration measurement
from=18; to=19; % EDIT HERE
dist_mm=2;  %mm
px2mm(c)=dist_mm/norm([xe(from),ye(from),ze(from)]-[xe(to),ye(to),ze(to)],2);


c=c+1; % first calibration measurement
from=17 ; to=18; % EDIT HERE
dist_mm=1;  %mm
px2mm(c)=dist_mm/norm([xe(from),ye(from),ze(from)]-[xe(to),ye(to),ze(to)],2);

c=c+1; % first calibration measurement
from=19 ; to=20; % EDIT HERE
dist_mm=1;  %mm
px2mm(c)=dist_mm/norm([xe(from),ye(from),ze(from)]-[xe(to),ye(to),ze(to)],2);

disp(['    got ' num2str(numel(px2mm)) ' measurements:']);
disp(px2mm);

px2mm=mean(px2mm);

xe=xe.*px2mm;
ye=ye.*px2mm;
ze=ze.*px2mm;

%% plot 3d results
figure(200); clf; hold on;

%plot3(x,y,z,'ro');

for i=1:reconstruct.Npoints
    if ~strcmp(reconstruct.labels{i},'N')
        plot3(xe(i),ye(i),ze(i),'go');
        text(xe(i),ye(i),ze(i),reconstruct.labels{i})
        for j=1:reconstruct.Npoints
            if strcmp(reconstruct.labels{j}, strcat(reconstruct.labels{i},'b'))
             plot3(xe([i j]),ye([i j]),ze([i j]),'b');
            end
        end
        
    end
     if strcmp(reconstruct.labels{i},'c')
         text(xe(i),ye(i),ze(i)+.05,num2str(i))
         
     end
    
end;

edges=[];
c=0;
for i=1:reconstruct.Npoints
    if strcmp(reconstruct.labels{i},'e')
        c=c+1;
        edges(c,:)=[xe(i),ye(i),ze(i)];
    end
end;
ab=[1 2;1 3; 1 4; 2 1; 2 3; 2 4; 3 1; 3 2; 3 4];
for i=1:size(ab,1);
    a=ab(i,1);
    b=ab(i,2);
    plot3(edges([a b],1),edges([a b],2),edges([a b],3),'-','color',[1 1 1].*.3);
end;

edges=[];
c=0;
for i=1:reconstruct.Npoints
    if strcmp(reconstruct.labels{i},'c')
        c=c+1;
        edges(c,:)=[xe(i),ye(i),ze(i)];
    end
end;
ab=[1 2;1 3; 1 4; 2 1; 2 3; 2 4; 3 1; 3 2; 3 4];
for i=1:size(ab,1);
    a=ab(i,1);
    b=ab(i,2);
    plot3(edges([a b],1),edges([a b],2),edges([a b],3),'-','color',[1 1 1].*.7);
end;


daspect([1 1 1])
grid on

reconstructed_3d=[];
reconstructed_3d.xe=xe;
reconstructed_3d.ye=ye;
reconstructed_3d.ze=ze;


%% save
save(fullfile(readInDirectory,'reconstructed_points.mat'),'reconstructed_3d');


