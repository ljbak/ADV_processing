

%% velocities
    % frames=[];
     noisy=0;
     [u,v,x,y,t,tr]=Velocities(vtracks,[],noisy);
%% grid velocity data not phase averaged


close all

for i=1:length(vtracks)
    lens(i)=vtracks(i).len;
end

nx=100; %grid points in x
ny=100; %grid points in y
pix_y=1024; %pixels (same x & y in this case)
pix_x=1024;
%each pixel corresponds to the start (left or top) of it.
[X Y]=meshgrid(linspace(1,pix_x+1,nx),linspace(1,pix_y+1,ny));

    ind=y<150;
    u(ind)=nan;
    v(ind)=nan;
    
is=min(t):max(t);

ui=nan(ny,nx,length(is));
vi=ui;

n=2000;
for i=1:n %1:max(index)
    i/n*100
    %     i/max(t)*100
    %     ind=i==t;%tindex;
    lencutoff=5;
    % for i=1:max(index)
    %     ind=i==tindex;
    %     ind=ind & (lens(tr')'>=lencutoff)';
    
    ind=t==is(i);
    ind=ind & lens(tr')'>lencutoff;
    
    ui(:,:,i)=griddata(x(ind),y(ind),u(ind),X,Y);
    vi(:,:,i)=griddata(x(ind),y(ind),v(ind),X,Y);
    if yesplot
        %%PLOT
        contourf(X,Y,flipud(ui(:,:,i)),[-max_disp:.25:max_disp],'edgecolor','none')
        c=colorbar;
        colormap hot
        caxis([-max_disp max_disp])
        drawnow
        pause(.01)
        
    end
end
Umeani=nanmean(ui,3); %mean gridded data
Vmeani=nanmean(vi,3); %mean gridded data
%% fit fourier series & get reflection coefficients
i1=40;
i2=40;
uu=ui(i1,i2,:);
uu=uu(:)*fps/cm2pix/100;
uu=uu(~isnan(uu));
vv=vi(i1,i2,:);
vv=vv(:)*fps/cm2pix/100;
vv=vv(~isnan(vv));
tt=(1:length(vv))/fps*wavef*2*pi;

% 
% 
% 
% omega=2*pi*wavef;
% dat=ui(i1,i2,1:300);
% dat=dat(:)'*fps/cm2pix/100;
% Tphase=0:1/fps:(length(dat)-1)/fps;
% yu = max(dat);
% yl = min(dat);
% yr = (yu-yl)/2;                               % Range of ?y?
% yz = yl+(yr/2);
% ym = mean(dat);                               % Estimate offset
% fitt = @(b,Tphase)  b(1).*(cos(-omega*Tphase + b(2))) + b(3);    % Function to fit
% fcn = @(b) sum((fitt(b,Tphase) - dat).^2);                              % Least-Squares cost function
% s = fminsearch(fcn, [yr;  0;  ym]);                       % Minimise Least-Squares
% xp = linspace(min(Tphase),max(Tphase));
% [fitdata]=fitt(s,Tphase);


tt=-tt;%+s(2);

ind=1:1000;

f = fit(tt(ind)',uu(ind),'fourier1')
plot(f,tt,uu)

C=f.a1;
D=f.b1;

f = fit(tt(ind)',vv(ind),'fourier1')
plot(f,tt,vv)

G=f.a1;
H=f.b1;

E=G; %tau=0
F=H; %tau=0;


g=9.81;
H=42/100;

ks=[4.29 8.83 14.1];
k=ks(wc);

yy=(-[0:(ny-1)]*1024/ny- 1024/ny/2);
top=290; % ???
yy=(yy/cm2pix+top/cm2pix);
z=yy(i1)/100;  

Z=g*k/omega*cosh(k*H+k*z)/cosh(k*H);
Y=g*k/omega*sinh(k*H+k*z)/cosh(k*H);

a_i=(D/Z+E/Y)^2+(C/Z-F/Y)^2;
a_i=1/2*sqrt(a_i)/omega;
a_r=(D/Z-E/Y)^2+(C/Z+F/Y)^2;
a_r=1/2*sqrt(a_r)/omega;




