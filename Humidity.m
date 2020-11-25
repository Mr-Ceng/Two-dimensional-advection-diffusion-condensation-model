clear;

%------define domain and time stepping variables------%
ny = 100;
nz = 100;
Ly = pi;
Lz = pi;
dy = Ly/ny;
dz = Lz/nz;
yc = linspace(dy/2, Ly-dy/2, ny);
zc = linspace(dy/2, Ly-dy/2, ny);
nt = 50000;
dt = 0.001;
tau = 0.1;


%----------define V, T, kappa, qs; examing CFL----------%
vy = zeros(nz, ny);
vz = zeros(nz, ny);
for k=1:nz
    for j=1:ny
        vy(k,j) = -sin(yc(j))*cos(zc(k));
        vz(k,j) = cos(yc(j))*sin(zc(k));
    end
end
Pe = 1000;
kappa = 1*Ly/Pe
CFL = dy/(1*dt)
delT = 76;
dT = delT/nz;
T = linspace(223.15+dT/2, 223.15+delT-dT/2, nz);
T = flip(T);
qs = zeros(nz,ny);
for j=1:ny
    qs(:,j) = 6.1078*exp(17.2693882.*(T-273.16)./(T-35.86))*0.622/1000;
end


%----------------------------integrate q----------------------------%
q = zeros(nt,nz,ny);
q(1,1,:) = qs(1,:);


%first stepping, Heun2
qtrend = calc_qtrend(squeeze(q(1,:,:)), qs, vy, vz, kappa, tau, nz, ny, dz, dy);
qest = squeeze(q(1,:,:)) + dt * qtrend;
qesttrend = calc_qtrend(qest, qs, vy, vz, kappa, tau, nz, ny, dz, dy);
q(2,:,:) = squeeze(q(1,:,:)) + 0.5 * dt * (qtrend + qesttrend);
q(2,nz,:) = q(2,nz-1,:);
q(2,:,1) = q(2,:,2);
q(2,:,ny) = q(2,:,ny-1);

%time stepping, Adams-Bashforth 2
for i = 3:nt
    qtrend1 = calc_qtrend(squeeze(q(i-1,:,:)), qs, vy, vz, kappa, tau, nz, ny, dz, dy);
    qtrend2 = calc_qtrend(squeeze(q(i-2,:,:)), qs, vy, vz, kappa, tau, nz, ny, dz, dy);
    q(i,:,:) = AB2(squeeze(q(i-1,:,:)), qtrend1, qtrend2, dt);
    q(i,nz,:) = q(i,nz-1,:);
    q(i,:,1) = q(i,:,2);
    q(i,:,ny) = q(i,:,ny-1);
end

%----------------------------Plotting----------------------------%

q_end = squeeze(q(nt/2,:,:));
q_end(q_end<0.) = 0.000000000001;
logq = log(q_end);
q_trend = calc_qtrend(squeeze(q(nt,:,:)), qs, vy, vz, kappa, tau, nz, ny, dz, dy) * dt;
RH = q_end./qs;
Precip = (q_end - qs) /tau;
Precip(Precip<0.) = 0.;

yc = yc/pi;
zc = zc/pi;


figure(1);
colormap(flipud(bone));
subplot(2,2,1);
contourf(yc(2:ny-1),zc(2:nz-1),RH(2:nz-1,2:ny-1),'LineColor','none');
colo = colorbar;
%caxis([0 1*1e-8]);
title('Relative humidity');
xlabel('y');
ylabel('z');

subplot(2,2,2);
contourf(yc(2:ny-1),zc(2:nz-1),logq(2:nz-1,2:ny-1),'LineColor','none');
colo = colorbar;
%caxis([0 1*1e-8]);
title('Water vapor (log)');
xlabel('y');
ylabel('z');

subplot(2,2,3);
contourf(yc(2:ny-1),zc(2:nz-1),Precip(2:nz-1,2:ny-1),'LineColor','none');
colo = colorbar;
%caxis([0 1*1e-8]);
title('Precipitation');
xlabel('y');
ylabel('z');


subplot(2,2,4);
contourf(yc(2:ny-1),zc(2:nz-1),q_trend(2:nz-1,2:ny-1),'LineColor','none');
colo = colorbar;
%caxis([0 1*1e-8]);
title('Tendency');
xlabel('y');
ylabel('z');



%---------------------------Functions---------------------------%
%calculate trend
function trend = calc_qtrend(x, qs, vy, vz, kappa, tau, nz, ny, dz, dy)
    xy = parx(x, 2, nz, ny, dz, dy);
    xy2 = parx2(x, 2, nz, ny, dz, dy);
    xz = parx(x, 1, nz, ny, dz, dy);
    xz2 = parx2(x, 1, nz, ny, dz, dy);
    S = (x - qs) / tau;
    S(S<0.) = 0.;
    trend = -vy.*xy -vz.*xz + kappa.*(xy2+xz2) - S;
end

% partial x, order 2
function varnew = parx(var, dim, nz, ny, dz, dy)
    varnew = zeros(nz, ny);
    if (dim == 1)
        varnew(2:nz-1,:) = ( var(3:nz,:) - var(1:nz-2,:) ) / ( 2*dz );
    elseif (dim == 2)
        varnew(:,2:ny-1) = ( var(:,3:ny) - var(:,1:ny-2) ) / ( 2*dy );
    else
        'Error: undefined dimension'
    end
end

% partial^2 x, order 2
function varnew = parx2(var, dim, nz, ny, dz, dy)
    varnew = zeros(nz, ny);
    if (dim == 1)
        varnew(2:nz-1,:)=(var(3:nz,:) - 2*var(2:nz-1,:) + var(1:nz-2,:))/(dz^2);
    elseif (dim == 2)
        varnew(:,2:ny-1)=(var(:,3:ny) - 2*var(:,2:ny-1) + var(:,1:ny-2))/(dy^2);
    else
        'Error: undefined dimension'
    end
end

%AB2
function xnew = AB2(x, xtrend1, xtrend2, dt)
    xnew = x + 0.5*dt*( 3*xtrend1 - xtrend2);
end










