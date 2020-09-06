

% -------------------------------------------------------------------------
% Matlab code developed at National University of Singapore TMSI/PORL
% on 2013/05/06
% -------------------------------------------------------------------------
%     Supporter: Pavel Tkalich
%     Coder: Luu Quang Hung
%     Email: luuquanghung@gmail.com
% -------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE NAO99 TIDE DATA FOR TPX GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pathnao = 'F:\Singapo\Research\Works\Data\Tide\NAO\NAO99b\data\';
filetpx = 'F:\Singapo\Research\Works\Models\ROMS\RomsAgrif\Data\TPXO7\TPXO7.nc';
fileout = 'F:\Singapo\Research\Works\Models\ROMS\RomsAgrif\Data\TPXO7\NAO99.nc';
constituent = {'Sa'};
period = 365.2425; % days
% constituent = {'sa';'ssa';'mm';'msf';'mf'};
% period = [365.2425,182.6211,27.5546,14.7653,13.6608,];
constr = [];
nm = numel(constituent);
for n = 1:nm
    constr = [constr constituent{n} ' '];
end


% read input TPXO
ncid = netcdf.open(filetpx,'NC_NOWRITE');
[lon_r, lon_r_value] = netcdf.inqDim(ncid,0);
[lat_r, lat_r_value] = netcdf.inqDim(ncid,1);
[lon_u, lon_u_value] = netcdf.inqDim(ncid,2);
[lat_u, lat_u_value] = netcdf.inqDim(ncid,3);
[lon_v, lon_v_value] = netcdf.inqDim(ncid,4);
[lat_v, lat_v_value] = netcdf.inqDim(ncid,5);
netcdf.close(ncid);
nctides = netcdf(filetpx);
periods = nctides{'periods'}(:);
xr = nctides{'lon_r'}(:);
yr = nctides{'lat_r'}(:);
xu = nctides{'lon_u'}(:);
yu = nctides{'lat_u'}(:);
xv = nctides{'lon_v'}(:);
yv = nctides{'lat_v'}(:);
ur = nctides{'u_r'}(1,:,:);
ui = nctides{'u_i'}(1,:,:);
vr = nctides{'v_r'}(1,:,:);
vi = nctides{'v_i'}(1,:,:);
ssh = nctides{'ssh_r'}(1,:,:);
hz = nctides{'h'}(:);
close(nctides);


% read input NAO
filenao = [pathnao constituent{1} '.nao'];
[a,p,im,jm,dx,dy] = readnao(filenao);
lon = (0:im-1)*dx;
lat = (0:jm-1)*dy-90.;%[lon2,lat2]=meshgrid(lon,lat);
lon = squeeze(lon);
lat = squeeze(lat);
u_r = zeros(jm,im,nm);
u_i = zeros(jm,im,nm);
v_r = zeros(jm,im,nm);
v_i = zeros(jm,im,nm);
ssh_r = zeros(jm,im,nm);
ssh_i = zeros(jm,im,nm);
h = zeros(jm,im);
for n = 1:nm
    filenao = [pathnao constituent{n} '.nao'];
    [a,p] = readnao(filenao);        
    for i = 1:im
        for j = 1:jm
            ssh_r(j,i,n) = a(j,i)*cos(p(j,i));
            ssh_i(j,i,n) = a(j,i)*sin(p(j,i)*pi/180.);
            if n==1
                if a(j,i)>0 %topo is created using ssh grid, only need once time
                    h(j,i) = 1;
                else
                    h(j,i) = 0;
                end
            end
        end
    end
    for i = 1:im
        for j = 1:jm
            if(h(j,i)~=1)
                u_r(j,i,n) = NaN;
                u_i(j,i,n) = NaN;
                v_r(j,i,n) = NaN;
                v_i(j,i,n) = NaN;
            end
        end
    end
end


% write output structure
nw = netcdf(fileout, 'clobber');
nw.title = 'NAO99b tide model similar to TPXO7';
nw.date = '06-May-2013';
nw.components = constr;
nw('lon_r') = im;
nw('lat_r') = jm;
nw('lon_u') = im;
nw('lat_u') = jm;
nw('lon_v') = im;
nw('lat_v') = jm;
nw('periods') = nm;
nw{'lon_r'} = ncdouble('lon_r');
nw{'lon_r'}.long_name = ncchar('Longitude at SSH points');
nw{'lon_r'}.long_name = 'Longitude at SSH points';
nw{'lat_r'} = ncdouble('lat_r');
nw{'lat_r'}.long_name = ncchar('Latitude at SSH points');
nw{'lat_r'}.long_name = 'Latitude at SSH points';
nw{'lon_u'} = ncdouble('lon_u');
nw{'lon_u'}.long_name = ncchar('Longitude at U points');
nw{'lon_u'}.long_name = 'Longitude at U points';
nw{'lat_u'} = ncdouble('lat_u');
nw{'lat_u'}.long_name = ncchar('Latitude at U points');
nw{'lat_u'}.long_name = 'Latitude at U points';
nw{'lon_v'} = ncdouble('lon_v');
nw{'lon_v'}.long_name = ncchar('Longitude at V points');
nw{'lon_v'}.long_name = 'Longitude at V points';
nw{'lat_v'} = ncdouble('lat_v');
nw{'lat_v'}.long_name = ncchar('Latitude at V points');
nw{'lat_v'}.long_name = 'Latitude at V points';
nw{'periods'} = ncdouble('periods');
nw{'periods'}.long_name = ncchar('Tide periods');
nw{'periods'}.long_name = 'Tide periods';
nw{'h'} = ncdouble('lat_r','lon_r');
nw{'h'}.long_name = ncchar('Topography');
nw{'h'}.long_name = 'Topography';
nw{'ssh_r'} = ncdouble('periods','lat_r','lon_r');
nw{'ssh_r'}.long_name = ncchar('Elevation real part');
nw{'ssh_r'}.long_name = 'Elevation real part';
nw{'ssh_i'} = ncdouble('periods','lat_r','lon_r');
nw{'ssh_i'}.long_name = ncchar('Elevation imaginary part');
nw{'ssh_i'}.long_name = 'Elevation imaginary part';
nw{'u_r'} = ncdouble('periods','lat_u','lon_u');
nw{'u_r'}.long_name = ncchar('U-transport component real part');
nw{'u_r'}.long_name = 'U-transport component real part';
nw{'u_i'} = ncdouble('periods','lat_u','lon_u');
nw{'u_i'}.long_name = ncchar('U-transport component imaginary part');
nw{'u_i'}.long_name = 'U-transport component imaginary part';
nw{'v_r'} = ncdouble('periods','lat_v','lon_v');
nw{'v_r'}.long_name = ncchar('V-transport component real part');
nw{'v_r'}.long_name = 'V-transport component real part';
nw{'v_i'} = ncdouble('periods','lat_v','lon_v');
nw{'v_i'}.long_name = ncchar('V-transport component imaginary part');
nw{'v_i'}.long_name = 'V-transport component imaginary part';


% write output data
nw{'lon_r'}(:) = lon;
nw{'lat_r'}(:) = lat;
nw{'lon_u'}(:) = lon;
nw{'lat_u'}(:) = lat;
nw{'lon_v'}(:) = lon;
nw{'lat_v'}(:) = lat;
nw{'u_r'}(:,:,:) = u_r;
nw{'u_i'}(:,:,:) = u_i;
nw{'v_r'}(:,:,:) = v_r;
nw{'v_i'}(:,:,:) = v_i;
nw{'ssh_r'}(:,:,:) = ssh_r;
nw{'ssh_i'}(:,:,:) = ssh_i;
nw{'h'}(:,:) = h;
nw{'periods'}(:) = period*24.;
close(nw);



















