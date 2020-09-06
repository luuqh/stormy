function download_GFS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  download_GFS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig)
%
%
%  Extract a subgrid from GFS to get a ROMS forcing
%  Store that into monthly files (to limit the problems
%  of bandwith...).
%  Take care of the Greenwitch Meridian.
%  Transform variables in the ROMS format.
% 
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    9-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
romstools_param

%
% Put the date in 'Yorig' time
%
rundate=datenum(today)-datenum(Yorig,1,1);
%
% GFS output name
%
gfs_name=[FRCST_dir,'GFS_',num2str(rundate),'.nc'];
%
% start
%
disp([' '])
disp(['Get GFS data for ',datestr(today)])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])
%
% Create the directory
%
disp(['Making output data directory ',FRCST_dir])
eval(['!mkdir ',FRCST_dir])
%
% Get the GFS file name (first check if there is an available forecast)
%
gfs_run_time0=0;
gfs_date0=today-0.5;
gfs_run_time=gfs_run_time0;
gfs_date=gfs_date0;
foundfile=0;
while foundfile==0
  fname=get_GFS_fname(gfs_date,gfs_run_time,1);
  warning off
  try
    x=loaddap('-A -e +v ', fname);
    foundfile=1;
  catch
    foundfile=0;
  end
  if foundfile==1 & ~isempty(x)
    disp('  File found')
  else
    foundfile=0;
    disp(['  GFS : did not found ',fname])
    gfs_run_time=gfs_run_time-6;
    if gfs_run_time<0
      gfs_date=gfs_date-1;
      gfs_run_time=18;
      if gfs_date<gfs_date0-8;
        error(['  GFS: did not found anything'])
      end
    end 
  end
  warning on
end
%
% Get the grid and the indices of the subgrid for GFS
%
gfs_run_time_GFS=gfs_run_time; 
gfs_date_GFS=gfs_date;  
fname=get_GFS_fname(gfs_date_GFS,gfs_run_time_GFS,1); 
[i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat]=...
         get_GFS_subgrid(fname,lonmin,lonmax,latmin,latmax);
mask=getdap('',fname,'landsfc','[1:1]','',jrange,...
             i1min,i1max,i2min,i2max,i3min,i3max);
mask(mask==1)=NaN;
mask(isfinite(mask))=1;
%
% Initialisation
N=70;
[M,L]=size(mask);
tx=zeros(N,M,L);
ty=tx;
uwnd=tx;
vwnd=tx;
%
n=0;
gfstime=0*(1:N);
LON=lon;
LAT=lat;
%
%==================================================
% Get the variables for hindcast
%==================================================
%
% Get GDAS 1 deg grid
%
gfs_run_time=12;
gfs_date=today-(hdays+1);
fname=get_GFS_fname(gfs_date,gfs_run_time,0);
[i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat]=...
         get_GFS_subgrid(fname,lonmin,lonmax,latmin,latmax);
%mask=getdap('',fname,'landsfc','[1:1]','',jrange,...
mask=getdap('',fname,'landsfc','[0:0]','',jrange,...
             i1min,i1max,i2min,i2max,i3min,i3max);
mask(mask==1)=NaN;
mask(isfinite(mask))=1;
[M,L]=size(mask);
%
% Loop on GDAS analyses 
% (starts hdays day ago 18Z, 1 analysis every 6h
%  ends before yesterday 18Z).
%
gfs_run_time0=12;
gfs_date0=today-(hdays+1); 
for frcst=1:4*hdays-3             % number of files until before yesterday 18Z 
  gfs_run_time0=gfs_run_time0+6;  
  if gfs_run_time0>18
    gfs_date0=gfs_date0+1;
    gfs_run_time0=0;
  end
%
% 1.1: check if the GDAS is available.
%
  gfs_run_time=gfs_run_time0;
  gfs_date=gfs_date0;
  t1=1;
  t1dap=t1-1;
  fname=get_GFS_fname(gfs_date,gfs_run_time,0);
  warning off
  try
    x=loaddap('-A -e +v ', fname);
    foundfile=1;
  catch
    foundfile=0;
  end
  if foundfile==1 & ~isempty(x)
    disp('  File found')
    gfs_run_time1=gfs_run_time0;  % Increment Forecast 
    gfs_date1=gfs_date0;          %           start time
  else
    foundfile=0;
    disp(['  GDAS: file not found, try next file'])
%    error(['  GDAS: file not found'])
  end
  warning on
%
% 1.2: read file 
%
  if foundfile, 
   missvalue=x.uflxsfc.missing_value;
% 
   n=n+1;
   disp(['N=',num2str(n)])
   [a(n),tx0,ty0,pressfc]=...
   get_GDAS(fname,mask,t1dap,jrange,i1min,i1max,i2min,i2max,i3min,i3max,missvalue);
   TTX=interp2(lon,lat',tx0,LON,LAT');
   TTY=interp2(lon,lat',ty0,LON,LAT');
   SSLP=interp2(lon,lat',pressfc,LON,LAT');

   gfstime=a(n);
   TX(n,:,:)=TTX;
   TY(n,:,:)=TTY;
   SLP(n,:,:)=SSLP;
%
  end 
end
%
%==================================================================
% 2: Get the variables for Forecast starting  yesterday 00Z
%==================================================================
gfs_run_time0=gfs_run_time1+6;
gfs_date0=gfs_date1;
  if gfs_run_time0>18
    gfs_date0=gfs_date0+1;
    gfs_run_time0=0;
  end
gfs_run_time=gfs_run_time0;
gfs_date=gfs_date0;

% Get the grid and the indices of the subgrid for GFS
fname=get_GFS_fname(gfs_date,gfs_run_time,1);
[i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat]=...
         get_GFS_subgrid(fname,lonmin,lonmax,latmin,latmax);

mask=getdap('',fname,'landsfc','[0:0]','',jrange,...
             i1min,i1max,i2min,i2max,i3min,i3max);
mask(mask==1)=NaN;
mask(isfinite(mask))=1;
[M,L]=size(mask);
%
% 2.1: check if the GFS forecast has been done for this time.
% if not: take the previous one (but increment time index)
%
t1=1;
foundfile=0;
while foundfile==0
  fname=get_GFS_fname(gfs_date,gfs_run_time,1);
  warning off
  try
    x=loaddap('-A -e +v ', fname);
    foundfile=1;
  catch
    foundfile=0;
  end
  if foundfile==1 & ~isempty(x)
    disp('  File found')
  else
    foundfile=0;
    disp(['  GFS : did not found ',fname])
    %%t1=t1+it; % increment time index  POSSIBLE BUG
	          % if it=2 it is Ok, if it=1, (*--> 3h00)
			  % need to to add 2*it ...
	t1=t1+2;
    gfs_run_time=gfs_run_time-6;
    if gfs_run_time<0
      gfs_date=gfs_date-1;
      gfs_run_time=18;
      if gfs_date<gfs_date0-8;
        error(['  GFS: did not found anything'])
      end
    end
  end
  warning on
end
%
% 2.2: read 60 time steps (dt = 3h)
%

for tndx=t1:it:60
  tndxdap=tndx-1;
  n=n+1;
  disp(['tndxdap=',num2str(tndxdap)])
  disp(['tndx=',num2str(tndx)])
  disp(['n=',num2str(n)])
  disp(['N_gfs=',num2str(n)])
%  [gfstime(n),tx(n,:,:),ty(n,:,:),...
%  pressfc(n,:,:)]=...
%  [gfstime(n),tx(n,:,:), ty(n,:,:), pressfc(n,:,:)]=...
  [a(n),b,c,d]=...
   get_GFS(fname,mask,tndxdap,jrange,i1min,i1max,i2min,i2max,...
           i3min,i3max,missvalue);

   TX(n,:,:)=interp2(lon,lat',b,LON,LAT');
   TY(n,:,:)=interp2(lon,lat',c,LON,LAT');
   SLP(n,:,:)=interp2(lon,lat',d,LON,LAT');
%   keyboard
end
%
% Reduce the matrices
%
gfstime=a(1:n);
tx=TX(1:n,:,:);
ty=TY(1:n,:,:);
slp=SLP(1:n,:,:);
%
% Put the time in Yorig time
%
gfstime=gfstime-datenum(Yorig,1,1);
%
% Create the GFS output file and write everything down
%
mask(isnan(mask))=0;
write_GFS(gfs_name,Yorig,lon,lat,mask,gfstime,tx,ty,slp)
%
disp('Download GFS: done')
%
return
