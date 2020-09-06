function  create_forcing_clim_tide (frcname,grdname,title,smst,slpt,smsc,slpc,ymd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf forcing file
%       frcname: name of the forcing file
%       grdname: name of the grid file
%       title: title in the netcdf file  
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
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(grdname);
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
result=close(nc);
Lp=L+1;
Mp=M+1;
disp([L M Lp Mp])

nw = netcdf(frcname, 'clobber');
% ncid = netcdf.create (frcname, 'nc_write');
redef(nw);
%disp(result);

%
%  Create dimensions
%
%Ntides = 1;
%components = 'Sa';

nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('sms_time') = length(smst);
nw('slp_time') = length(slpt);
%netcdf.defDim(ncid,'xi_n',L);
%netcdf.defDim(ncid,'eta_u',Mp);
%netcdf.defDim(ncid,'xi_v',Lp);
%netcdf.defDim(ncid,'eta_v',M);
%netcdf.defDim(ncid,'xi_rho',Lp);
%netcdf.defDim(ncid,'eta_rho',Mp);
%netcdf.defDim(ncid,'xi_psi',L);
%netcdf.defDim(ncid,'eta_psi',M);
%netcdf.defDim(ncid,'sms_time',length(smst));
%netcdf.defDim(ncid,'slp_time',length(slpt));

%nw('tide_period')=Ntides;

%
%  Create variables and attributes
%
nw{'sms_time'} = ncdouble('sms_time');
nw{'sms_time'}.long_name = ncchar('surface momentum stress time');
nw{'sms_time'}.long_name = 'surface momentum stress time';
nw{'sms_time'}.units = ncchar('days');
nw{'sms_time'}.units = 'days';
nw{'sms_time'}.cycle_length = smsc;

% smst_var=netcdf.defVar(ncid,'sms_time','double',[eta_u, xi_u]);
% netcdf.endDef(ncid)

nw{'slp_time'} = ncdouble('slp_time');
nw{'slp_time'}.long_name = ncchar('sealevel pressure time');
nw{'slp_time'}.long_name = 'sealevel pressure time';
nw{'slp_time'}.units = ncchar('days');
nw{'slp_time'}.units = 'days';
nw{'slp_time'}.cycle_length = slpc;


nw{'sustr'} = ncdouble('sms_time', 'eta_u', 'xi_u');
nw{'sustr'}.long_name = ncchar('surface u-momentum stress');
nw{'sustr'}.long_name = 'surface u-momentum stress';
nw{'sustr'}.units = ncchar('Newton meter-2');
nw{'sustr'}.units = 'Newton meter-2';

nw{'svstr'} = ncdouble('sms_time', 'eta_v', 'xi_v');
nw{'svstr'}.long_name = ncchar('surface v-momentum stress');
nw{'svstr'}.long_name = 'surface v-momentum stress';
nw{'svstr'}.units = ncchar('Newton meter-2');
nw{'svstr'}.units = 'Newton meter-2';

nw{'slp'} = ncdouble('slp_time', 'eta_rho', 'xi_rho');
nw{'slp'}.long_name = ncchar('sealevel pressure');
nw{'slp'}.long_name = 'sealevel pressure';
nw{'slp'}.units = ncchar('mb');
nw{'slp'}.units = 'mb';


%
%  Add variables and attributes
%
nw{'tide_period'} = ncdouble('tide_period');
nw{'tide_period'}.long_name = ncchar('Tide angular period');
nw{'tide_period'}.long_name = 'Tide angular period';
nw{'tide_period'}.units = ncchar('Hours');
nw{'tide_period'}.units = 'Hours';

nw{'tide_Ephase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Ephase'}.long_name = ncchar('Tidal elevation phase angle');
nw{'tide_Ephase'}.long_name = 'Tidal elevation phase angle';
nw{'tide_Ephase'}.units = ncchar('Degrees');
nw{'tide_Ephase'}.units = 'Degrees';

nw{'tide_Eamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Eamp'}.long_name = ncchar('Tidal elevation amplitude');
nw{'tide_Eamp'}.long_name = 'Tidal elevation amplitude';
nw{'tide_Eamp'}.units = ncchar('Meter');
nw{'tide_Eamp'}.units = 'Meter';

nw{'tide_Cmin'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cmin'}.long_name = ncchar('Tidal current ellipse semi-minor axis');
nw{'tide_Cmin'}.long_name = 'Tidal current ellipse semi-minor axis';
nw{'tide_Cmin'}.units = ncchar('Meter second-1');
nw{'tide_Cmin'}.units = 'Meter second-1';

nw{'tide_Cmax'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cmax'}.long_name = ncchar('Tidal current, ellipse semi-major axis');
nw{'tide_Cmax'}.long_name = 'Tidal current, ellipse semi-major axis';
nw{'tide_Cmax'}.units = ncchar('Meter second-1');
nw{'tide_Cmax'}.units = 'Meter second-1';

nw{'tide_Cangle'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cangle'}.long_name = ncchar('Tidal current inclination angle');
nw{'tide_Cangle'}.long_name = 'Tidal current inclination angle';
nw{'tide_Cangle'}.units = ncchar('Degrees between semi-major axis and East');
nw{'tide_Cangle'}.units = 'Degrees between semi-major axis and East';

nw{'tide_Cphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nw{'tide_Cphase'}.long_name = ncchar('Tidal current phase angle');
nw{'tide_Cphase'}.long_name = 'Tidal current phase angle';
nw{'tide_Cphase'}.units = ncchar('Degrees');
nw{'tide_Cphase'}.units = 'Degrees';

endef(nw);

%
% Create global attributes
%

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.type = ncchar('ROMS forcing file');
nw.type = 'ROMS forcing file';

nw.start_tide_mjd = ymd;
nw.components = ncchar(components);
nw.components = components;

%
% Write time variables
%


nw{'sms_time'}(:) = smst;
nw{'slp_time'}(:) = slpt;
close(nw);


