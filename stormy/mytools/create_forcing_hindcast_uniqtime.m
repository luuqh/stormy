function  create_forcing_hindcast_uniqtime (frcname,grdname,title,smst,smsc)
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
result = redef(nw);

%
%  Create dimensions
%



nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('sms_time') = netcdf.getConstant('NC_UNLIMITED');
%nw('sms_time') = length(smst);
%nw('slp_time') = length(slpt);
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

%nw{'slp_time'} = ncdouble('slp_time');
%nw{'slp_time'}.long_name = ncchar('sealevel pressure time');
%nw{'slp_time'}.long_name = 'sealevel pressure time';
%nw{'slp_time'}.units = ncchar('days');
%nw{'slp_time'}.units = 'days';
%nw{'slp_time'}.cycle_length = slpc;


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

%nw{'slp'} = ncdouble('slp_time', 'eta_rho', 'xi_rho');
nw{'slp'} = ncdouble('sms_time', 'eta_rho', 'xi_rho');
nw{'slp'}.long_name = ncchar('sealevel pressure');
nw{'slp'}.long_name = 'sealevel pressure';
nw{'slp'}.units = ncchar('mb');
nw{'slp'}.units = 'mb';

result = endef(nw);

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

%
% Write time variables
%

nw{'sms_time'}(1:length(smst)) = smst;
%nw{'sms_time'}(:) = smst;
%nw{'slp_time'}(:) = slpt;

close(nw);


