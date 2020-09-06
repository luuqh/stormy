

% -------------------------------------------------------------------------
% Matlab code developed at National University of Singapore TMSI/PORL
% on 2013/04/09
% -------------------------------------------------------------------------
%     Supporter: Pavel Tkalich
%     Coder: Luu Quang Hung
%     Email: luuquanghung@gmail.com
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
% CREATE FORCING FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function  create_forcing_robust (frcname,grdname,time)


% --------------------------------------   
% get grid
% --------------------------------------   

nc = netcdf(grdname);
L = length(nc('xi_psi'));
M = length(nc('eta_psi'));
close(nc);
Lp = L+1;
Mp = M+1;
disp([L M Lp Mp])


% --------------------------------------   
% create dimensions
% --------------------------------------   

disp(['creating a new forcing file: ' frcname]);

%nw = netcdf(frcname,'write');
%result = redef(nw);
%nw('xi_u') = L;
%nw('eta_u') = Mp;
%nw('xi_v') = Lp;
%nw('eta_v') = M;
%nw('xi_rho') = Lp;
%nw('eta_rho') = Mp;
%nw('xi_psi') = L;
%nw('eta_psi') = M;
%nw('sms_time') = length(time);
%nw('slp_time') = length(time);


% --------------------------------------   
%  create variables and attributes
% --------------------------------------   


%nw{'sms_time'} = ncdouble('sms_time');
%nw{'sms_time'}.long_name = ncchar('surface momentum stress time');
%nw{'sms_time'}.long_name = 'surface momentum stress time';
%nw{'sms_time'}.units = ncchar('days');
%nw{'sms_time'}.units = 'days';
%nw{'sms_time'}.cycle_length = time;

%nw{'slp_time'} = ncdouble('slp_time');
%nw{'slp_time'}.long_name = ncchar('sealevel pressure time');
%nw{'slp_time'}.long_name = 'sealevel pressure time';
%nw{'slp_time'}.units = ncchar('days');
%nw{'slp_time'}.units = 'days';
%nw{'slp_time'}.cycle_length = time;

%nw{'sustr'} = ncdouble('sms_time', 'eta_u', 'xi_u');
%nw{'sustr'}.long_name = ncchar('surface u-momentum stress');
%nw{'sustr'}.long_name = 'surface u-momentum stress';
%nw{'sustr'}.units = ncchar('Newton meter-2');
%nw{'sustr'}.units = 'Newton meter-2';

%nw{'svstr'} = ncdouble('sms_time', 'eta_v', 'xi_v');
%nw{'svstr'}.long_name = ncchar('surface v-momentum stress');
%nw{'svstr'}.long_name = 'surface v-momentum stress';
%nw{'svstr'}.units = ncchar('Newton meter-2');
%nw{'svstr'}.units = 'Newton meter-2';

%nw{'slp'} = ncdouble('slp_time', 'eta_rho', 'xi_rho');
%nw{'slp'}.long_name = ncchar('sealevel pressure');
%nw{'slp'}.long_name = 'sealevel pressure';
%nw{'slp'}.units = ncchar('mb');
%nw{'slp'}.units = 'mb';

%result = endef(nw);


% --------------------------------------   
% create dimensions
% --------------------------------------   

ncid = netcdf.create (frcname, 'nc_write');
dimid_xi_u = netcdf.defDim(ncid,'xi_u',L);
dimid_eta_u = netcdf.defDim(ncid,'eta_u',Mp);
dimid_xi_v = netcdf.defDim(ncid,'xi_v',Lp);
dimid_eta_v = netcdf.defDim(ncid,'eta_v',M);
dimid_xi_rho = netcdf.defDim(ncid,'xi_rho',Lp);
dimid_eta_rho = netcdf.defDim(ncid,'eta_rho',Mp);
dimid_xi_psi = netcdf.defDim(ncid,'xi_psi',L);
dimid_eta_psi = netcdf.defDim(ncid,'eta_psi',M);
dimid_smst = netcdf.defDim(ncid,'sms_time',length(time));
dimid_slpt = netcdf.defDim(ncid,'slp_time',length(time));

% --------------------------------------   
%  create variables and attributes
% --------------------------------------   

varid_smst = netcdf.defVar(ncid,'sms_time','double',dimid_smst);
netcdf.putAtt(ncid,varid_smst,'long_name','surface momentum stress time');
netcdf.putAtt(ncid,varid_smst,'units','days');
netcdf.putAtt(ncid,varid_smst,'cycle_length',length(time));

varid_slpt = netcdf.defVar(ncid,'slp_time','double',dimid_slpt);
netcdf.putAtt(ncid,varid_slpt,'long_name','sealevel pressure time');
netcdf.putAtt(ncid,varid_slpt,'units','days');
netcdf.putAtt(ncid,varid_slpt,'cycle_length',length(time));

varid_sustr = netcdf.defVar(ncid,'sustr','double',[dimid_smst dimid_eta_u dimid_xi_u]);
netcdf.putAtt(ncid,varid_sustr,'long_name','surface u-momentum stress');
netcdf.putAtt(ncid,varid_sustr,'units','Newton meter-2');

varid_svstr = netcdf.defVar(ncid,'svstr','double',[dimid_smst dimid_eta_v dimid_xi_v]);
netcdf.putAtt(ncid,varid_svstr,'long_name','surface v-momentum stress');
netcdf.putAtt(ncid,varid_svstr,'units','Newton meter-2');

varid_slp = netcdf.defVar(ncid,'slp','double',[dimid_slpt dimid_eta_rho dimid_xi_rho]);
netcdf.putAtt(ncid,varid_slp,'long_name','sealevel pressure');
netcdf.putAtt(ncid,varid_slp,'units','mb');

% --------------------------------------   
%  assign time
% --------------------------------------   

netcdf.endDef(ncid)
netcdf.putVar(ncid,varid_smst,time);
netcdf.putVar(ncid,varid_slpt,time);
netcdf.close(ncid) 





