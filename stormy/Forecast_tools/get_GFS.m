function [t,tx,ty,pressfc]=...
         get_GFS(fname,mask,tndx,jrange,i1min,i1max,...
	 i2min,i2max,i3min,i3max,missvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Download one full subset of GFS for ROMS bulk for 1 time step
% Put them in the ROMS units
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
trange=['[',num2str(min(tndx)),':',num2str(max(tndx)),']'];
%
% Get GFS variables for 1 time step
%
disp(' ')
disp('====================================================')
%disp('time...')
%disp(['TNDX=',num2str(tndx)])
%disp(['TRANGE=',num2str(trange)])
t=readdap(fname,'time',trange);
disp(['TRANGE=',num2str(trange)])
disp(['GFS raw time=',sprintf('%5.3f',t)])
%t=t+364.75 % put it in "matlab" time. PM
t=t+365; % put it in "matlab" time. GC
disp(['GFS: ',datestr(t)])
disp('====================================================')

%disp('ty...')
ty=mask.*getdap('',fname,'vflxsfc',trange,'',jrange,...
                i1min,i1max,i2min,i2max,i3min,i3max);
ty(abs(ty)>=missvalue)=NaN;

%disp('tx...')
tx=mask.*getdap('',fname,'uflxsfc',trange,'',jrange,...
                i1min,i1max,i2min,i2max,i3min,i3max);
tx(abs(tx)>=missvalue)=NaN;


%disp('tair...')
pressfc=mask.*getdap('',fname,'pressfc',trange,'',jrange,...
                i1min,i1max,i2min,i2max,i3min,i3max);
pressfc(abs(pressfc)>=missvalue)=NaN;

%
% Transform the variables
%
%
% 1: Surface Pressure: Convert from Pa to Millibar
%
pressfc=pressfc/100;
%
% % 7: Compute the stress following large and pond
% %
% [Cd,uu]=cdnlp(wspd,10.);
% rhoa=air_dens(tair,rhum*100);
% tx=Cd.*rhoa.*u.*wspd;
% ty=Cd.*rhoa.*v.*wspd;
%
return
