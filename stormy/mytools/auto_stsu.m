

% -------------------------------------------------------------------------
% Matlab code developed at National University of Singapore TMSI/PORL
% on 2012/11/06
% -------------------------------------------------------------------------
%     Leader: Pavel Tkalich
%     Coder: Luu Quang Hung
%     Email: luuquanghung@gmail.com
% -------------------------------------------------------------------------


function auto_stsu (frcname,id,torig,trun,outpath,inpath)


idstr = num2str(id);
timestr = num2str(trun);
hisstr = ['roms_his_' idstr '.nc'];
avgstr = ['roms_avg_' idstr '.nc'];
inistr = ['roms_ini_' idstr '.nc'];
rststr = ['roms_rst_' idstr '.nc'];


% --------------------------------------    
% generate various roms.in
% --------------------------------------
    
froutname = ['roms' idstr '.in'];
frin = fopen('roms.in','rt');
frout = fopen([outpath froutname],'wt');
while ~feof(frin)
    s = fgetl(frin);
    s = strrep(s, 'AAAAAAAA', timestr);
    s = strrep(s, 'BBBBBBBB', frcname);
    s = strrep(s, 'CCCCCCCC', hisstr);
    s = strrep(s, 'DDDDDDDD', avgstr);
    s = strrep(s, 'EEEEEEEE', inistr);
    s = strrep(s, 'FFFFFFFF', rststr);
    fprintf(frout,'%s\n',s);
    % disp(s);
end
fclose(frin);
fclose(frout);


% --------------------------------------    
% generate various roms.qsub
% --------------------------------------
    
fqin = fopen('roms.qsub','rt');
fqout = fopen([outpath 'roms' idstr '.qsub'],'wt');
while ~feof(fqin)
    s = fgetl(fqin);
    % s = strrep(s, 'AAAAAAAA', initstr);
    s = strrep(s, 'BBBBBBBB', froutname);
    s = strrep(s, 'CCCCCCCC', ['status' idstr '.lst']);
    s = strrep(s, 'DDDDDDDD', ['std' idstr '.out']);
    s = strrep(s, 'EEEEEEEE', ['std' idstr '.err']);    
    fprintf(fqout,'%s\n',s);
    % disp(s);
end
fclose(fqin);
fclose(fqout);


% --------------------------------------    
% scrum time in init file
% --------------------------------------

copyfile([inpath 'roms_ini.nc'],[outpath inistr]);
nc = netcdf([outpath inistr],'w');
nc{'scrum_time'}(:)=torig*24*3600;
close(nc);







