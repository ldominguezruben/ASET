function ASET_CreateTecplotFile(guiparams,savefile,handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a file to read the following variables in Tecplot(R)

% by Dominguez Ruben,L. FICH-UNL

% TecPlot Variable List 7
% +=======================================================================+
% |   NAME             |   DESCRIPTION                                    |
% +=======================================================================+
% |   Dist             |   dist across XS(m or ft)                        |
% |   Depth_uMag       |   Velocity depth (m or ft)                       |
% |   UMag             |   vel magnitude (need better desc.) (cm/s)       |
% |   Depth_Ms2        |   Ms2 depth (m or ft)                            |
% |   Ms2              |   Sand Concentration (mg/l or tn/ft3)            |
% |   Depth_Qss        |   Ms2 depth (m or ft)                            |
% |   Qss              |   Sand Transport (kg/s or tn/s)                  |
% +=======================================================================+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

% Create block style matrix of all processed data
tecdata = [];

% Sort the Distances such that when plotting in 2D (Dist. vs. Depth), 
% you are looking upstream into the transect
Dist = sort(guiparams.Array.Dist,2,'descend');

% Build tecplot data matrix
for k = 1:size(guiparams.Array.Depth_uMag,2)
    for i = 1:size(guiparams.Array.Depth_uMag,1)
        tempvec = [guiparams.xutm(k,1) guiparams.yutm(k,1) Dist(1,k) ...
            guiparams.Array.Depth_uMag(i,k) guiparams.Array.uMag(i,k) ...
            guiparams.Array.Depth_Css(i,k) guiparams.Array.Css(i,k)*1000 ...
            guiparams.Array.Depth_Gss(i,k) guiparams.Array.Gss(i,k)];
        tecdata = [tecdata; tempvec];
    end
end

% Replace NaNs with a no data numeric value
nodata =-999;
n = find(isnan(tecdata));
tecdata(n) = nodata;

% Name of output file (needs to be modified to take handle args from GUI)
outfile = [savefile(1:end-4) '_Tecplot.dat'];


% Print out a TECPLOT FILE
fid = fopen(outfile,'w');
fprintf(fid, 'TITLE     = "AVEXSEC_TECOUT"\n');
if handles.units==0%SI
    fprintf(fid, 'VARIABLES = "X [m]"\n');
    fprintf(fid, '"Y [m]"\n');
    fprintf(fid, '"Distance [m]"\n');
    fprintf(fid, '"Depth_U [m]"\n');
    fprintf(fid, '"U [cm/s]"\n');
    fprintf(fid, '"Depth_Ms2 [m]"\n');
    fprintf(fid, '"Ms2 [mg/l]"\n');
    fprintf(fid, '"Depth_Qss [m]"\n');
    fprintf(fid, '"Qss [kg/s]"\n');
elseif handles.units==1%English units
    fprintf(fid, 'VARIABLES = "X [ft]"\n');
    fprintf(fid, '"Y [ft]"\n');
    fprintf(fid, '"Distance [ft]"\n');
    fprintf(fid, '"Depth_U [ft]"\n');
    fprintf(fid, '"U [fps]"\n');
    fprintf(fid, '"Depth_Ms2 [ft]"\n');
    fprintf(fid, '"Ms2 [t/ft3]"\n');
    fprintf(fid, '"Depth_Qss [ft]"\n');
    fprintf(fid, '"Qss [t/s]"\n');
end
fprintf(fid, 'ZONE T="Zone_1"\n');
fprintf(fid, ' I=%d  J=1',i);
fprintf(fid, '  K=%d',k);
fprintf(fid, ' F=POINT\n');
fprintf(fid, 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n');
for m = 1:size(tecdata,1)
    fprintf(fid,'%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n',tecdata(m,:));
end

fclose(fid);
format short

