function ASET_CreateTecplotFileBat(guiparams,savefile,handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the file of the bottom to Tecplot (R)

% by Dominguez Ruben, L FICH-UNL

%%Bottom file

% TecPlot Variable List
% +=======================================================================+
% |   NAME             |   DESCRIPTION                                    |
% +=======================================================================+
% |   BedDepth         |   Bed depth (m or ft)                            |
% |   Dist             |   Dist across XS, oriented looking u s-1 (m or ft)|
% |   BedElev          |   Bed elevation (m or ft)                        |
% +=======================================================================+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outfile = [savefile(1:end-4) '_TecplotBottom.dat'];

% Sort the Distances such that when plotting in 2D (Dist. vs. Depth), 
% you are looking upstream into the transect
Dist = sort(guiparams.Array.Dist,2,'descend');

% Build tecplot data matrix
tecdata = [guiparams.xutm(:,1) guiparams.yutm(:,1) Dist(1,:)' guiparams.Array.Bed' ...
    guiparams.Array.Bed'];

%size(tecdata)
% Replace NaNs with a no data numeric value
nodata = -999;
n = find(isnan(tecdata));
tecdata(n) = nodata;

% Print out a TECPLOT FILE
fid = fopen(outfile,'w');
fprintf(fid, 'TITLE     = "AVEXSEC_TECOUT"\n');
if handles.units==0%SI
    fprintf(fid, 'VARIABLES = "X [m]"\n');
    fprintf(fid, '"Y [m]"\n');
    fprintf(fid, '"Distance [m]"\n');
    fprintf(fid, '"Depth_U [m]"\n');
    fprintf(fid, '"Depth_Ms2 [m]"\n');
elseif handles.units==1%english
    fprintf(fid, 'VARIABLES = "X [ft]"\n');
    fprintf(fid, '"Y [ft]"\n');
    fprintf(fid, '"Distance [ft]"\n');
    fprintf(fid, '"Depth_U [ft]"\n');
    fprintf(fid, '"Depth_Ms2 [ft]"\n');
end
fprintf(fid, 'ZONE T="ZONE 1"\n');
fprintf(fid, ' I=%d  J=1',size(tecdata,1));
fprintf(fid, '  K=1');
fprintf(fid, ' F=POINT\n');
fprintf(fid, 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n');
for m = 1:size(tecdata,1)
    fprintf(fid,'%10.8f %10.8f %10.8f %6.8f %10.8f\n',tecdata(m,:));
end
fclose(fid);
format short