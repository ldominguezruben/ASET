function [VarInp]=ASET_ReadInputFile(path,file,lastPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function read the ASCII or MAT file dependeing the mpodule read the
% file format. Also read the pd0 file .

% by Dominguez Ruben, L., FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read the format of file
format=mat2str(file(end-2:end));
if format(2)=='T'%Read ASCII file

    VarInp.readformatfile=2;

    % screenData 0 leaves bad data as -32768, 1 converts to NaN
    screenData=1;
    ignoreBS = 0;

    % read file
    [A]=tfile(fullfile(path,file),screenData,ignoreBS);           

    % Variables predefined
    VarInp.V.mcsBack=nanmean(A.Wat.backscatter,3);

    % Replace NaNs with a no data numeric value
    nodata = nan;
    n = find(VarInp.V.mcsBack==255);
    VarInp.V.mcsBack(n) = nodata;

    VarInp.V.mcsMag=A.Wat.vMag;
    VarInp.V.mcsDepth=A.Wat.binDepth;
    A.Nav.depth(A.Nav.depth<-32678)=nan;
    VarInp.V.mcsBed=nanmean(A.Nav.depth');
    VarInp.V.mcsDist=A.Nav.length';
    VarInp.V.navref=A.Sup.vRef;
    VarInp.R.draft = double(A.Sup.draft_cm)./100;
    VarInp.A=A;


   %Find file and determinate left and right distance for discharge 
    VarInp.Edge.leftdistance=A.Q.startDist(1,1);%
    VarInp.Edge.rightdistance=A.Q.endDist(1,1);%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Distance

    for i=1:length(A.Nav.length)
        if i==1
             VarInp.V.ensDeltaTime(i)=0;
        else
             VarInp.V.ensDeltaTime(i)=A.Sup.timeElapsed_sec(i,1)-...
                 A.Sup.timeElapsed_sec(i-1,1);%Deltat
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Position only ASCII

    if isnan(A.Nav.lat_deg(1,1))
       VarInp.xutm=A.Nav.totDistEast;
       VarInp.yutm=A.Nav.totDistNorth;
    else
       [VarInp.xutm,VarInp.yutm]=deg2utm(A.Nav.lat_deg,A.Nav.long_deg);%m
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the cross product used in the other computations
% -------------------------------------------------------------------------

    wVelx=A.Wat.vEast;
    wVely=A.Wat.vNorth;
    bVelx=A.Nav.bvEast';
    bVely=A.Nav.bvNorth';

    VarInp.V.u = ((wVelx.*repmat(bVely,size(wVelx,1),1))-...
    (wVely.*repmat(bVelx,size(wVely,1),1)))./100;%cm/s

    if nanmean(nanmean(VarInp.V.u)) < 0 % the wrong start bank was used so the velocities have the wrong sign    
       VarInp.V.u = -VarInp.V.u;
    end

    %General Parameters
    VarInp.FileName = file;
    VarInp.PathName = path;
    VarInp.SaveFile = fullfile(path,file);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ MAT FILE (from VMT USGS)
elseif format(2)=='m'
 %Read .mat file        
    VarInp.readformatfile=1;
    % Make sure the selected file is a valid file:
    % --------------------------------------------
    varnames = fieldnames(load(fullfile(path,file)));
    if isequal(sort(varnames),{'A' 'Map' 'V' 'z'}')
        load(fullfile(path,file))
        VarInp.z = z;
        VarInp.A = A;
        VarInp.V = V;       
        VarInp.V.navref=A(1).Sup.vRef;
         
      %Find file and determinate left and right distance
        VarInp.Edge.leftdistance = A(1).Q.startDist(1,1);%
        VarInp.Edge.rightdistance = A(1).Q.endDist(1,1);%
        VarInp.R.draft = double(A(1).Sup.draft_cm)./100;
        VarInp.R.multi = 0; 
        VarInp.xutm = V.mcsX';
        VarInp.yutm = V.mcsY';
                
    else % Not a valid file
        errordlg('The selected file is not a valid ADCP data MAT file.', ...
            'Invalid File...')
    end
    
        %General Parameters
    VarInp.FileName = file;
    VarInp.PathName = path;
    VarInp.SaveFile = fullfile(path,file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ PD0 file
message=['Select PD0 File corresponding to ' file];

%Start code
if lastPath == 0
[filePD0,pathPD0] = uigetfile({'*.PD0',...
    'ASET Files (*.PD0)';'*.*',  'All Files (*.*)'}, message, path);
else %remember the lastpath
[filePD0,pathPD0] = uigetfile({'*.PD0',...
    'ASET Files (*.PD0)';'*.*',  'All Files (*.*)'}, message, lastPath);
end

%Parameters 
M=1;
draft=[];
[~,VarInp.Inst,VarInp.Cfg,VarInp.Sensor]=ASET_ReadPd0(fullfile(pathPD0,filePD0),draft,M);

%General Parameters
VarInp.FileNamePD0 = filePD0;
VarInp.PathNamePD0 = pathPD0;
VarInp.SaveFilePD0 = fullfile(pathPD0,filePD0);
    
VarInp.first = 0;%to ADCP data GUI not the first read file
VarInp.Calibration = 0;% to ADCP data GUI not calibration


