function [R,Inst,Cfg,Sensor]=ASET_ReadPd0(infile,draft,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function reads the the pd0s files and obtained the intesity of signal 
% or electric variables for Calculation Module or Calibration Module.
%Dominguez Ruben, 12/2017, FICH-UNL, Argentina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M==1%calculate of Qss read variables

    [Inst, Cfg, Sensor]=...
    pd0ReadMody(infile);

    %Control if Freq diferent to 1200 or 600kHz
     if Inst.freq(1)==1200 | Inst.freq(1)==600
     else
         wrong_data=warndlg('The ADCP frequency is not applicable into ASET.',...
          'WARNING');
     end

     R=[];

elseif M==2 %Calibration read variables
%Start code   
[~, Inst, Cfg, Sensor, ~, Wt, Bt, ~, ~, Surface, AutoMode]=...
    pd0Read(infile);
% 
    R.draft= draft;%
%     

    % Retrieve data to compute sidelobe
    % ---------------------------------
    lag=Cfg.lag_cm'./100;
    pulseLen=Cfg.xmitPulse_cm'./100;
    regCellSize=Cfg.ws_cm'./100;
    regCellSize(regCellSize==0)=nan;
    obj.beamDepths=Bt.depth_m;
    obj.beamDepths(obj.beamDepths<0.01)=nan; % Screen bad depths reported as zero
    cell1Dist=Cfg.distBin1_cm'./100;
    obj.numEns=size(Wt.vel_mps(:,:,1),2);

    % Beam angle is used to ID RiverRay data with variable modes and
    % lags
    % --------------------------------------------------------------
    if Inst.beamAng(1)>21
        lag(AutoMode.Beam1.mode<2 | AutoMode.Beam1.mode>4)=0;
    end
    
    % surf* data are to accomodate RiverRay. readpd0rr2 sets these
    % values to nan when reading Rio Grande or StreamPro data
    % -----------------------------------------------------------
    surfCells_idx=find(isnan(Surface.no_cells));
    noSurfCells=Surface.no_cells';
    noSurfCells(surfCells_idx)=0;
    maxSurfCells=nanmax(noSurfCells);
    surfCellSize=Surface.cell_size_cm'./100;
    surfCell1Dist=Surface.dist_bin1_cm'./100;
    numRegCells=size(Wt.vel_mps,1);
    obj.maxCells=maxSurfCells+numRegCells;
    
    % Beam angle is used to ID RiverRay data with variable modes and
    % lags
    % --------------------------------------------------------------
    %
    % Compute side lobe interference limit
    % ------------------------------------
    lagEffect_m=(lag+pulseLen+regCellSize)./2;
    depthmin=nanmin(obj.beamDepths);
    lastcell=depthmin.*cosd(Inst.beamAng(1)')-(lagEffect_m);

    obj.numCells=max([floor(((lastcell-cell1Dist)./regCellSize)+1); zeros(size(lastcell))],[],1);
    obj.numCells(obj.numCells>numRegCells)=numRegCells;
    if nanmax(noSurfCells)>0
        obj.numCells=obj.numCells+noSurfCells;
    end
    %
    % Create matrix with only cells above sidelobe
    % ---------------------------------------------
    obj.cellsAboveSL=ones(obj.maxCells,obj.numEns);
    for jj=1:obj.numEns
        for l=obj.numCells(jj)+1:obj.maxCells
            obj.cellsAboveSL(l,jj)=nan;
        end
    end


    % Develop distance to center of top cell for each ensemble
    % --------------------------------------------------------
    depth=obj.beamDepths';
    if nanmax(noSurfCells)>0
        dist_cell1_m=surfCell1Dist;
        dist_cell1_m(surfCells_idx)=cell1Dist(surfCells_idx);
    else
        dist_cell1_m=cell1Dist;
    end

    % Combine cell size and cell range from transducer for both
    % surface and regular cells
    % ---------------------------------------------------
    obj.cellDepth=nan(obj.maxCells,obj.numEns);
    cellSizeAll=nan(obj.maxCells,obj.numEns);
    for ii=1:obj.numEns
        if nanmax(noSurfCells)>0
            numRegCells=obj.maxCells-noSurfCells(ii);
        else
            numRegCells=obj.maxCells;
        end
        %
        % Surface cell are present
        % ------------------------
        if nanmax(noSurfCells)>0 && noSurfCells(ii)>0
            obj.cellDepth(1:noSurfCells(ii),ii)=dist_cell1_m(ii)+(0:surfCellSize(ii):(noSurfCells(ii)-1)*surfCellSize(ii))';
            obj.cellDepth(noSurfCells(ii)+1:end,ii)=obj.cellDepth(noSurfCells(ii),ii)+0.5*surfCellSize(ii)+0.5.*regCellSize(ii)+(0:regCellSize(ii):(numRegCells-1)*regCellSize(ii));
            cellSizeAll(1:noSurfCells(ii),ii)=repmat(surfCellSize(ii),noSurfCells(ii),1);
            cellSizeAll(noSurfCells(ii)+1:end,ii)=repmat(regCellSize(ii),numRegCells,1);
            %
            % No surface cells
            % ----------------
        else
            %
            obj.cellDepth(1:numRegCells,ii)=dist_cell1_m(ii)+[0:regCellSize(ii):(numRegCells-1)*regCellSize(ii)];
            cellSizeAll(1:end,ii)=repmat(regCellSize(ii),numRegCells,1);
        end
    end

        % Compute weighted mean depth
            % ---------------------------
            w=1-depth./repmat(nansum(depth,2),1,4);
            obj.depthEns=nansum((depth.*w)./repmat(nansum(w,2),1,4),2)';

    % Compute cell depths from water surface
    % --------------------------------------
     R.Depth_m=obj.cellDepth+R.draft;

     R.Bed=obj.depthEns+R.draft;


    %%
    %Backscatter
    R.Back.Beam1(:,:,1)=Wt.rssi(:,:,1).*(isfinite(obj.cellsAboveSL));
    R.Back.Beam2(:,:,1)=Wt.rssi(:,:,2).*(isfinite(obj.cellsAboveSL));
    R.Back.Beam3(:,:,1)=Wt.rssi(:,:,3).*(isfinite(obj.cellsAboveSL));
    R.Back.Beam4(:,:,1)=Wt.rssi(:,:,4).*(isfinite(obj.cellsAboveSL));

    clear obj wVel wVel_reGGA wVel_reVTG vMag_norm wVelx wVely bVelx bVely
end