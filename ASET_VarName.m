%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function describe the variables used in ASET

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Struct
% PARAMETERS            DESCRIPTION
%% ADCP                 ADCP data
% ADCP.Inst             Data structure for variables related to the instrument
% ADCP.rssi_beam1       RSSI Slope beam1
% ADCP.rssi_beam2       RSSI Slope beam2
% ADCP.rssi_beam3       RSSI Slope beam3
% ADCP.rssi_beam4       RSSI Slope beam4
% ADCP.R.draft          Value of draft
% ADCP.R.multi          YES or NO enable multiselect option
% ADCP.kcEr             kcEr calibration value
% ADCP.draft            YES or NO data of draft value
% ADCP.noenablekcEr     0 or 1 data of kcEr(Calibration or calculate)
% ADCP.first            Is the first or no calculate (0 or 1)
% ADCP.length           Number of file
% ADCP.Calibration      Calibration or calculate (0 or 1)
% ADCP.ccoef            C coefficient depending od ADCP

%% Array                Array final 
% Array.indexCssbottom  Start and end index of bottom region for Ms2 or Css
% Array.indexCssmeas    Start and end index of measured region for Ms2 or Css
% Array.indexCsssurf    Start and end index of surface region for Ms2 or Css
% Array.Css             Array of Css [kg/l]
% Array.Depth_Css       Array of depth corresponding to Css value [m]
% Array.indexUbottom    Stard and end index of bottom region in Array.u
% Array.indexCssmeas    Start and end index of measured region in Array.u
% Array.indexCsssurf    Start and end index of surface region in Array.u
% Array.u               Array of u [cm/s]
% Array.Depth_u         Array of depth correponding to u [m]
% Array.Bed             Bed elevation [m]
% Array.Dist            Distance travelled [m]
% Array.Gss             Array of Gss [kg/s] (coarse material transport)
% Array.Gw              Array of Gw [kg/s] (washload)
% Array.Depth_Gss       Depth of Gss values [m]
% Array.Depth_Gw        Depth of Gw values [m]
% Array.Css_ave         Adimensional Css (Css/max(Css)) for ensemble
% Array.u_ave           Adimensional u (u/max(u)) for ensemble
% Array.Depth_u_ave     Adimensional depth (Array.Depth_u/Array.Bed) f Ens.
% Array.Depth_Css_ave   Adimensional depth (Array.Depth_Css/Arra.Bed)
% Array.Vel20           Velocity corresponding to 20% of depth for graph
% Array.Vel40           Velocity corresponding to 40% of depth for graph
% Array.Vel60           Velocity corresponding to 60% of depth for graph
% Array.Vel80           Velocity corresponding to 80% of depth for graph
% Array.Velmin          Minimum of Array.Vel20,40,60 and 80%
% Array.Velmax          Maximum of Array.Vel20,40,60 and 80%
% Array.Css20           Css corresponding to 20% of depth for graph
% Array.Css40           Css corresponding to 40% of depth for graph
% Array.Css60           Css corresponding to 60% of depth for graph
% Array.Css80           Css corresponding to 80% of depth for graph
% Array.Cssmin          Minimum of Array.Css20,40,60 and 80%
% Array.Cssmax          Maximum of Array.Css20,40,60 and 80%
% Array.GssT            Gss total, the sum of different zones [kg/s]
% Array.GwT             Gw total, the sum of different zones [kg/s]

%% C_back               Variables of size array for the Css 
% C_back.zfrombottom    Depth cell from bottom for all regions to Css [m]
% C_back.zbottomback    Depth cell from bottom for bottom region to Css [m]
% C_back.zbottom        Depth of first cell from bottom [m]
% C_back.zpredbottomback Depth for predictable variable [m]
% C_back.cellbottomBack Distance between cell in the bottom region [m]
% C_back.zsurfback      Depth cell from bottom fot surface region to Css
% C_back.zsurf          Depth of last cell from bottom [m]
% C_back.zpredsurfback  Depth for predictable variable [m]
% C_back.cellsurfBack   Distance between cell in the surface region [m]

%% C_vel               Variables of size array for Velocity
% C_vel.zfrombottom    Depth cell from bottom for all regions to Vel [m]
% C_vel.zbottomvel     Depth cell from bottom for bottom region to Vel [m]
% C_vel.zbottom        Depth of first cell from bottom [m]
% C_vel.zpredbottomvel Depth for predictable variable [m]
% C_vel.cellbottomvel  Distance between cell in the bottom region [m]
% C_vel.zsurfvel       Depth cell from bottom fot surface region to Vel [m]
% C_vel.zsurf          Depth of last cell from bottom [m]
% C_vel.zpredsurfvel   Depth for predictable variable [m]
% C_vel.cellsurfvel    Distance between cell in the surface region [m]

%% Css                 Coarse suspended sediment concentration [kg/l]
% Css                  Array of coarse suspended sediment concetration [kg/l]

%% CssExtra            Data obtain for the extrapolation of the Css
% CssExtra.pred        Concentration for predictable equation fit
% CssExtra.bottom      Concentration extrapolated for bottom region
% CssExtra.surf        Concetration extrapolated for surface region
% CssExtra.nCa         Concetration in depth reference (0.5h)
% CssExtra.NRouse      Rouse Number
% CssExtra.a           Reference depth (0.5h)
% CssExtra.zpred       Depth from surface [m]
% CssExtra.acoezpred   Depth fix from Rouse equation. See Rouse (1957)
% CssExtra.acoezbottom Depth fix from Rouse equation to bottom region
% CssExtra.acoezsurf   Depth fix from Rouse equation to surface region
% CssExtra.rcuad       rcuad statistics parameter (fit parameter)
% CssExtra.ustarrouse  ustar from Rouse equation

%% CssExtraM           Extrapolation method choose for Css

%% Dis                 Discharge liquid struct
% Dis.Meas             Array of discharge liquid of measured region
% Dis.MeasT            Discharge liquid total, the sum of the measured region
% Dis.bootm            Array of discharge liquid of bottom region
% Dis.bottomT          Discharge liquid, the sum of the bottom region
% Dis.surf             Array of dicharge liquid of surface region
% Dis.surfT            Discharge liquid, the sum of the surface region
% Dis.Total            Discharge total, the sum of the all regions

%% Meas                Data of Gss and Gw for Measured region
% Meas.coarse          Gss for the measured region per each cell [kg/s] 
% Meas.fine            Gw for measured region per each cell [kg/s] 
% Meas.GssT            Total Gss for measured region [kg/s]
% Meas.GwT             Total Gw for measured region [kg/s]

%% NoMeas              Data of Gss and Gww for unmeasured regions 
% NoMeas.Surf.coarse   Gss for surface region per each cell [kg/s]
% NoMeas.Surf.fine     Gw for surface region per each cell [kg/s]
% NoMeas.Surf.GssT     Total Gss for surface region [kg/s]
% NoMeas.Surf.GwT      Gw for surface region per each cell [kg/s]
% NoMeas.Bottom.coarse Gss for bottom region per each cell [kg/s]
% NoMeas.Bottom.fine   Gw for bottom region per each cell [kg/s]
% NoMeas.Bottom.GssT   Total Gss for bottom region [kg/s]
% NoMeas.Bottom.GwT    Gw for bottom region per each cell [kg/s]

%% P                   Parameters for calculate
%% P.Ac                Acoustic Parameters
% P.Ac.FREQ            ADCP Frequency
% P.Ac.Draft_ADCP      Draft distance [m or ft]
% P.Ac.ALPHAt          ADD
% P.Ac.beam1           RSSI Slope beam1
% P.Ac.beam2           RSSI Slope beam2
% P.Ac.beam3           RSSI Slope beam3
% P.Ac.beam4           RSSI Slope beam4
% P.Ac.mean_beam       Average value of RSSI Slope
% P.Ac.kcEr            Calibration value
% P.Ac.ccoef           Coeficcient depend of the ADCP
% P.Ac.vson            Sound velocity 1500 [m/s]
% P.Ac.lambda          Wave longitud
% P.Ac.k               k coefficient
% P.Ac.ALPHAw          water attenuation
% P.Ac.gama            frequency seconds
% P.Ac.sf              ADD
% P.Ac.sc              ADD
% P.Ac.tauf            ADD
% P.Ac.tauc            ADD
% P.Ac.ALPHAf1         1 part of equation of ALPHAf (fine material)
% P.Ac.ALPHAf2         2 part of equation of ALPHAf (fine material)
% P.Ac.ALPHAf3         3 part of equation of ALPHAf (fine material)
% P.Ac.ALPHAf          Equation of ALPHAf (fine material)
% P.Ac.ALPHAg1         1 part of equation of ALPHAg (coarse material)
% P.Ac.ALPHAg2         2 part of equation of ALPHAg (coarse material)
% P.Ac.ALPHAg3         3 part of equation of ALPHAg (coarse material)
% P.Ac.ALPHAg          Equation of ALPHAg (coarse material)
% P.Ac.r               Projection 20 degrees beam depth [m]
% P.Ac.R               Correction phi [m]
% P.Ac.psi             Near field correction
% P.Ac.ALPHAgc         ADD
%% P.El                Electric Parameters
% P.El.TC              Transmit current
% P.El.TV              Transmit voltage
% P.El.L               Longitude
% P.ElPt               ADD
%% P.Ph                Physics Parameters 
% P.Ph.Temp            Average Temperature [degrees]
% P.Ph.Rhos            Sediment density (cuarzo)[mg/l]
% P.Ph.x               ADD
% P.Ph.f               ADD
% P.Ph.ks              ks**2
% P.Ph.Rhof            Water density [mg/l]
% P.Ph.nu              Kinetic water viscosity [m2/s]  
% P.Ph.s               RHO relation especific value

%% Pos                 GPS data. the coordinate X and Y
% Pos.xutm             X coordinate position in UTM [m]
% Pos.yutm             Y coordinate position in UTM [m]

%% QEdge               Data of Gss, Gw and Q close to the bank
% QEdge.Gss.left       Gss close to the left [kg/s]. See Teledyne 2005
% QEdge.Gss.left       Gss close to the left [kg/s]. See Teledyne 2005
% QEdge.Gw.left        Gw close to the left [kg/s]. See Teledyne 2005
% QEdge.Gw.left        Gw close to the left [kg/s]. See Teledyne 2005
% QEdge.Q.left         Q close to the left [kg/s]. See Teledyne 2005
% QEdge.Q.left         Q close to the left [kg/s]. See Teledyne 2005

%% S                   Sediment information input by ASET user's
% S.Df                 Diameter of fne material [m]
% S.Dc                 Diameter of coarse material [m]
% S.Csf                Concetration of fine material [mg/l]
% S.Css                Concetration os coarse material [mg/l]. Default data

%% St                  Struct input dimension
% St.nvel              Rows number of velocity array
% St.mvel              Column number of velocity array
% St.nback             Rows number of intensity array
% St.mback             Column number of intensity array
% St.cell              Cell size [m]

%% V                   Data read by the input file
% V.mcsBack            Array of Intensity [db]
% V.mcsMag             Velocity magnitude [cm/s]
% V.mcsDepth           Depth from surface [m]
% V.mcsBed             Bed [m]
% V.mcsDist            Distance from left bank [m]
% V.navref             Navigation reference
% V.ensDeltaTime       Time traveled [s]
% V.u                  Norma velocity of path [cm/s]

%% VelExtra            Data obtain for the extrapolation of the Velocity
% VelExtra.pred        Concentration for predictable equation fit
% VelExtra.bottom      Concentration extrapolated for bottom region
% VelExtra.surf        Concetration extrapolated for surface region
% VelExtra.zpred       Depth from surface [m]
% VelExtra.rcuad       rcuad statistics parameter (fit parameter)
% VelExtra.ustarwl     ustar from Law of th wall equation

%% VelExtraM           Extrapolation method choose for Vel


