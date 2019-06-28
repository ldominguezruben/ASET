function [Hdr, Inst, Cfg, Sensor, Gps, Wt, Bt, Nmea, Gps2, Surface, AutoMode]=readpd0rrss(fullName,varargin)
% Read RDI PD0 File Format (WinRiver)
% The function reads the RDI PD0 format raw data as recorded by WinRiver.
% It accounts for the special locations used by WinRiver to store GPS data.
% It should read standard PD0 formatted files but was written specifically
% for PD0 files created by WinRiver.
%
% Data are stored in Matlab data structures to provide easy and efficient
% transfer of data between functions. For numeric data the arrays are
% preallocated with "nan" (not a number) so that no value is assumed for bad
% or missing data. Data coded with -32768 are not stored and default to "nan".
% No filtering of the data is provided, although some units conversion is
% provided for consistency and convience. All data structures begin with an
% upper case letter. All variable names begin with a lower case letter. If
% a specific dimension is associated with a variable the variable name will
% include an underscore "_" followed by the dimension. The dimension
% abrevations used are as follows:
%   cm      - centimeters
%   deg     - degrees
%   degc    - degrees Celsius
%   dm      - decimeters
%   m       - meters
%   mm      - millimeters
%   mmps    - millimeters per second
%   mps     - meters per second
%   msc     - minutes seconds hundreths of a second
%   pascal  - pascals
%   ppt     - parts per thousand
%   sec     - seconds
%
% Data that vary with depth are stored so that the bin (depth) varies by
% row and the ensemble varies by column. Most other data are stored so that
% the ensemble varies by row. Data associated with various beams are stored
% with a 3rd dimension reflecting the beam number.
%
% David S. Mueller (dmueller@usgs.gov)
% U.S. Geological Survey
% Office of Surface Water
% July 11, 2005
%
% MODIFIED 6/23/2006 dsm
% 1) Fixed Gps. structure to be consisent with (iensembles,n).
% 2) Added check for checksum
%
% MODIFIED 6/28/2006 dsm
% 3) Added data type '5000' to allow reading of Streampro raw data
%
% MODIFIED 7/25/2006 dsm
% 4) Modified code to correct for bad checksum and to refind a valid
%    ensemble and read remainder of file.
%
% MODIFIED 4/18/2007 dsm
% 5) Added code to search for a '7F7F' leader_id before proceeding into the
% remainder of the code. Some files have the instrument defined or other
% information prior to the first '7F7F'.
%
% 6) Changed the way the end of file is determined.
%
% 7) Fixed the "other" case to handle undefined data types
%
% 8) Commented out data type '5000' since the code now handles unknown data
% types properly and we don't have the coding information for '5000'.
%
% MODIFIED 5/2/2007 dsm
% 9) Corrected use of bytePerEnsemble so that the correct value was used to
% compute the checksum.
%
% 10) Added support for the new 2022 data type with sub datatypes for the
% VTG, GGA, DBT, and HDG datatypes. These data are stored in Gps2.
%
% 11) Modified the way datatypes 21xx were read from the data file.
%
% MODIFIED 10/18/2007 dsm
% 12) Modified the way GPS data are recorded to include all decoded NMEA
% strings recorded during an ensemble. These data are stored in GPS2.
%

%
% MODIFIED 11/2/2007 dsm
% 14) Modified GPS 2022 data structure to read sub-data types 104-107.
%
% MODIFIED 11/14/2008 dsm
% 15) Change the method for computing number of ensembles used to
% preallocate arrays. Number of bytes per ensemble is now variable when GPS
% data is stored.
%
% 16) Added numerous comments and improved code effeciency based on MLint
% suggestions.
%
% MODIFIED 2/19/2010 dsm
% 17) Section by Section renumbers the ensembles so the current approach to
% preallocating the arrays does not work. If statements were modified where
% -32768 was screen to store nan.
%
% MODIFIED 7/12/2010 dsm
% 18) Fixed decimal on temperature data (divided by 100).
%
% 19) Added ability to decode RiverRay surface data.
%
% MODIFIED 6/8/2011 dsm
% 20) Added code to check for lost ensembles. If the ensembles numbers are
% not consecutive the ensemble number is computed and stored in Sensor.num
% and the iEns in incremented. Thus lost ensembles now store nan in all
% data field in the correct order represented by the ensemble number.
%
% MODIFIED 8/5/2011 dsm
% 21) Modified code to determine number of ensembles to check the checksum
% for a valid checksum. This eliminated the problem where there were
% multiple 7F7F combinations that really weren't the beginning of ensemble
% data.

%
% Open File
% =========
hwait=waitbar(0.0,'Reading PD0 File');


%
% Open File
% =========
    %
    % Get and display file information
    % --------------------------------
    FileInfo=dir(fullName);
info=dir(fullName)
%%
if info.bytes>0
    fid=fopen(fullName,'r','l');
    %
    % Read Selected Parameters
    % ========================
    %
    % Selected portions of the raw data file are read to determine the number
    % of ensembles, the number of bins, and the number of beams. These
    % variables are used to preallocate the arrays in the data structures used
    % to store the data in Matlab. Although preallocation is not required it
    % does improve the performance of Matlab.
    % -------------------------------------------------------------------------
    leader_id=dec2hex(fread(fid,1,'uint16'),4);
    if ~strcmp(leader_id,'7F7F')
        while ~strcmp(leader_id,'7F7F')
            fseek(fid,-1,'cof');
            leader_id=dec2hex(fread(fid,1,'uint16'),4);
        end
    end
    initialPos=ftell(fid)-2;
    %
    % Set the position in the file and read the number of bytes per ensemble.
    % -----------------------------------------------------------------------
    bytesPerEns=fread(fid,1,'uint16');
    %
    % Set the position in the file and read the number of data types in the
    % file.
    % ----------------------------------------------------------------------
    fseek(fid,1,'cof');
    nTypes=fread(fid,1,'uint8');
    offset=fread(fid,1,'uint16');
    
    % Set the position in the file and read the number of beams and the number
    % of bins.
    % ------------------------------------------------------------------------
    fseek(fid,initialPos+offset+8,'bof');
    ftell(fid);
    nBeams=fread(fid,1,'uint8');
    nBins=fread(fid,1,'uint8');
    %
    % Compute the number of ensembles
    % -------------------------------
    [nEnsembles]=numberofensembles(fid,FileInfo);
    waitbar(0.1);
    %
    % Define Data Storage Structure.
    % ==============================
    %
    % A data storage structure is used to provide an easy and efficient method
    % for passing large amounts of data between functions. The arrays within
    % the structure are preallocated for efficiency.
    % -------------------------------------------------------------------------
    %
    % Clear variables to be used.
    % ---------------------------
    clear Hdr Inst Cfg Sensor Gps Wt Bt  Nmea;
    %
    % Data structure for the Binary Header Data
    % -----------------------------------------
    Hdr=struct( 'bytesPerEns', zeros(nEnsembles,1),...
        'dataOffsets', zeros(nEnsembles,nTypes),...
        'nDataTypes', zeros(nEnsembles,1),...
        'dataOK',zeros(nEnsembles,1));
    %
    % Data structure for variables related to the instrument
    % ------------------------------------------------------
    Inst=struct('beamAng', zeros(nEnsembles,1),...
        'beams', zeros(nEnsembles,1),...
        'dataType', repmat(blanks(4),nEnsembles,1),...
        'firmVer', zeros(nEnsembles,1),...
        'freq', zeros(nEnsembles,1),...
        'pat', repmat(blanks(7),nEnsembles,1),...
        'resRDI', zeros(1),...
        'sensorCfg', nan(nEnsembles,1),...
        'xducer', repmat(blanks(12),nEnsembles,1));
    %
    % Data structure for direct commands and other configuration information
    % ----------------------------------------------------------------------
    Cfg=struct( 'ba', nan(nEnsembles,1),...
        'bc', nan(nEnsembles,1),...
        'be_mmps', nan(nEnsembles,1),...
        'bg', nan(nEnsembles,1),...
        'bm', nan(nEnsembles,1),...
        'bp', nan(nEnsembles,1),...
        'bx_dm', nan(nEnsembles,1),...
        'codeReps', nan(nEnsembles,1),...
        'coordSys', repmat(blanks(5),nEnsembles,1),...
        'cpuSerNo', nan(nEnsembles,8),...
        'cq', nan(nEnsembles,1),...
        'cx', nan(nEnsembles,1),...
        'distBin1_cm', nan(nEnsembles,1),...
        'ea_deg', nan(nEnsembles,1),...
        'eb_deg', nan(nEnsembles,1),...
        'ec', repmat(blanks(8),nEnsembles,1),...
        'ex', repmat(blanks(8),nEnsembles,1),...
        'ez', repmat(blanks(8),nEnsembles,1),...
        'headSrc', repmat(blanks(11),nEnsembles,1),...
        'lag_cm', nan(nEnsembles,1),...
        'mapBins', repmat(blanks(3),nEnsembles,1),...
        'nBeams', nan(nEnsembles,1),...
        'pitchSrc', repmat(blanks(11),nEnsembles,1),...
        'refLayEndCell', nan(nEnsembles,1),...
        'refLayStrCell', nan(nEnsembles,1),...
        'rollSrc', repmat(blanks(11),nEnsembles,1),...
        'salSrc', repmat(blanks(9),nEnsembles,1),...
        'wm', nan(nEnsembles,1),...
        'sosSrc', repmat(blanks(11),nEnsembles,1),...
        'tempSrc', repmat(blanks(11),nEnsembles,1),...
        'tp_sec', nan(nEnsembles,1),...
        'use3beam', repmat(blanks(3),nEnsembles,1),...
        'usePR', repmat(blanks(3),nEnsembles,1),...
        'wa', nan(nEnsembles,1),...
        'wb', nan(nEnsembles,1),...
        'wc', nan(nEnsembles,1),...
        'we_mmps', nan(nEnsembles,1),...
        'wf_cm', nan(nEnsembles,1),...
        'wg_per', nan(nEnsembles,1),...
        'wj', nan(nEnsembles,1),...
        'wn', nan(nEnsembles,1),...
        'wp', nan(nEnsembles,1),...
        'ws_cm', nan(nEnsembles,1),...
        'xdcrDepSrs', repmat(blanks(9),nEnsembles,1),...
        'xmitPulse_cm', nan(nEnsembles,1));
    %
    % Data structure for data obtained from the various internal and external
    % sensors.
    % -------------------------------------------------------------------------
    Sensor=struct(  'ambientTemp', nan(nEnsembles,1),...
        'attitudeTemp', nan(nEnsembles,1),...
        'attitude', nan(nEnsembles,1),...
        'bitTest', nan(nEnsembles,1),...
        'contamSensor', nan(nEnsembles,1),...
        'date', nan(nEnsembles,4),...
        'dateNotY2k', nan(nEnsembles,3),...
        'dateY2k', nan(nEnsembles,4),...
        'errorStatusWord', repmat(blanks(8),[nEnsembles,1,4]),...
        'headingStdDev', nan(nEnsembles,1),...
        'heading_deg', nan(nEnsembles,1),...
        'mpt_msc', nan(nEnsembles,3),...
        'num', nan(nEnsembles,1),...
        'numFact', nan(nEnsembles,1),...
        'orient', repmat(blanks(4),nEnsembles,1),...
        'pitchStdDev', nan(nEnsembles,1),...
        'pitch_deg', nan(nEnsembles,1),...
        'pressureNeg', nan(nEnsembles,1),...
        'pressurePos', nan(nEnsembles,1),...
        'pressureVar_pascal', nan(nEnsembles,1)/100,... % modified by SAMoore 2013-09-30
        'pressure_pascal', nan(nEnsembles,1)/100,...  % modified by SAMoore 2013-09-30
        'rollStdDev_deg', nan(nEnsembles,1),...
        'roll_deg', nan(nEnsembles,1),...
        'salinity_ppt', nan(nEnsembles,1),...
        'sos_mps', nan(nEnsembles,1),...
        'temperature_degc', nan(nEnsembles,1),...
        'time', nan(nEnsembles,4),...
        'timeY2k', nan(nEnsembles,4),...
        'xdcrDepth_dm', nan(nEnsembles,1),...
        'xmitCurrent', nan(nEnsembles,1),...
        'xmitVoltage', nan(nEnsembles,1));
    %
    % Data structure for the water track data. Data are stored in 3-dimensional
    % arrays with the 1st dimension (rows) being the bin number, the second
    % dimension (column) being the ensemble index, and the 3rd dimension being
    % the beam number.
    % -------------------------------------------------------------------------
    Wt=struct(  'corr', nan(nBins,nEnsembles,nBeams),...
        'pergd', nan(nBins,nEnsembles,nBeams),...
        'rssi', nan(nBins,nEnsembles,nBeams),...
        'vel_mps', nan(nBins,nEnsembles,nBeams));
    %
    % Data structure for the bottom track data. Data are stored in 2-dimensional
    % arrays with the 1st dimension (rows) being the beam number and the 2nd
    % dimension being the ensemble index.
    % -------------------------------------------------------------------------
    Bt=struct(  'corr', nan(nBeams,nEnsembles),...
        'depth_m', nan(nBeams,nEnsembles),...
        'evalAmp', nan(nBeams,nEnsembles),...
        'extDepth_cm', nan(nEnsembles,1),...
        'pergd', nan(nBeams,nEnsembles),...
        'rssi', nan(nBeams,nEnsembles),...
        'vel_mps', nan(nBeams,nEnsembles));
    %
    % Data structure for WinRiver 10.06 and previous GPS data
    % -------------------------------------------------------
    Gps=struct( 'alt_m', nan(nEnsembles,1),...
        'ggaDiff', nan(nEnsembles,1),...
        'ggaHdop', nan(nEnsembles,1),...
        'ggaNStats', nan(nEnsembles,1),...
        'ggaVelE_mps', nan(nEnsembles,1),...
        'ggaVelN_mps', nan(nEnsembles,1),...
        'gsaPdop', nan(nEnsembles,1),...
        'gsaSat', nan(nEnsembles,6),...
        'gsaVdop', nan(nEnsembles,1),...
        'lat_deg', nan(nEnsembles,1),...
        'long_deg', nan(nEnsembles,1),...
        'vtgVelE_mps', nan(nEnsembles,1),...
        'vtgVelN_mps', nan(nEnsembles,1));
    %
    % Data structure for WinRiver 2 GPS data
    % --------------------------------------
    % SAM 2013/07/31
    % I have no idea why there are supposed to be 20 elements
    %===============================
    Gps2=struct( 'ggaDeltaTime', nan(nEnsembles,20),...
        'ggaHeader', repmat(blanks(1),[nEnsembles 20 7]),...
        'utc', nan(nEnsembles,20),...
        'lat_dm', nan(nEnsembles,20),...
        'latRef', repmat(blanks(1),nEnsembles,20),...
        'lon_dm', nan(nEnsembles,20),...
        'lonRef', repmat(blanks(1),nEnsembles,20),...
        'corrQual', nan(nEnsembles,20),...
        'numSats', nan(nEnsembles,20),...
        'hdop', nan(nEnsembles,20),...
        'alt', nan(nEnsembles,20),...
        'altUnit', repmat(blanks(1),nEnsembles,20),...
        'geoid', nan(nEnsembles,20),...
        'geoidUnit', nan(nEnsembles,20),...
        'dgpsAge', nan(nEnsembles,20),...
        'refStatID', nan(nEnsembles,20),...
        'vtgDeltaTime', nan(nEnsembles,20),...
        'vtgHeader', repmat(blanks(1),[nEnsembles 20 7]),...
        'courseTrue', nan(nEnsembles,20),...
        'trueIndicator', repmat(blanks(1),nEnsembles,20),...
        'courseMag', nan(nEnsembles,20),...
        'magIndicator', repmat(blanks(1),nEnsembles,20),...
        'speedKnots', nan(nEnsembles,20),...
        'knotsIndicator', repmat(blanks(1),nEnsembles,20),...
        'speedKmph', nan(nEnsembles,20),...
        'kmphIndicator', repmat(blanks(1),nEnsembles,20),...
        'modeIndicator', repmat(blanks(1),nEnsembles,20), ...
        'ggaVelE_mps', nan(nEnsembles,1),...
        'ggaVelN_mps', nan(nEnsembles,1),...
        'vtgVelE_mps', nan(nEnsembles,1),...
        'vtgVelN_mps', nan(nEnsembles,1));
    %
    % Data structure for RiverRay Surface Cells
    % -----------------------------------------
    Surface=struct( 'no_cells', nan(nEnsembles,1),...
        'cell_size_cm', nan(nEnsembles,1),...
        'dist_bin1_cm', nan(nEnsembles,1),...
        'vel_mps', nan(5,nEnsembles,4),...
        'corr', nan(5,nEnsembles,4),...
        'rssi', nan(5,nEnsembles,4));
    %
    % Data structure for RiverRay autoconfiguration data
    % --------------------------------------------------
    AutoMode=struct('beam_count', nan(nEnsembles,1),...
        'Beam1',struct( 'mode', nan(nEnsembles,1),...
        'depth_cm', nan(nEnsembles,1),...
        'ping_count', nan(nEnsembles,1),...
        'ping_type', nan(nEnsembles,1),...
        'cell_count', nan(nEnsembles,1),...
        'cell_size_cm', nan(nEnsembles,1),...
        'cell_mid_cm', nan(nEnsembles,1),...
        'code_repeat', nan(nEnsembles,1),...
        'trans_length_cm', nan(nEnsembles,1),...
        'lag_length_cm', nan(nEnsembles,1),...
        'transmit_bw', nan(nEnsembles,1),...
        'receive_bw', nan(nEnsembles,1),...
        'ping_interval_ms', nan(nEnsembles,1)),...
        'Beam2',struct( 'mode', nan(nEnsembles,1),...
        'depth_cm', nan(nEnsembles,1),...
        'ping_count', nan(nEnsembles,1),...
        'ping_type', nan(nEnsembles,1),...
        'cell_count', nan(nEnsembles,1),...
        'cell_size_cm', nan(nEnsembles,1),...
        'cell_mid_cm', nan(nEnsembles,1),...
        'code_repeat', nan(nEnsembles,1),...
        'trans_length_cm', nan(nEnsembles,1),...
        'lag_length_cm', nan(nEnsembles,1),...
        'transmit_bw', nan(nEnsembles,1),...
        'receive_bw', nan(nEnsembles,1),...
        'ping_interval_ms', nan(nEnsembles,1)),...
        'Beam3',struct( 'mode', nan(nEnsembles,1),...
        'depth_cm', nan(nEnsembles,1),...
        'ping_count', nan(nEnsembles,1),...
        'ping_type', nan(nEnsembles,1),...
        'cell_count', nan(nEnsembles,1),...
        'cell_size_cm', nan(nEnsembles,1),...
        'cell_mid_cm', nan(nEnsembles,1),...
        'code_repeat', nan(nEnsembles,1),...
        'trans_length_cm', nan(nEnsembles,1),...
        'lag_length_cm', nan(nEnsembles,1),...
        'transmit_bw', nan(nEnsembles,1),...
        'receive_bw', nan(nEnsembles,1),...
        'ping_interval_ms', nan(nEnsembles,1)),...
        'Beam4',struct( 'mode', nan(nEnsembles,1),...
        'depth_cm', nan(nEnsembles,1),...
        'ping_count', nan(nEnsembles,1),...
        'ping_type', nan(nEnsembles,1),...
        'cell_count', nan(nEnsembles,1),...
        'cell_size_cm', nan(nEnsembles,1),...
        'cell_mid_cm', nan(nEnsembles,1),...
        'code_repeat', nan(nEnsembles,1),...
        'trans_length_cm', nan(nEnsembles,1),...
        'lag_length_cm', nan(nEnsembles,1),...
        'transmit_bw', nan(nEnsembles,1),...
        'receive_bw', nan(nEnsembles,1),...
        'ping_interval_ms', nan(nEnsembles,1)),...
        'Reserved', nan(nEnsembles,1));
    
    
    % Data structure for raw NMEA data strings. The strings are stored
    % in a character array with rows being the ensemble index and columns
    % the characters in the respective string.
    % --------------------------------------------------------------------
    Nmea=struct('gga', repmat(blanks(97),nEnsembles,1),...
        'gsa', repmat(blanks(60),nEnsembles,1),...
        'vtg', repmat(blanks(50),nEnsembles,1),...
        'raw', repmat(blanks(100),nEnsembles,1));
    %
    % Read Raw Data
    % =============
    % Data in the PD0 format are organized by leader_id. The leader_id is read
    % and the appropriate statements executed to read the data defined by the
    % leader_id. All data are read and stored, even data that should not change
    % between ensembles. The data are read until the end of file character is
    % encountered.
    % -------------------------------------------------------------------------
    %
    % Reset file pointer
    % ------------------
    waitbar(0.2);
    fseek(fid,initialPos,'bof');
    %
    % Initialize variables
    % --------------------
    iEns=0;
    disp('Reading file')
    endFileCheck=0;
    endFile=FileInfo.bytes;
    idataTypes=0;
    nDataTypes=1;
    %
    % Begin processing the file
    % =========================
    while endFileCheck<endFile
        %
        % Read leader_id
        % --------------
        leader_id=dec2hex(fread(fid,1,'uint16'),4);
        %
        % Fix for situation where last ensemble is invalid or incomplete
        % --------------------------------------------------------------
        if idataTypes>=nDataTypes && ~strcmp(leader_id,'7F7F')
            leader_id='9999'
        end
        
        %
        % Select appropriate code to read data based on leader_id
        % =======================================================
        switch leader_id
            %
            % Read Binary Header Data
            % =======================
            case '7F7F'
                i2022=0;
                fileloc=ftell(fid)-2;
                %
                % Check for end of file
                % ---------------------
                if fileloc+bytesPerEns > endFile
                    endFileCheck=endFile+1;
                else
                    %
                    % Check and adjust for lost ensembles
                    % ===================================
                    idataTypes=0;
                    storefileloc=ftell(fid);
                    bytesPerEns=fread(fid,1,'uint16');
                    %
                    % Check checksum for valid ensemble
                    % ---------------------------------
                    fseek(fid,fileloc,'bof');
                    Testb = fread(fid,bytesPerEns,'uchar');
                    check=sum(Testb);
                    checkh=dec2hex(check);
                    checkh=checkh(end-3:end);
                    fseek(fid,fileloc+bytesPerEns,'bof');
                    checksum=fread(fid,1,'uint16');
                    if hex2dec(checkh)==checksum
                        %
                        % If checksum is valid read header data
                        % -------------------------------------
                        dataOK=1;
                        fseek(fid,fileloc+5,'bof');
                        nDataTypes=fread(fid,1,'uint8');
                        dataOffsets(1:nDataTypes)=fread(fid,nDataTypes,'uint16');
                        %
                        % Find variable leader id
                        % -----------------------
                        while (idataTypes+1)<=nDataTypes && ~strcmp(leader_id,'0080')
                            fseek(fid,(dataOffsets(idataTypes+1)+fileloc),'bof');
                            leader_id=dec2hex(fread(fid,1,'uint16'),4);
                            idataTypes=idataTypes+1;
                        end
                        %
                        % Check for consecutive ensembles numbers
                        % ---------------------------------------
                        if iEns>0 && strcmp(leader_id,'0080')
                            ensNum=fread(fid,1,'uint16');
                            ensNumdiff=ensNum-Sensor.num(iEns);
                            %
                            % Insert lost ensemble numbers and increment iEns
                            % -----------------------------------------------
                            if ensNumdiff>1
                                for nn=1:ensNumdiff-1
                                    Sensor.num(iEns+1)=Sensor.num(iEns)+1;
                                    iEns=iEns+1;
                                end
                            end
                        end
                    else
                        %
                        % If checksum is invalid look for next ensemble
                        % ---------------------------------------------
                        disp('Bad Checksum New Code');
                        search_id='    ';
                        searchloc=fileloc+2;
                        while search_id~=hex2dec('7F7F');
                            searchloc=searchloc+1;
                            fseek(fid,searchloc,'bof');
                            search_id=fread(fid,1,'uint16');
                        end
                        fseek(fid,searchloc,'bof');
                        idataTypes=-1;
                    end
                    %
                    % Reset file location
                    % -------------------
                    fseek(fid,storefileloc,'bof');
                    %
                    % Initialize variables
                    % --------------------
                    idataTypes=0;
                    j100=0;
                    j101=0;
                    j102=0;
                    j103=0;
                    iEns=iEns+1;
                    waitbar(0.2+0.7.*iEns/nEnsembles);
                    %
                    % Read bytes in this ensemble
                    % ---------------------------
                    Hdr.bytesPerEns(iEns)=fread(fid,1,'uint16');
                    bytesPerEns=Hdr.bytesPerEns(iEns);
                    %
                    % Check checksum for valid ensemble
                    % ---------------------------------
                    fseek(fid,fileloc,'bof');
                    Testb = fread(fid,Hdr.bytesPerEns(iEns),'uchar');
                    check=sum(Testb);
                    checkh=dec2hex(check);
                    checkh=checkh(end-3:end);
                    fseek(fid,fileloc+Hdr.bytesPerEns(iEns),'bof');
                    checksum=fread(fid,1,'uint16');
                    if hex2dec(checkh)==checksum
                        %
                        % If checksum is valid read header data
                        % -------------------------------------
                        Hdr.dataOK(iEns)=1;
                        fseek(fid,fileloc+5,'bof');
                        Hdr.nDataTypes(iEns)=fread(fid,1,'uint8');
                        Hdr.dataOffsets(iEns,1:Hdr.nDataTypes(iEns))=fread(fid,Hdr.nDataTypes(iEns),'uint16');
                        if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                            fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                        else
                            fseek(fid,fileloc+bytesPerEns-2);
                        end
                    else
                        %
                        % If checksum is invalid look for next ensemble
                        % ---------------------------------------------
                        disp('Bad Checksum');
                        search_id='    ';
                        searchloc=fileloc+2;
                        while search_id~=hex2dec('7F7F');
                            searchloc=searchloc+1;
                            fseek(fid,searchloc,'bof');
                            search_id=fread(fid,1,'uint16');
                        end
                        fseek(fid,searchloc,'bof');
                        idataTypes=-1;
                    end
                end
                %
                % Read Binary Fixed Leader Data
                % =============================
            case '0000'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read and decode firmware version
                % --------------------------------
                Inst.firmVer(iEns)=fread(fid,1,'uint8');
                Inst.firmVer(iEns)=Inst.firmVer(iEns)+fread(fid,1,'uint8')/100;
                %
                % Read and decode instrument characteristics
                % ------------------------------------------
                bitls=fread(fid,1,'uint8');
                bitls=dec2base(bitls,2,8);
                
                switch base2dec(bitls(6:8),2)
                    case 0;     Inst.freq(iEns)=75;
                    case 1;     Inst.freq(iEns)=150;
                    case 2;     Inst.freq(iEns)=300;
                    case 3;     Inst.freq(iEns)=600;
                    case 4;     Inst.freq(iEns)=1200;
                    case 5;     Inst.freq(iEns)=2400;
                    otherwise;  Inst.freq(iEns)=nan;
                end;
                
                switch base2dec(bitls(5),2)
                    case 0;     Inst.pat(iEns,:)='Concave';
                    case 1;     Inst.pat(iEns,:)='Convex ';
                    otherwise;   Inst.pat(iEns,:)='n/a    ';
                end;
                
                Inst.sensorCfg(iEns)=base2dec(bitls(3:4),2)+1;
                
                switch base2dec(bitls(2),2)
                    case 0;     Inst.xducer(iEns,:)='Not Attached';
                    case 1;     Inst.xducer(iEns,:)='Attached    ';
                    otherwise;  Inst.xducer(iEns,:)='n/a         ';
                end;
                
                switch base2dec(bitls(1),2)
                    case 0;     Sensor.orient(iEns,:)='Down';
                    case 1;     Sensor.orient(iEns,:)='Up  ';
                    otherwise;  Sensor.orient(iEns,:)='n/a ';
                end;
                
                bitms=fread(fid,1,'uint8');
                bitms=dec2base(bitms,2,8);
                
                switch base2dec(bitms(7:8),2)
                    case 0;     Inst.beamAng(iEns)=15;
                    case 1;     Inst.beamAng(iEns)=20;
                    case 2;     Inst.beamAng(iEns)=30;
                    case 3;     Inst.beamAng(iEns)=nan;
                    otherwise;  Inst.beamAng(iEns)=nan;
                end;
                
                switch base2dec(bitms(1:4),2)
                    case 4
                        Inst.beams(iEns)=4;
                    case 5
                        Inst.beams(iEns)=5;
                        Inst.demod(iEns)=1;
                    case 15
                        Inst.beams(iEns)=5;
                        Inst.demod(iEns)=2;
                    otherwise
                        Inst.beams(iEns)=nan;
                        Inst.demod(iEns)=nan;
                end;
                
                switch fread(fid,1,'uint8')
                    case 0;     Inst.dataType(iEns,:)='Real';
                    otherwise;  Inst.dataType(iEns,:)='Simu';
                end;
                %
                % Position file pointer and read configuration information
                % --------------------------------------------------------
                fseek(fid,1,'cof');
                Cfg.nBeams(iEns)=fread(fid,1,'uint8');
                Cfg.wn(iEns)=fread(fid,1,'uint8');
                Cfg.wp(iEns)=fread(fid,1,'uint16');
                Cfg.ws_cm(iEns)=fread(fid,1,'uint16');
                Cfg.wf_cm(iEns)=fread(fid,1,'uint16');
                Cfg.wm(iEns)=fread(fid,1,'uint8');
                Cfg.wc(iEns)=fread(fid,1,'uint8');
                Cfg.codeReps(iEns)=fread(fid,1,'uint8');
                Cfg.wg_per(iEns)=fread(fid,1,'uint8');
                Cfg.we_mmps(iEns)=fread(fid,1,'uint16');
                Cfg.tp_sec(iEns,:)=sum(fread(fid,3,'uint8').*[60 1 0.01]');
                Cfg.ex(iEns,:)=dec2base(fread(fid,1,'uint8'),2,8);
                
                switch base2dec(Cfg.ex(iEns,4:5),2)
                    case 0;     Cfg.coordSys(iEns,:)='Beam ';
                    case 1;     Cfg.coordSys(iEns,:)='Inst ';
                    case 2;     Cfg.coordSys(iEns,:)='Ship ';
                    case 3;     Cfg.coordSys(iEns,:)='Earth';
                    otherwise;  Cfg.coordSys(iEns,:)='n/a  ';
                end;
                
                switch base2dec(Cfg.ex(iEns,6),2)
                    case 0;     Cfg.usePR(iEns,:)='No ';
                    case 1;     Cfg.usePR(iEns,:)='Yes';
                    otherwise;  Cfg.usePR(iEns,:)='n/a';
                end;
                
                switch base2dec(Cfg.ex(iEns,7),2)
                    case 0;     Cfg.use3beam(iEns,:)='No ';
                    case 1;     Cfg.use3beam(iEns,:)='Yes';
                    otherwise;  Cfg.use3beam(iEns,:)='n/a';
                end;
                
                switch base2dec(Cfg.ex(iEns,8),2)
                    case 0;     Cfg.mapBins(iEns,:)='No ';
                    case 1;     Cfg.mapBins(iEns,:)='Yes';
                    otherwise;  Cfg.mapBins(iEns,:)='n/a';
                end;
                
                Cfg.ea_deg(iEns)=fread(fid,1,'int16')*0.01;
                Cfg.eb_deg(iEns)=fread(fid,1,'uint16')*0.01;
                Cfg.ez(iEns,:)=dec2base(fread(fid,1,'uint8'),2,8);
                
                switch base2dec(Cfg.ez(iEns,1:2),2)
                    case 0;     Cfg.sosSrc(iEns,:)='Manual EC  ';
                    case 1;     Cfg.sosSrc(iEns,:)='Calculated ';
                    case 3;     Cfg.sosSrc(iEns,:)='SVSS Sensor';
                    otherwise;  Cfg.sosSrc(iEns,:)='n/a        ';
                end
                
                switch base2dec(Cfg.ez(iEns,3),2)
                    case 0;     Cfg.xdcrDepSrc(iEns,:)='Manual ED';
                    case 1;     Cfg.xdcrDepSrc(iEns,:)='Sensor   ';
                    otherwise;  Cfg.xdcrDepSrc(iEns,:)='n/a      ';
                end
                
                switch base2dec(Cfg.ez(iEns,4),2)
                    case 0;     Cfg.headSrc(iEns,:)='Manual EH  ';
                    case 1;     Cfg.headSrc(iEns,:)='Int. Sensor';
                    otherwise;  Cfg.headSrc(iEns,:)='n/a        ';
                end
                
                switch base2dec(Cfg.ez(iEns,5),2)
                    case 0;     Cfg.pitchSrc(iEns,:)='Manual EP  ';
                    case 1;     Cfg.pitchSrc(iEns,:)='Int. Sensor';
                    otherwise;  Cfg.pitchSrc(iEns,:)='n/a        ';
                end
                
                switch base2dec(Cfg.ez(iEns,6),2)
                    case 0;     Cfg.rollSrc(iEns,:)='Manual ER  ';
                    case 1;     Cfg.rollSrc(iEns,:)='Int. Sensor';
                    otherwise;  Cfg.rollSrc(iEns,:)='n/a        ';
                end
                
                switch base2dec(Cfg.ez(iEns,7),2)
                    case 0;     Cfg.salSrc(iEns,:)='Manual ES';
                    case 1;     Cfg.salSrc(iEns,:)='n/a      ';
                    otherwise;  Cfg.salSrc(iEns,:)='n/a      ';
                end
                
                switch base2dec(Cfg.ez(iEns,8),2)
                    case 0;     Cfg.tempSrc(iEns,:)='Manual ET  ';
                    case 1;     Cfg.tempSrc(iEns,:)='Int. Sensor';
                    otherwise;  Cfg.tempSrc(iEns,:)='n/a        ';
                end
                
                Cfg.sensorAvail(iEns,:)=dec2base(fread(fid,1,'uint8'),2,8); % changed by SAMoore 2013/07/16
                Cfg.distBin1_cm(iEns)=fread(fid,1,'uint16');
                Cfg.xmitPulse_cm(iEns)=fread(fid,1,'uint16');
                Cfg.refLayStrCell(iEns)=fread(fid,1,'uint8'); % if zero, than there is not water reference level averaging
                Cfg.refLayEndCell(iEns)=fread(fid,1,'uint8');
                Cfg.wa(iEns)=fread(fid,1,'uint8');
                Cfg.cx(iEns)=fread(fid,1,'uint8');
                Cfg.lag_cm(iEns)=fread(fid,1,'uint16');
                Cfg.cpuSerNo(iEns,:)=fread(fid,8,'uint8');
                Cfg.wb(iEns)=fread(fid,1,'uint8');
                Cfg.cq(iEns)=fread(fid,1,'uint8');
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Variable Leader Data
                % =========================
            case '0080'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read instrument clock and sensor data
                % -------------------------------------
                Sensor.num(iEns)=fread(fid,1,'uint16');
                Sensor.dateNotY2k(iEns,:)=fread(fid,3,'uint8');
                Sensor.time(iEns,:)=fread(fid,4,'uint8');
                Sensor.numFact(iEns)=fread(fid,1,'uint8');
                Sensor.numTot(iEns)=Sensor.num(iEns)+Sensor.numFact(iEns)*65535;
                Sensor.bitTest(iEns)=fread(fid,1,'uint16');
                Sensor.sos_mps(iEns)=fread(fid,1,'uint16');
                Sensor.xdcrDepth_dm(iEns)=fread(fid,1,'uint16');
                Sensor.heading_deg(iEns)=fread(fid,1,'uint16')/100;
                Sensor.pitch_deg(iEns)=fread(fid,1,'int16')/100;
                Sensor.roll_deg(iEns)=fread(fid,1,'int16')/100;
                Sensor.salinity_ppt(iEns)=fread(fid,1,'uint16');
                Sensor.temperature_degc(iEns)=fread(fid,1,'uint16')./100;
                Sensor.mpt_msc(iEns,:)=fread(fid,3,'uint8');
                Sensor.headingStdDev_deg(iEns)=fread(fid,1,'uint8');
                Sensor.pitchStdDev_deg(iEns)=fread(fid,1,'uint8')/10;
                Sensor.rollStdDev_deg(iEns)=fread(fid,1,'uint8')/10;
                Sensor.xmitCurrent(iEns)=fread(fid,1,'uint8');
                Sensor.xmitVoltage(iEns)=fread(fid,1,'uint8');
                Sensor.ambientTemp(iEns)=fread(fid,1,'uint8');
                Sensor.pressurePos(iEns)=fread(fid,1,'uint8');
                Sensor.pressureNeg(iEns)=fread(fid,1,'uint8');
                Sensor.attitudeTemp(iEns)=fread(fid,1,'uint8');
                Sensor.atttitude(iEns)=fread(fid,1,'uint8');
                Sensor.contamSensor(iEns)=fread(fid,1,'uint8');
                Sensor.errorStatusWord(iEns,:,1)=dec2base(fread(fid,1,'uint8'),2,8);
                Sensor.errorStatusWord(iEns,:,2)=dec2base(fread(fid,1,'uint8'),2,8);
                Sensor.errorStatusWord(iEns,:,3)=dec2base(fread(fid,1,'uint8'),2,8);
                Sensor.errorStatusWord(iEns,:,4)=dec2base(fread(fid,1,'uint8'),2,8);
                fseek(fid,2,'cof');
                % SAMoore 2013-09-30 % pressure is given in deca pascals
                % see p. 137 of the Nov 2007 edition of the workhorse
                % commands and output data format
                Sensor.pressure_pascal(iEns)=fread(fid,1,'uint32')/100; % 
                Sensor.pressureVar_pascal(iEns)=fread(fid,1,'uint32')/100; %
                
                fseek(fid,1,'cof');
                Sensor.dateY2k(iEns,:)=fread(fid,4,'uint8');
                Sensor.timeY2k(iEns,:)=fread(fid,4,'uint8');
                Sensor.date(iEns,:)=Sensor.dateNotY2k(iEns);
                Sensor.date(iEns,1)=Sensor.dateY2k(iEns,1)*100+Sensor.dateY2k(iEns,2);
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Water-Tracking Velocity Data
                % =================================
            case '0100'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read water-tracking velocity data
                % ---------------------------------
                for iBin=1:Cfg.wn(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'int16');
                        if dummy~=-32768
                            Wt.vel_mps(iBin,iEns,iBeam)=dummy/1000;
                        else
                            Wt.vel_mps(iBin,iEns,iBeam)=nan;
                        end;
                    end;
                end;
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Correlation Magnitude
                % --------------------------
            case '0200'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read water-tracking correlation data
                % ------------------------------------
                for iBin=1:Cfg.wn(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Wt.corr(iBin,iEns,iBeam)=dummy;
                        else
                            Wt.corr(iBin,iEns,iBeam)=nan;
                        end;
                    end;
                end;
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Echo Intensity
                % ===================
            case '0300'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read water-tracking RSSI
                % ------------------------
                for iBin=1:Cfg.wn(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Wt.rssi(iBin,iEns,iBeam)=dummy;
                        else
                            Wt.rssi(iBin,iEns,iBeam)=nan;
                        end;
                    end;
                end;
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Percent-Good Data
                % ======================
            case '0400'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read water-tracking percent good data
                % -------------------------------------
                for iBin=1:Cfg.wn(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Wt.pergd(iBin,iEns,iBeam)=dummy;
                        else
                            Wt.pergd(iBin,iEns,iBeam)=nan;
                        end;
                    end;
                end;
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Status
                % ===========
                %        case '5000'
                %       Data for this data type are not decoded.
                %       ----------------------------------------
                %
                %           Update data types counter
                %           -------------------------
                %            idataTypes=idataTypes+1;
                %            if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                %                fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                %            else
                %                fseek(fid,fileloc+bytesPerEns-2,'bof');
                %            end
                %
                % Read Bottom Track Data
                % ======================
            case '0600'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Read bottom track configuration data
                % ------------------------------------
                Cfg.bp(iEns)=fread(fid,1,'uint16');
                long1=fread(fid,1,'uint16');
                Cfg.bc(iEns)=fread(fid,1,'uint8');
                Cfg.ba(iEns)=fread(fid,1,'uint8');
                Cfg.bg(iEns)=fread(fid,1,'uint8');
                Cfg.bm(iEns)=fread(fid,1,'uint8');
                Cfg.be_mmps(iEns)=fread(fid,1,'uint16');
                %
                % Read WinRiver 10.06 format GPS data
                % -----------------------------------
                Gps.lat_deg(iEns,1)=(fread(fid,1,'int32')/2^31)*180; %
            
                % Read least significant byte for beam depths
                % -------------------------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'uint16');
                    if dummy~=-32768
                        Bt.depth_m(iBeam,iEns)=dummy/100;
                    else
                        Bt.depth_m(iBeam,iEns)=dummy/100;
                    end;
                end;
                %
                % Read bottom-track velocities, in coordinate system specified
                % in Cfg.coordsys
                % ----------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'int16');
                    if dummy~=-32768
                        Bt.vel_mps(iBeam,iEns)=dummy/1000;
                    else
                        Bt.vel_mps(iBeam,iEns)=nan;
                    end;
                end;
                %
                % Read bottom-track correlations
                % -------------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'uint8');
                    if dummy~=-32768
                        Bt.corr(iBeam,iEns)=dummy;
                    else
                        Bt.corr(iBeam,iEns)=nan;
                    end;
                end;
                %
                % Read bottom-track evaluation amplitude
                % --------------------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'uint8');
                    if dummy~=-32768
                        Bt.evalAmp(iBeam,iEns)=dummy;
                    else
                        Bt.evalAmp(iBeam,iEns)=nan;
                    end;
                end;
                %
                % Read bottom-track percent good
                % ------------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'uint8');
                    if dummy~=-32768
                        Bt.pergd(iBeam,iEns)=dummy;
                    else
                        Bt.pergd(iBeam,iEns)=nan;
                    end
                end;
                %
                % Read WinRiver 10.06 format GPS data
                % -----------------------------------
                dummy=fread(fid,1,'uint16');
                if dummy~=-32768
                    Gps.alt_m(iEns,1)=(dummy-32768)/10;
                else
                    Gps.alt_m(iEns,1)=nan;
                end
                long2=fread(fid,1,'uint16');
                Gps.long_deg(iEns,1)=((long1+long2*2^16)/2^31)*180;
                if Gps.long_deg(iEns,1) > 180
                    Gps.long_deg(iEns,1)=Gps.long_deg(iEns,1)-360;
                end;
                Bt.extDepth_cm(iEns)=fread(fid,1,'int16');
                dummy=fread(fid,1,'int16');
                if dummy~=-32768
                    Gps.ggaVelE_mps(iEns,1)=dummy*-1/1000;
                else
                    Gps.ggaVelE_mps(iEns,1)=nan;
                end
                dummy=fread(fid,1,'int16');
                if dummy~=-32768
                    Gps.ggaVelN_mps(iEns,1)=dummy*-1/1000;
                else
                    Gps.ggaVelN_mps(iEns,1)=nan;
                end
                dummy=fread(fid,1,'int16');
                if dummy~=-32768
                    Gps.vtgVelE_mps(iEns,1)=dummy*-1/1000;
                else
                    Gps.vtgVelE_mps(iEns,1)=nan;
                end
                dummy=fread(fid,1,'int16');
                if dummy~=-32768
                    Gps.vtgVelN_mps(iEns,1)=dummy*-1/1000;
                else
                    Gps.vtgVelN_mps(iEns,1)=nan;
                end
                dummy=fread(fid,1,'uint8');
                if dummy~=0
                    Gps.gsaVdop(iEns,1)=dummy;
                end
                dummy=fread(fid,1,'uint8');
                if dummy~=0
                    Gps.gsaPdop(iEns,1)=dummy;
                end
                dummy=fread(fid,1,'uint8');
                if dummy~=0
                    Gps.ggaNStats(iEns,1)=dummy;
                end
                fseek(fid,1,'cof');
                Gps.gsaSat(iEns,5)=fread(fid,1,'uint8');
                Gps.gsaSat(iEns,6)=fread(fid,1,'uint8');
                Gps.ggaDiff(iEns,1)=fread(fid,1,'uint8');
                dummy=fread(fid,1,'uint8');
                if dummy~=0
                    Gps.ggaHdop(iEns,1)=dummy/10;
                end
                Gps.gsaSat(iEns,1)=fread(fid,1,'uint8');
                Gps.gsaSat(iEns,2)=fread(fid,1,'uint8');
                Gps.gsaSat(iEns,3)=fread(fid,1,'uint8');
                Gps.gsaSat(iEns,4)=fread(fid,1,'uint8');
                %
                % Read bx configuration setting
                % -----------------------------
                Cfg.bx_dm(iEns)=fread(fid,1,'uint16');
                %
                % Read bottom-tracking RSSI
                % -------------------------
                Bt.rssi(1,iEns)=fread(fid,1,'uint8');
                Bt.rssi(2,iEns)=fread(fid,1,'uint8');
                Bt.rssi(3,iEns)=fread(fid,1,'uint8');
                Bt.rssi(4,iEns)=fread(fid,1,'uint8');
                %
                % Read wj configuration setting
                % -----------------------------
                Cfg.wj(iEns)=fread(fid,1,'uint8');
                %
                % Read most significant byte and compute beam depths
                % --------------------------------------------------
                for iBeam=1:Inst.beams(iEns)
                    dummy=fread(fid,1,'uint8');
                    if dummy~=-32768
                        Bt.depth_m(iBeam,iEns)=Bt.depth_m(iBeam,iEns)+...
                            (dummy*2^16)/100;
                    end
                end
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read  General NMEA Structure
                % ============================
                % Data type '2022' contains sub data types the identify specfic NMEA
                % 0183 data types that will be decoded. There may be multiple values
                % for a single ensemble.
                % -------------------------------------------------------------------
            case '2022'
                %
                % Update data types counter
                % -------------------------
                i2022=i2022+1;
                idataTypes=idataTypes+1;
                specificID=fread(fid,1,'int16');
                msgSize=fread(fid,1,'int16');
                deltaTime=fread(fid,1,'double'); %8,uchar
                switch specificID
                    %
                    % Read GGA data
                    % -------------
                    case 100
                        j100=j100+1;
                        Gps2.ggaDeltaTime(iEns,j100)=deltaTime; % Time diff between end of ensemble and GPS
                        Gps2.ggaHeader(iEns,j100,:)=fread(fid,10,'*char')';
                        Gps2.utc(iEns,j100,:)=str2double(fread(fid,10,'*char')');
                        Gps2.lat_dm(iEns,j100)=fread(fid,1,'float64');
                        Gps2.latRef(iEns,j100)=fread(fid,1,'*char');
                        Gps2.lon_dm(iEns,j100)=fread(fid,1,'float64');
                        Gps2.lonRef(iEns,j100)=fread(fid,1,'*char');
                        Gps2.corrQual(iEns,j100)=fread(fid,1,'uint8');
                        Gps2.numSats(iEns,j100)=fread(fid,1,'uint8');
                        Gps2.hdop(iEns,j100)=fread(fid,1,'float32');
                        Gps2.alt(iEns,j100)=fread(fid,1,'float32');
                        Gps2.altUnit(iEns,j100)=fread(fid,1,'*char');
                        Gps2.geoid(iEns,j100)=fread(fid,1,'float32');
                        Gps2.geoidUnit(iEns,j100)=fread(fid,1,'*char');
                        Gps2.dgpsAge(iEns,j100)=fread(fid,1,'float32');
                        Gps2.refStatID(iEns,j100)=fread(fid,1,'int16');
                        %
                        % Read VTG data
                        % -------------
                    case 101
                        j101=j101+1;
                        Gps2.vtgDeltaTime(iEns,j101)=deltaTime;
                        Gps2.vtgHeader(iEns,j101,:)=fread(fid,10,'*char')';
                        Gps2.courseTrue(iEns,j101)=fread(fid,1,'float32');
                        Gps2.trueIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.courseMag(iEns,j101)=fread(fid,1,'float32');
                        Gps2.magIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.speedKnots(iEns,j101)=fread(fid,1,'float32');
                        Gps2.knotsIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.speedKmph(iEns,j101)=fread(fid,1,'float32');
                        Gps2.kmphIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.modeIndicator(iEns,j101)=fread(fid,1,'*char');
                        %
                        % Read DBT data
                        % -------------
                    case 102
                        j102=j102+1;
                        Gps2.dbtDeltaTime(iEns,j102)=deltaTime;
                        Gps2.dbtHeader(iEns,j102,:)=fread(fid,10,'*char')';
                        Gps2.depth_ft(iEns,j102)=fread(fid,1,'float32');
                        Gps2.ftIndicator(iEns,j102)=fread(fid,1,'*char');
                        Gps2.depth_m(iEns,j102)=fread(fid,1,'float32');
                        Gps2.mIndicator(iEns,j102)=fread(fid,1,'*char');
                        Gps2.depth_fath(iEns,j102)=fread(fid,1,'float32');
                        Gps2.fathIndicator(iEns,j102)=fread(fid,1,'*char');
                        %
                        % Read HDT data
                        % -------------
                    case 103
                        j103=j103+1;
                        Gps2.hdtDeltaTime(iEns,j103,:)=deltaTime;
                        Gps2.hdtHeader(iEns,j103,:)=fread(fid,10,'*char')';
                        Gps2.heading_deg(iEns,j103)=fread(fid,1,'float32');
                        Gps2.hTrueIndicactor(iEns,j103)=fread(fid,1,'*char');
                        %
                        % Read GGA data
                        % -------------
                    case 104
                        j100=j100+1;
                        Gps2.ggaDeltaTime(iEns,j100)=deltaTime;
                        Gps2.ggaHeader(iEns,j100,1:7)=fread(fid,7,'*char')';
                        Gps2.utc(iEns,j100)=str2double(fread(fid,10,'*char')');
                        Gps2.lat_dm(iEns,j100)=fread(fid,1,'float64');
                        Gps2.latRef(iEns,j100)=fread(fid,1,'*char');
                        Gps2.lon_dm(iEns,j100)=fread(fid,1,'float64');
                        Gps2.lonRef(iEns,j100)=fread(fid,1,'*char');
                        Gps2.corrQual(iEns,j100)=fread(fid,1,'uint8');
                        Gps2.numSats(iEns,j100)=fread(fid,1,'uint8');
                        Gps2.hdop(iEns,j100)=fread(fid,1,'float32');
                        Gps2.alt(iEns,j100)=fread(fid,1,'float32');
                        Gps2.altUnit(iEns,j100)=fread(fid,1,'*char');
                        Gps2.geoid(iEns,j100)=fread(fid,1,'float32');
                        Gps2.geoidUnit(iEns,j100)=fread(fid,1,'*char');
                        Gps2.dgpsAge(iEns,j100)=fread(fid,1,'float32');
                        Gps2.refStatID(iEns,j100)=fread(fid,1,'int16');
                        %
                        % Read VTG data
                        % -------------
                    case 105
                        j101=j101+1;
                        Gps2.vtgDeltaTime(iEns,j101)=deltaTime;
                        Gps2.vtgHeader(iEns,j101,:)=fread(fid,7,'*char')';
                        Gps2.courseTrue(iEns,j101)=fread(fid,1,'float32');
                        Gps2.trueIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.courseMag(iEns,j101)=fread(fid,1,'float32');
                        Gps2.magIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.speedKnots(iEns,j101)=fread(fid,1,'float32');
                        Gps2.knotsIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.speedKmph(iEns,j101)=fread(fid,1,'float32');
                        Gps2.kmphIndicator(iEns,j101)=fread(fid,1,'*char');
                        Gps2.modeIndicator(iEns,j101)=fread(fid,1,'*char');
                        %
                        % Read DBT data
                        % -------------
                    case 106
                        j102=j102+1;
                        Gps2.dbtDeltaTime(iEns,j102)=deltaTime;
                        Gps2.dbtHeader(iEns,j102,:)=fread(fid,7,'*char')';
                        Gps2.depth_ft(iEns,j102)=fread(fid,1,'float32');
                        Gps2.ftIndicator(iEns,j102)=fread(fid,1,'*char');
                        Gps2.depth_m(iEns,j102)=fread(fid,1,'float32');
                        Gps2.mIndicator(iEns,j102)=fread(fid,1,'*char');
                        Gps2.depth_fath(iEns,j102)=fread(fid,1,'float32');
                        Gps2.fathIndicator(iEns,j102)=fread(fid,1,'*char');
                        %
                        % Read HDT data
                        % -------------
                    case 107
                        j103=j103+1;
                        Gps2.hdtDeltaTime(iEns,j103,:)=deltaTime;
                        Gps2.hdtHeader(iEns,j103,:)=fread(fid,7,'*char')';
                        Gps2.heading_deg(iEns,j103)=fread(fid,1,'float32');
                        Gps2.hTrueIndicactor(iEns,j103)=fread(fid,1,'*char');
                end
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read DBT NMEA String
                % ====================
            case '2100'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Reposition file pointer
                % -----------------------
                fseek(fid,(Hdr.dataOffsets(iEns,idataTypes)+fileloc+2),'bof');
                %
                % Read DBT sentence
                % -----------------
                Nmea.dbt(iEns,:)=fread(fid,38,'*char');
                dummy=zeros(1,38);
                endstr=find(Nmea.dbt(iEns,:)==char(13));
                dummy(1,1:endstr)=Nmea.dbt(iEns,1:endstr);
                Nmea.dbt(iEns,:)=dummy(1,:);
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read GGA NMEA String
                % ====================
            case '2101'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Reposition file pointer
                % -----------------------
                fseek(fid,(Hdr.dataOffsets(iEns,idataTypes)+fileloc+4),'bof');
                %
                % Determine number of characters to read
                % --------------------------------------
                if idataTypes < Hdr.nDataTypes(iEns)
                    num2Read=Hdr.dataOffsets(iEns,idataTypes+1)-Hdr.dataOffsets(iEns,idataTypes)-4;
                else
                    num2Read=bytesPerEns-Hdr.dataOffsets(iEns,idataTypes)-6;
                end
                %
                % Read GGA sentence
                % -----------------
                dummy=fread(fid,num2Read,'*char');
                dummy1=blanks(97);
                endstr=size(dummy,1);
                dummy1(1,1:endstr)=dummy(1:endstr);
                Nmea.gga(iEns,:)=dummy1(1,:);
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read VTG NMEA String
                % ====================
            case '2102'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Reposition file pointer
                % -----------------------
                fseek(fid,(Hdr.dataOffsets(iEns,idataTypes)+fileloc+4),'bof');
                %
                % Determine number of characters to read
                % --------------------------------------
                if idataTypes < Hdr.nDataTypes(iEns)
                    num2Read=Hdr.dataOffsets(iEns,idataTypes+1)-Hdr.dataOffsets(iEns,idataTypes)-4;
                else
                    num2Read=bytesPerEns-Hdr.dataOffsets(iEns,idataTypes)-6;
                end
                %
                % Read VTG sentence
                % -----------------
                dummy=fread(fid,num2Read,'*char');
                dummy1=blanks(50);
                endstr=size(dummy,1);
                dummy1(1,1:endstr)=dummy(1:endstr);
                Nmea.vtg(iEns,:)=dummy1(1,:);
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read GSA NMEA String
                % ====================
            case '2103'
                %
                % Update data types counter
                % -------------------------
                idataTypes=idataTypes+1;
                %
                % Reposition file pointer
                % -----------------------
                fseek(fid,(Hdr.dataOffsets(iEns,idataTypes)+fileloc+4),'bof');
                %
                % Determine number of characters to read
                % --------------------------------------
                if idataTypes < Hdr.nDataTypes(iEns)
                    num2Read=Hdr.dataOffsets(iEns,idataTypes+1)-Hdr.dataOffsets(iEns,idataTypes)-4;
                else
                    num2Read=bytesPerEns-Hdr.dataOffsets(iEns,idataTypes)-6;
                end
                %
                % Read GSA sentence
                % -----------------
                dummy=fread(fid,num2Read,'*char');
                dummy1=blanks(60);
                endstr=size(dummy,1);
                dummy1(1,1:endstr)=dummy(1:endstr);
                Nmea.gsa(iEns,:)=dummy1(1,:);
                %
                % Check if more data types need to be read and position file
                % pointer
                % -----------------------------------------------------------
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read RiverRay Surface Layer Leader
                % ==================================
            case '0010'
                idataTypes=idataTypes+1;
                Surface.no_cells(iEns)=fread(fid,1,'uint8');
                Surface.cell_size_cm(iEns)=fread(fid,1,'uint16');
                Surface.dist_bin1_cm(iEns)=fread(fid,1,'uint16');
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Surface Velocity Data
                % ==========================
            case '0110'
                idataTypes=idataTypes+1;
                for iBin=1:Surface.no_cells(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'int16');
                        if dummy~=-32768
                            Surface.vel_mps(iBin,iEns,iBeam)=dummy/1000;
                        end;
                    end;
                end;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Surface Correlation Magnitude
                % ==================================
            case '0210'
                idataTypes=idataTypes+1;
                for iBin=1:Surface.no_cells(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Surface.corr(iBin,iEns,iBeam)=dummy;
                        end;
                    end;
                end;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Surface Echo Intensity
                % ===========================
            case '0310'
                idataTypes=idataTypes+1;
                for iBin=1:Surface.no_cells(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Surface.rssi(iBin,iEns,iBeam)=dummy;
                        end;
                    end;
                end;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Surface Percent-Good Data
                % ==============================
            case '0410'
                idataTypes=idataTypes+1;
                for iBin=1:Surface.no_cells(iEns)
                    for iBeam=1:Inst.beams(iEns)
                        dummy=fread(fid,1,'uint8');
                        if dummy~=-32768
                            Surface.pergd(iBin,iEns,iBeam)=dummy;
                        end;
                    end;
                end;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % Read Status - Not implemented
                % =============================
            case '0510'
                idataTypes=idataTypes+1;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
                %
                % RiverRay autoconfiguration data.
                % ================================
            case '4401'
                idataTypes=idataTypes+1;
                AutoMode.beam_count(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.mode(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.depth_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.ping_count(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.ping_type(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.cell_count(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.cell_size_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.cell_mid_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.code_repeat(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.trans_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.lag_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam1.transmit_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.receive_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam1.ping_interval_ms(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.mode(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.depth_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.ping_count(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.ping_type(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.cell_count(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.cell_size_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.cell_mid_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.code_repeat(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.trans_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.lag_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam2.transmit_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.receive_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam2.ping_interval_ms(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.mode(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.depth_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.ping_count(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.ping_type(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.cell_count(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.cell_size_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.cell_mid_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.code_repeat(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.trans_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.lag_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam3.transmit_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.receive_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam3.ping_interval_ms(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.mode(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.depth_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.ping_count(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.ping_type(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.cell_count(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.cell_size_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.cell_mid_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.code_repeat(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.trans_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.lag_length_cm(iEns)=fread(fid,1,'uint16');
                AutoMode.Beam4.transmit_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.receive_bw(iEns)=fread(fid,1,'uint8');
                AutoMode.Beam4.ping_interval_ms(iEns)=fread(fid,1,'uint16');
                AutoMode.Reserved(iEns)=fread(fid,1,'uint8');
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    fseek(fid,fileloc+bytesPerEns-2,'bof');
                end
            otherwise
                %
                % If leader_id is not recognized, process next datatype or next
                % ensemble
                % -------------------------------------------------------------
                Hdr.invalid(iEns,:)=leader_id;
                idataTypes=idataTypes+1;
                if (idataTypes+1)<=Hdr.nDataTypes(iEns)
                    %
                    % Reposition file pointer for next data type
                    % ------------------------------------------
                    fseek(fid,(Hdr.dataOffsets(iEns,idataTypes+1)+fileloc),'bof');
                else
                    %
                    % Locate next ensemble if no more datatypes
                    % -----------------------------------------
                    if (idataTypes+1)>Hdr.nDataTypes(iEns)+1
                        currentloc=ftell(fid);
                        srchString=fread(fid,endFile-currentloc,'uchar');
                        srchString=dec2Hex(srchString);
                        srchString=srchString';
                        srchString=srchString(1:end);
                        nextEns=findstr(srchString,'7F7F');
                        if size(nextEns,2)>0
                            nextEns=(nextEns(1)-1)/2;
                            fseek(fid,currentloc+nextEns,'bof');
                            idataTypes=0;
                        else
                            endFileCheck=endFile+1;
                        end
                    else
                        fseek(fid,fileloc+bytesPerEns-2,'bof');
                    end
                end
        end
        %
        % If all datatypes have been read, read last two bytes of ensemble
        % ----------------------------------------------------------------
        if idataTypes>=Hdr.nDataTypes(iEns)
            Inst.resRDI=fread(fid,1,'uint16');
            checksum=fread(fid,1,'uint16');
        end
        %
        % Check end of file
        % -----------------
        if endFileCheck<endFile
            endFileCheck=ftell(fid);
        end
    end
    
    %
    % If vtg data are available compute north and east components
    % -----------------------------------------------------------
    if Gps2.vtgHeader(1,1,1)=='$'
        %
        % Find minimum of absolute value of delta time from raw data
        % ----------------------------------------------------------
        vtgDeltaTime=abs(Gps2.vtgDeltaTime);
        vtgMin=nanmin(vtgDeltaTime,[],2);
        %
        % Compute the velocity components in m/s
        % --------------------------------------
        for i=1:size(vtgDeltaTime,1)
            idx=find(vtgDeltaTime(i,:)==vtgMin(i),1,'first');
            [Gps2.vtgVelE_mps(i),Gps2.vtgVelN_mps(i)]=...
                pol2cart((90-Gps2.courseTrue(i,idx))*pi/180,...
                Gps2.speedKmph(i,idx).*0.2777778);
        end
    end
    %
    % If gga data are available compute north and east components
    % -----------------------------------------------------------
    if Gps2.ggaHeader(1,1,1)=='$'
        %
        % Initialize variables
        % --------------------
        clear idx
        eRadius=6378137;
        coeff=eRadius.*pi/180;
        ellip=1./298.257223563;
        %
        % Find minimum of absolute value of delta time from raw data
        % ----------------------------------------------------------
        ggaDeltaTime=abs(Gps2.ggaDeltaTime);
        ggaMin=nanmin(ggaDeltaTime,[],2);
        for i=1:size(ggaDeltaTime,1)
            idx(i)=find(ggaDeltaTime(i,:)==ggaMin(i),1,'first');
            if i>1
                %
                % Use TRDI's algorithm to compute velocity from gga positions
                % -----------------------------------------------------------
                L=(Gps2.lat_dm(i,idx(i))+Gps2.lat_dm(i-1,idx(i-1)))./2;
                sL=sind(L);
                RE=coeff.*(1+ellip.*sL.*sL);
                RN=coeff.*(1-2*ellip+3*ellip.*sL.*sL);
                dx=RE.*(Gps2.lon_dm(i,idx(i))-Gps2.lon_dm(i-1,idx(i-1))).*cosd(L);
                dy=RN.*(Gps2.lat_dm(i,idx(i))-Gps2.lat_dm(i-1,idx(i-1)));
                dt=Gps2.utc(i,idx(i))-Gps2.utc(i-1,idx(i-1));
                Gps2.ggaVelE_mps(i)=dx./dt;
                Gps2.ggaVelN_mps(i)=dy./dt;
            else
                Gps2.ggaVelE_mps(i)=nan;
                Gps2.ggaVelN_mps(i)=nan;
            end
        end
    end
    
    
    disp('File complete')
    waitbar(1);
    close(hwait);
else
    disp('Empty file')
    Hdr=nan;
    Inst=nan;
    Cfg=nan;
    Sensor=nan;
    Gps=nan;
    Wt=nan;
    Bt=nan;
    Nmea=nan;
    Gps2=nan;
end

% ************************************************************************
% **************************** FUNCTIONS *********************************
% ************************************************************************
%==========================================================================
function [numens]=numberofensembles(fid,FileInfo)
%
% This function finds the ensemble number of the first and last ensemble
% and computes the number of ensembles in the data file. This can no longer
% be done using the file size and number of bytes in the first ensemble
% because the number of bytes in each ensemble are variable and depend on
% the stored GPS data.
% =========================================================================
%
% Initialize variables
% --------------------
i=0;
leader_id='0000';
%
% Find first ensemble
% -------------------
while ~strcmp(leader_id,'7F7F') && i<FileInfo.bytes
    fseek(fid,i,'bof');
    i=i+1;
    leader_id=dec2hex(fread(fid,1,'uint16'),4);
end
%
% Call findeensno to report ensemble number of first ensemble
% -----------------------------------------------------------
firstnum=findensno(fid);
%
% Reinitialize variables
% ----------------------
i=0;
leader_id='0000';
lastnum=-1;
while lastnum<0
    %
    % Find last ensemble
    % ------------------
    while ~strcmp(leader_id,'7F7F') && i<FileInfo.bytes
        i=i+1;
        fseek(fid,-i,'eof');
        leader_id=dec2hex(fread(fid,1,'uint16'),4);
    end
    %
    % Call findeensno to report ensemble number of last ensemble
    % -----------------------------------------------------------
    lastnum=findensno(fid);
    leader_id='0000';
end
%
% Compute number of ensembles
% ---------------------------
numens=lastnum-firstnum+1;
%==========================================================================
function [ensno]=findensno(fid)
%
% This function assumes the current position of the file pointer is just
% after '7F7F'. The function then reads the ensemble header and
% works through the data offsets until the 00800 data type is found. The
% ensemble number is then read.
% =========================================================================
%
% Reposition file pointer to top of ensemble.
% -------------------------------------------
fileloc=ftell(fid)-2;
%
% Read number of bytes in the ensemble
% ------------------------------------
bytesPerEns=fread(fid,1,'uint16');
%
% Read the entire ensemble
% ------------------------
fseek(fid,fileloc,'bof');
Testb = fread(fid,bytesPerEns,'uchar');
%
% Compute checksum
% ----------------
check=sum(Testb);
checkh=dec2hex(check);
if length(checkh)>3
    checkh=checkh(end-3:end);
    %
    % Read checksum
    % -------------
    fseek(fid,fileloc+bytesPerEns,'bof');
    checksum=fread(fid,1,'uint16');
    %
    % Check validity of the checksum
    % ------------------------------
    if hex2dec(checkh)==checksum
        %
        % If checksum is valid read header information
        % --------------------------------------------
        fseek(fid,fileloc+5,'bof');
        nDataTypes=fread(fid,1,'uint8');
        dataOffsets=fread(fid,nDataTypes,'uint16');
        %
        % Initialize variables
        % --------------------
        i=0;
        leader_id='0000';
        %
        % Search for data type '0080'
        % ---------------------------
        while ~strcmp(leader_id,'0080') && i<=nDataTypes
            i=i+1;
            fseek(fid,(dataOffsets(i)+fileloc),'bof');
            leader_id=dec2hex(fread(fid,1,'uint16'),4);
        end
        %
        % Read enemble number from data type '0080'
        % -----------------------------------------
        if strcmp(leader_id,'0080')
            ensno=fread(fid,1,'uint16');
        else
            disp('Ensemble number not found')
        end
    else
        ensno=-1;
    end
else
    ensno=-1;
end