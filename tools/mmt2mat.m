function [MMT,MMT_Site_Info,MMT_Transects,MMT_Field_Config,...
          MMT_Active_Config,MMT_Summary_None,MMT_Summary_BT,...
          MMT_Summary_GGA,MMT_Summary_VTG,MMT_QAQC,MMT_MB_Transects,...
          MMT_MB_Field_Config,MMT_MB_Active_Config] = mmt2mat(infile)
            %
            % This function accepts and input file name of an MMT file created by TRDI
            % WinRiver II and parses the data from the MMT file's xml structure and
            % stores the data in Matlab structures.
            %
            % This code assumes that the first configuration file is the field
            % configuration. If the user has rearranged the order of the configuration
            % files the field configuration will not be correct.
            %
            % This code is compatible with WinRiver II 2.07 but allowances have been
            % made and it should be compatible with WinRiver II mmt files. Please
            % report any issues to dmueller@usgs.gov.
            %
            % David S. Mueller, Hydrologist, Office of Surface Water
            % Nathan Reynolds, Student, Illinois Water Science Center
            %
            % 10/25/2010
            % -----------------------------------------------------------------------
            % 09/24/2012 - JGL edit
            % Used new function, xml2struct, in order to improve processing speed.
            % Changed code to reflect the new structure returned by xml2struct.
            %
            % 6/6/2013 - DSM
            % Fixed areas in configuration section where YES/NO text fields
            % were not being decoded properly.
            %
            % 11/18/2013 - DSM
            % Fixed configuration DS to be compatible with changes made to
            % the mmt file in WinRiver II 2.12
            % ========================================================================
            %
            % Read mmt file into a generic Matlab data structure
            % --------------------------------------------------
            set(0, 'DefaulttextInterpreter', 'none');
            idx = find(infile=='\',1,'last');    

            % Display a waitbar
            fileName = infile(idx+1:end);
            waitmessage = ['Reading ' fileName];
            hwait = waitbar(0.05,waitmessage);
            waitbar(0.2);

            % Read file using xml2struct
            s = xml2struct1(infile);

            % Initialize data structures to null values
            % -----------------------------------------
            MMT =[];
            MMT_Active_Config =[];
            MMT_Field_Config =[];
            MMT_MB_Active_Config =[];
            MMT_MB_Field_Config =[];
            MMT_MB_Transects =[];
            MMT_QAQC =[];
            MMT_Site_Info =[];
            MMT_Summary_BT =[];
            MMT_Summary_GGA =[];
            MMT_Summary_None =[];
            MMT_Summary_VTG =[];
            MMT_Transects =[];

            waitbar(0.3);

            % MMT - Name, Version, Locked
            % ----------------------------
            MMT.Name = s.WinRiver.Project.Attributes.Name;
            MMT.Version = s.WinRiver.Project.Attributes.Version;

            if isstruct(s.WinRiver.Project.Locked)
                MMT.Locked = s.WinRiver.Project.Locked.Text;
            else
                MMT.Locked = nan;
            end

            waitbar(0.4);

            % MMT_Site_Information
            % --------------------
            % Remove 'Attributes' from structure because it is not used in
            % MMT_Site_Info
            s.WinRiver.Project.Site_Information =...
                rmfield(s.WinRiver.Project.Site_Information,'Attributes');

            % Get fieldnames of Site_Informationj
            fields_siteinfo = fieldnames(s.WinRiver.Project.Site_Information);

            % Get Site_Information data; 
            % all data is stored as text except for 'Water_Temperature' 
            s_end = 'Text';
            for i = 1:length(fields_siteinfo)   
                if ~isempty(s.WinRiver.Project.Site_Information.(fields_siteinfo{i}).(s_end))
                    if strcmp(fields_siteinfo{i},'Water_Temperature')
                        MMT_Site_Info.(fields_siteinfo{i}) =...
                            str2double(s.WinRiver.Project.Site_Information.(fields_siteinfo{i}).(s_end));
                    else
                        MMT_Site_Info.(fields_siteinfo{i}) =...
                            s.WinRiver.Project.Site_Information.(fields_siteinfo{i}).(s_end);
                    end % if water temp
                else
                    MMT_Site_Info.(fields_siteinfo{i}) = nan;
                end % if ~isempty
            end % for i

            waitbar(0.5);

            % Assign local variable to trans
            trans = s.WinRiver.Project.Site_Discharge.Transect;

            % If there is only 1 transect, convert it to a cell to process
            % properly.
            [trans] = convert2cell(trans);

            % Site_Discharge includes: Transect, Discharge_Summary
            % MMT_Transects
            for i = 1:length(trans)

                % MMT_Transects
                % -------------
                % .Checked
                MMT_Transects.Checked(i) =...
                    str2double(trans{i}.Attributes.Checked);

                % Files
                % -----
                % Assign local variable
                files = trans{i}.File;

                % If there is only 1 files, convert it to a cell to process
                % properly.
                [files] = convert2cell(files);

                for f = 1:length(files)

                    % .Path -> make a column cell; numFiles
                    MMT_Transects.Path{i,f} = files{f}.Attributes.PathName; 

                    % .Files -> make a column cell; numFiles
                    MMT_Transects.Files{i,f} = files{f}.Text;  

                    % .Number -> make a column cell; numFiles
                    MMT_Transects.Number{i,f} = files{f}.Attributes.TransectNmb;
                end % for f

                % Notes
                % -----
                if isfield(trans{i},'Note')

                    % Assign local variable
                    note = trans{i}.Note;

                    % If there is only 1 note, convert it to a cell to process
                    % properly.
                    [note] = convert2cell(note);
                    for n = 1:length(note)
                        MMT_Transects.NoteDate{i,n} = note{n}.Attributes.TimeStamp;
                        MMT_Transects.NoteText{i,n} = note{n}.Attributes.Text;
                    end % for n
                end % if Note

                % Configuration
                % --------------------
                % Assign local variable
                config = trans{i}.Configuration;

                % If there is only 1 configuration, convert it to a cell to process
                % properly.
                [config] = convert2cell(config);

                % Find the field and active config; 
                % attributes - > field = 0; active = 1
                act = nan(length(config),1);
                for j = 1:length(config)
                    a = str2double(config{j}.Attributes.Checked);
                    act(j) = a;
                end % for j

                % Find indices and pick the first occurance for field and active
                idx_field = find(act == 0,1,'first');
                idx_active = find(act == 1,1,'first');

                % If there is no field configuration, then treat the field
                % configuration as an active configuration.
                % This is what is done in the previous mmt2mat.m.
                if isempty(idx_field)
                    idx_field = 1;
                end % if isempty

                % MMT_Field_Config
                % ----------------
                field_config = config{idx_field};
                [MMT_Field_Config] =...
                    mmtconfig(field_config, i, MMT_Field_Config);

                % MMT_Active_Config
                % -----------------
                active_config = config{idx_active};
                [MMT_Active_Config] =...
                    mmtconfig(active_config, i, MMT_Active_Config);

            end

            waitbar(0.6);

            % Discharge Summary
            % -------------------
            if isfield(s.WinRiver.Project.Site_Discharge, 'Discharge_Summary')

                % Assign local variable
                dischargesummary =...
                    s.WinRiver.Project.Site_Discharge.Discharge_Summary;

                % Summary_None
                [MMT_Summary_None] =...
                    mmtqsum(dischargesummary.None,MMT_Summary_None);

                % Summary BT
                [MMT_Summary_BT] =...
                    mmtqsum(dischargesummary.BottomTrack,MMT_Summary_BT);

                % Summary GGA
                [MMT_Summary_GGA] =...
                    mmtqsum(dischargesummary.GGA,MMT_Summary_GGA);

                % Summary VTG
                [MMT_Summary_VTG] =...
                    mmtqsum(dischargesummary.VTG,MMT_Summary_VTG);
            end % if isfield

            waitbar(0.7);

            % Check to see if field exists
            if isfield(s.WinRiver.Project,'QA_QC')

                % QA_QC
                % -------------------
                % Assign local variable
                QAQC = s.WinRiver.Project.QA_QC;

                % Get fields
                fields_qaqc = fieldnames(QAQC);

                % Get data
                for q = 1:length(fields_qaqc)

                    % RG_Test
                    % -------------
                    % If field is RG_Test, get TestResult data
                    if strcmp(fields_qaqc{q},'RG_Test')
                        [QAQC.(fields_qaqc{q}).TestResult] =...
                            convert2cell(QAQC.(fields_qaqc{q}).TestResult);
                        for irg = 1:length(QAQC.(fields_qaqc{q}))
                            for tr = 1:length(QAQC.(fields_qaqc{q}).TestResult)
                                MMT_QAQC.RG_Test{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.Text.Text;
                                MMT_QAQC.RG_Test_TimeStamp{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.TimeStamp.Text;
                            end % for tr
                        end % for irq            
                    end % if strcmp

                    % Compass_Calibration
                    % -------------------
                    % If field is Compass_Calibration, get  data
                    if strcmp(fields_qaqc{q},'Compass_Calibration')
                        [QAQC.(fields_qaqc{q}).TestResult] =...
                            convert2cell(QAQC.(fields_qaqc{q}).TestResult);
                        for irg = 1:length(QAQC.(fields_qaqc{q}))
                            for tr = 1:length(QAQC.(fields_qaqc{q}).TestResult)
                                MMT_QAQC.Compass_Cal_Test{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.Text.Text;
                                MMT_QAQC.Compass_Cal_Timestamp{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.TimeStamp.Text;
                            end % for tr
                        end % for irq       
                    end % if strcmp

                    % Compass_Evaluation
                    % ------------------
                    % If field is Compass_Calibration, get  data
                    if strcmp(fields_qaqc{q},'Compass_Evaluation')
                        [QAQC.(fields_qaqc{q}).TestResult] =...
                            convert2cell(QAQC.(fields_qaqc{q}).TestResult);
                        for irg = 1:length(QAQC.(fields_qaqc{q}))
                            for tr = 1:length(QAQC.(fields_qaqc{q}).TestResult)
                                MMT_QAQC.Compass_Eval_Test{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.Text.Text;
                                MMT_QAQC.Compass_Eval_Timestamp{irg,tr} =...
                                    QAQC.(fields_qaqc{q}).TestResult{tr}.TimeStamp.Text;
                            end % for tr
                        end % for irq         
                    end % if strcmp

                    % Moving_Bed Transects
                    % --------------------
                    % If field is Moving_Bed_Test, get  data
                    if strcmp(fields_qaqc{q},'Moving_Bed_Test')
                        if isfield(QAQC.Moving_Bed_Test,'Transect')

                            % Assign local data variables
                            transects_mb = QAQC.Moving_Bed_Test.Transect;

                            % If there is only 1 transect, 
                            % convert it to a cell to process properly
                            [transects_mb] = convert2cell(transects_mb);

                            % Get data
                            for i = 1:length(transects_mb)

                                % MMT_Transects
                                % -------------------
                                % .Checked
                                MMT_MB_Transects.Checked(i) =...
                                    str2double(transects_mb{i}.Attributes.Checked);

                                % Files
                                % -------------------
                                files = transects_mb{i}.File;

                                % If there is only 1 files, 
                                % convert it to a cell to process
                                % properly.
                                [files] = convert2cell(files);

                                for f = 1:length(files)

                                    % .Path -> make a column cell; 
                                    % numTransects x numFiles
                                    MMT_MB_Transects.Path{i,f} =...
                                        files{f}.Attributes.PathName;

                                    % .Files -> make a column cell; 
                                    % numTransects x numFiles
                                    MMT_MB_Transects.Files{i,f} = files{f}.Text;

                                    % .Number -> make a column cell; 
                                    % numTransects x numFiles
                                    MMT_MB_Transects.Number{i,f} =...
                                        files{f}.Attributes.TransectNmb;    
                                end % for f

                                % Configuration
                                % --------------------
                                % Assign local variable
                                config_mb = transects_mb{i}.Configuration;

                                % If there is only 1 configuration, 
                                % convert it to a cell to process
                                % properly.
                                [config_mb] = convert2cell(config_mb);

                                % Find the field and active config; 
                                % attributes - > field = 0; active = 1
                                act_mb = nan(length(config_mb),1);
                                for j = 1:length(config_mb)
                                    a = str2double(config_mb{j}.Attributes.Checked);
                                    act_mb(j) = a;
                                end % for j

                                % Find indices and pick the first 
                                % occurance for field and active
                                idx_field_mb = find(act_mb == 0,1,'first');
                                idx_active_mb = find(act_mb == 1,1,'first');

                                % If there is no field configuration in the Moving Bed Test,
                                % then treat the field configuration as an active configuration.
                                % This is what is done in the previous mmt2mat.m.
                                if isempty(idx_field_mb)
                                    idx_field_mb = 1;
                                end % if isempty

                                % MMT_MB_Field_Config
                                field_config_mb = config_mb{idx_field_mb};
                                [MMT_MB_Field_Config] =...
                                    mmtconfig(field_config_mb, i, MMT_MB_Field_Config);

                                % MMT_MB_Active_Config
                                active_config_mb = config_mb{idx_active_mb};
                                [MMT_MB_Active_Config] =...
                                    mmtconfig(active_config_mb, i, MMT_MB_Active_Config);  
                            end % for i
                        end % if isfield
                    end % if strcmp
                end % for q
            end % if isfield
        waitbar(1);
        close(hwait);
        end % function mmt2mat

        %% function mmtqsum
        function [MMT_Summary] = mmtqsum(mmtDataChild,MMT_Summary)
        % This function is called from mmt2mat and accepts the data structure containing 
        % a single discharge summary for a transect. Data are passed to the function
        % through mmtDataChild and then back to the calling program through
        % MMT_Summary. MMT_Summary is also an input argument which allows arrays
        % of summary variables in the structure to be constructed through
        % multiple calls.

            % field_mmtDataChild are Index_0,Index_1,...Index_5
            fields_mmtDataChild = fieldnames(mmtDataChild);

            % For each index field, get the data
            for iDS = 1:length(fields_mmtDataChild)

                % Assign local variable
                indexfield = mmtDataChild.(fields_mmtDataChild{iDS});
                
                % Field of the index field are: UseInSummary, BeginLeft,
                % IsSubSectioned, FileName, TransectNmb, etc.
                fields_index = fieldnames(indexfield);

                % Get data - can use loop when child names are the same or similar
                for i = 1:length(fields_index)
                    
                    % Some field names do not match MMT_Summary field names and
                    % FileName filed needs to be a column cell array of 
                    if strcmp(fields_index{i}, 'UseInSummary')
                        MMT_Summary.Use(iDS) =...
                            str2double(indexfield.(fields_index{i}).Text);   
                    elseif strcmp(fields_index{i}, 'BeginLeft')
                        MMT_Summary.Begin_Left(iDS) =...
                            str2double(indexfield.(fields_index{i}).Text);
                    elseif strcmp(fields_index{i}, 'FileName')
                        MMT_Summary.FileName{iDS,1} =...
                            indexfield.(fields_index{i}).Text;
                    elseif strcmp(fields_index{i}, 'LeftEdgeSlopeCoeff')
                        MMT_Summary.LeftEdgeCoeff(iDS) =...
                            str2double(indexfield.(fields_index{i}).Text);
                    elseif strcmp(fields_index{i}, 'RightEdgeSlopeCoeff')
                        MMT_Summary.RightEdgeCoeff(iDS) =...
                            str2double(indexfield.(fields_index{i}).Text);
                    else
                        MMT_Summary.(fields_index{i})(iDS) =...
                            str2double(indexfield.(fields_index{i}).Text);
                    end % if strcmp
                end % for i
            end % for iDS
        end % function


        %% function mmtconfig
        function [MMT_Config] = mmtconfig(mmtDataChild,iTrans,MMT_Config)
        % This function is called from mmt2mat and accepts the data structure containing 
        % a single configuration for a transect. Data are passed to the function
        % through mmtDataChild and then back to the calling program through
        % MMT_Config. MMT_Config is also and input argument which allows arrays
        % of configuration variables in the structure to be constructed through multiple calls.


            % Field = 0, active = 1
            % ---------------------
            MMT_Config.Active(iTrans)=str2double(mmtDataChild.Attributes.Checked);
            
            % COMMANDS  
            % --------
            % Assign local variable
            commands = mmtDataChild.Commands;

            % Get data 
            fields_commands = fieldnames(commands);
            
            % Dynamically loop through the fields of the commands; 'Fixed
            % Commands','Fixed_Commands_StreamPro',etc. as get proper data
            for i=1:length(fields_commands)
                
                % If fields list has 'Attributes' remove it and only get the
                % Commands
                if isfield(commands.(fields_commands{i}),'Attributes');
                    commands.(fields_commands{i}) = rmfield(commands.(fields_commands{i}),'Attributes');
                end % if isfield

                % Length of each field command; 'Command_0','Command_1', etc.
                f = fieldnames(commands.(fields_commands{i}));
                if strcmp(f,'Text')
                    MMT_Config.(fields_commands{i}){iTrans} = commands.(fields_commands{i}).Text;
                else
                    for k = 1:length(f)
                        MMT_Config.(fields_commands{i}){iTrans,k} = commands.(fields_commands{i}).(f{k}).Text;
                    end % for k
                end % if strcmp
            end % for i   


            % DEPTH SOUNDER
            % -------------
            % Assign local variable
            ds = mmtDataChild.Depth_Sounder;
            
            % Get data 
            if isfield(ds,'Use_Depth_Sounder_In_Processing')
                if strcmpi(ds.Use_Depth_Sounder_In_Processing.Text,'YES');
                     MMT_Config.DS_Use_Process(iTrans)=1;
                else
                     MMT_Config.DS_Use_Process(iTrans)=0;
                end
            else
                MMT_Config.DS_Use_Process(iTrans)=-1;
            end
            MMT_Config.DS_Transducer_Depth(iTrans) =...
                str2double(ds.Depth_Sounder_Transducer_Depth.Text);
            MMT_Config.DS_Transducer_OFF(iTrans) =...
                str2double(ds.Depth_Sounder_Transducer_Offset.Text);
            
            if strcmpi(ds.Depth_Sounder_Correct_Speed_of_Sound.Text,'YES')    
                MMT_Config.DS_Cor_Spd_Sound(iTrans) =1;
            else
                MMT_Config.DS_Cor_Spd_Sound(iTrans) =0;
            end
            
            MMT_Config.DS_Scale_Factor(iTrans) =...
                str2double(ds.Depth_Sounder_Scale_Factor.Text);

            %% EXTERNAL HEADING
            % -----------------
            % Assign local variable
            eh = mmtDataChild.Ext_Heading;
            
            % Get data
            MMT_Config.Ext_Heading_Offset(iTrans) = str2double(eh.Offset.Text);

            % GPS 
            % ---
            % Assign local variable
            gps = mmtDataChild.GPS;
            
            % Get data
            MMT_Config.GPS_Time_Delay(iTrans) = str2double(gps.Time_Delay.Text);  

            % DISCHARGE 
            % ---------
            % Assign local variable
            dis = mmtDataChild.Discharge; 
            
            % Get data
            MMT_Config.Q_Top_Method(iTrans ) =...
                str2double(dis.Top_Discharge_Estimate.Text);
            MMT_Config.Q_Bottom_Method(iTrans) =...
                str2double(dis.Bottom_Discharge_Estimate.Text);
            MMT_Config.Q_Power_Curve_Coeff(iTrans) =...
                str2double(dis.Power_Curve_Coef.Text);
            MMT_Config.Q_Cut_Top_Bins(iTrans) =...
                str2double(dis.Cut_Top_Bins.Text);
            MMT_Config.Q_Bins_Above_Sidelobe(iTrans) =...
                str2double(dis.Cut_Bins_Above_Sidelobe.Text);
            MMT_Config.Q_Left_Edge_Type(iTrans) =...
                str2double(dis.River_Left_Edge_Type.Text);
            MMT_Config.Q_Left_Edge_Coeff(iTrans) =...
                str2double(dis.Left_Edge_Slope_Coeff.Text);
            MMT_Config.Q_Right_Edge_Type(iTrans) =...
                str2double(dis.River_Right_Edge_Type.Text);
            MMT_Config.Q_Right_Edge_Coeff(iTrans) =...
                str2double(dis.Right_Edge_Slope_Coeff.Text);
            MMT_Config.Q_Shore_Pings_Avg(iTrans)=...
                str2double(dis.Shore_Pings_Avg.Text);

            % EDGE ESTIMATES
            % --------------
            % Assign local variable
            ee = mmtDataChild.Edge_Estimates; 
            
            % Get data
            MMT_Config.Edge_Begin_Shore_Distance(iTrans) =...
                str2double(ee.Begin_Shore_Distance.Text);
            if strcmpi(ee.Begin_Left_Bank.Text,'YES');
                MMT_Config.Edge_Begin_Left_Bank(iTrans)=1;
            else
                MMT_Config.Edge_Begin_Left_Bank(iTrans)=0;
            end
            MMT_Config.Edge_End_Shore_Distance(iTrans)=...
                str2double(ee.End_Shore_Distance.Text);

            % OFFSETS
            % -------
            % Assign local variable
            off = mmtDataChild.Offsets;  
            fields_off = fieldnames(off);
            prefix = 'Offsets_';
            
            % Get data - can use loop when child names are the same or similar
            for i = 1:length(fields_off)
               if strcmp(fields_off{i},'ADCP_Transducer_Depth')
                   child = strcat(prefix,'Transducer_Depth');
               else
                   child = strcat(prefix,fields_off{i});
               end % if strcmp
               try
                   MMT_Config.(child)(iTrans) = str2double(off.(fields_off{i}).Text);
               catch
                   MMT_Config.(child)(iTrans) = str2double(off.(fields_off{i}).Attributes.Status);
               end % try catch
            end % for i

            % PROCESSING
            % ----------
            % Assign local variable
            pr = mmtDataChild.Processing; 
            fields_pr = fieldnames(pr);
            prefix = 'Proc_';
            
            % Get data - can use loop when child names are the same or similar
            for i = 1:length(fields_pr)
               if strcmp(fields_pr{i},'Use_3_Beam_Solution_For_BT')
                   child = strcat(prefix,'Use_3_Beam_BT');
               elseif strcmp(fields_pr{i},'Use_3_Beam_Solution_For_WT')
                   child = strcat(prefix,'Use_3_Beam_WT');
               elseif strcmp(fields_pr{i},'BT_Error_Velocity_Threshold')
                   child = strcat(prefix,'BT_Error_Vel_Threshold');
               elseif strcmp(fields_pr{i},'WT_Error_Velocity_Threshold')
                   child = strcat(prefix,'WT_Error_Vel_Threshold');
               elseif strcmp(fields_pr{i},'BT_Up_Velocity_Threshold')
                   child = strcat(prefix,'BT_Up_Vel_Threshold');
               elseif strcmp(fields_pr{i},'WT_Up_Velocity_Threshold')
                   child = strcat(prefix,'WT_Up_Vel_Threshold');
               elseif strcmp(fields_pr{i},'Fixed_Speed_Of_Sound')
                   child = strcat(prefix,'Fixed_Speed_of_Sound');
               elseif strcmp(fields_pr{i},'Mark_Below_Botom_Bad')
                   child = strcat(prefix,'Mark_Below_Bottom_Bad');
               elseif strcmp(fields_pr{i},'Use_Weighted_Mean_Depth')
                   child = strcat(prefix,'Use_Weighted_Mean');
               elseif strcmp(fields_pr{i},'Absorption')
                   child = strcat(prefix,'Absorbtion');
               else
                   child = strcat(prefix,fields_pr{i});
               end % if strcmp

               if ~isnan(str2double(pr.(fields_pr{i}).Text))
                   MMT_Config.(child)(iTrans) = str2double(pr.(fields_pr{i}).Text);
               else
                   if strcmpi(pr.(fields_pr{i}).Text,'YES')                      
                       MMT_Config.(child)(iTrans) =1;
                   else
                       MMT_Config.(child)(iTrans) =0;
                   end
                    
               end % if ~isnan
            end % for i

            % RECORDING
            % ---------
            % Assign local variable  
            rec = mmtDataChild.Recording;
            
            % Get data; put data in a row x 1 column cell array
            MMT_Config.Rec_Filename_Prefix{iTrans,1} = rec.Filename_Prefix.Text;
            MMT_Config.Rec_Output_Directory{iTrans,1} = rec.Output_Directory.Text;
            
            % Some versions prior to 2.07 do not have the Root Directory field
            if ~isempty(rec.Root_Directory.Text)
                MMT_Config.Rec_Root_Directory{iTrans,1} = rec.Root_Directory.Text;
            else
                MMT_Config.Rec_Root_Directory{iTrans,1} = nan;
            end % if ~isempty
            MMT_Config.Rec_MeasNmb{iTrans,1} =...
                rec.MeasurmentNmb.Text;
            MMT_Config.Rec_GPS{iTrans,1} =...
                rec.GPS_Recording.Text;
            MMT_Config.Rec_DS{iTrans,1} =...
                rec.DS_Recording.Text;
            MMT_Config.Rec_EH{iTrans,1} =...
                rec.EH_Recording.Text;
            MMT_Config.Rec_ASCII_Output{iTrans,1} =...
                rec.ASCII_Output_Recording.Text;
            MMT_Config.Rec_Max_File_Size{iTrans,1} =...
                str2double(rec.Maximum_File_Size.Text);
            MMT_Config.Rec_Next_Transect_Number{iTrans,1} =...
                str2double(rec.Next_Transect_Number.Text);
            MMT_Config.Rec_Add_Date_Time{iTrans,1} =...
                str2double(rec.Add_Date_Time.Text);
            MMT_Config.Rec_Use_Delimiter{iTrans,1} =...
                str2double(rec.Use_Delimiter.Text);
            MMT_Config.Rec_Delimiter{iTrans,1} =...
                rec.Custom_Delimiter.Text;
            MMT_Config.Rec_Prefix{iTrans,1} =...
                rec.Use_Prefix.Text;
            MMT_Config.Rec_Use_MeasNmb{iTrans,1} =...
                rec.Use_MeasurementNmb.Text;
            MMT_Config.Rec_Use_TransectNmb{iTrans,1}=...
                rec.Use_TransectNmb.Text;
            MMT_Config.Rec_Use_SequenceNmb{iTrans,1} =...
                rec.Use_SequenceNmb.Text;

            % WIZARD
            % ------
            % Assign local variable  
            wiz = mmtDataChild.Wizard_Info;
            
            % Get data 
            MMT_Config.Wiz_ADCP_Type(iTrans) =...
                str2double(wiz.ADCP_Type.Text);
            MMT_Config.Wiz_Firmware(iTrans) =...
                str2double(wiz.ADCP_FW_Version.Text);
            MMT_Config.Wiz_Use_Ext_Heading{iTrans} =...
                wiz.Use_Ext_Heading.Text;
            MMT_Config.Wiz_Use_GPS{iTrans} =...
                wiz.Use_GPS.Text;
            MMT_Config.Wiz_Use_DS{iTrans} =...
                wiz.Use_Depth_Sounder.Text;
            MMT_Config.Wiz_Max_Water_Depth(iTrans) =...
                str2double(wiz.Max_Water_Depth.Text);
            MMT_Config.Wiz_Max_Water_Speed(iTrans) =...
                str2double(wiz.Max_Water_Speed.Text);
            MMT_Config.Wiz_Max_Boat_Speed(iTrans )=...
                str2double(wiz.Max_Boat_Speed.Text);
            MMT_Config.Wiz_Material(iTrans) =...
                str2double(wiz.Material.Text);
            MMT_Config.Wiz_Water_Mode(iTrans) =...
                str2double(wiz.Water_Mode.Text);
            MMT_Config.Wiz_Bottom_Mode(iTrans) =...
                str2double(wiz.Bottom_Mode.Text);
            MMT_Config.Wiz_Beam_Angle(iTrans) =...
                str2double(wiz.Beam_Angle.Text);
            MMT_Config.Wiz_Pressure_Sensor{iTrans} =...
                wiz.Pressure_Sensor.Text;
            MMT_Config.Wiz_Water_Mode_13(iTrans) =...
                str2double(wiz.Water_Mode_13_Avail.Text);
            MMT_Config.Wiz_StreamPro_Default(iTrans) =...
                str2double(wiz.Use_StreamPro_Def_Cfg.Text);
            MMT_Config.Wiz_StreamPro_Bin_Size(iTrans) =...
                str2double(wiz.StreamPro_Bin_Size.Text);
            MMT_Config.Wiz_StreamPro_Bin_Number(iTrans) =...
                str2double(wiz.StreamPro_Bin_Num.Text);

            % Check to see if the following 2 fields exist.
            if isfield(wiz,'Use_GPS_Internal')
                MMT_Config.Wiz_Use_GPS_Internal(iTrans) = str2double(wiz.Use_GPS_Internal.Text);
            end % if isfield
            if isfield(wiz,'Internal_GPS_Baud_Rate_Index')
                MMT_Config.Wiz_Internal_GPS_Baud_Rate_Index(iTrans) = str2double(wiz.Internal_GPS_Baud_Rate_Index.Text);
            end % if isfield
        end % function

        %% function convert2cell
        function [dataStruct] = convert2cell(data)
            % If there is only dataStruct, convert it to a cell to process
            % properly.
            if size(data,2) == 1
                dataStruct = {data};
            else
                dataStruct = data;
            end % if size
        end % function




    
    
   