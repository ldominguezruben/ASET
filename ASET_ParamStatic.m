function [St]=ASET_ParamStatic(M,V,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determinate the size of the array for the different read
% format file and module used.

% by Dominguez Ruben L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M==1 %Single section, Multisection
    %Size
    [St.nvel,St.mvel]=size(V.mcsMag);

    [St.nback,St.mback]=size(V.mcsBack);

    %cell size
    for j=1:St.mvel
        for i=1:St.nvel
            if i==1
                St.cell(i,j)=V.mcsDepth(i+1,j)-V.mcsDepth(i,j);%take the same value that the second row
            else
                St.cell(i,j)=V.mcsDepth(i,j)-V.mcsDepth(i-1,j);
            end
        end
    end

elseif M==2%Calibration 

        [St.nvel,St.mvel]=size(R.Depth_m);

        [St.nback,St.mback]=size(R.Depth_m);

        %cell size
        for j=1:St.mvel
            for h=1:St.nvel
                if h==1
                    St.cell(h,j)=R.Depth_m(h+1,j)-R.Depth_m(h,j);%take the same value that the second row
                else
                    St.cell(h,j)=R.Depth_m(h,j)-R.Depth_m(h-1,j);
                end
            end
        end

end





