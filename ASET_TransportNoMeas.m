function [NoMeas]=ASET_TransportNoMeas(V,S,CssExtra,VelExtra,C_back,C_vel,St,Dis,Cut_vel,Cut_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculatse the Gss and Gw  unmeasured zones

%Dominguez Ruben UNL, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(VelExtra) | isempty(CssExtra)
   NoMeas.Surf.coarse = 0;
   NoMeas.Surf.fine = 0;
   NoMeas.Bottom.coarse = 0;
   NoMeas.Bottom.fine = 0;
else
    %Check the min array. Evaluate th backscatter array with velocity array
    if St.nback<St.nvel
        St.n=St.nback;
    else
        St.n=St.nvel;
    end

    if St.mback<St.mvel
        St.m=St.mback;
    else
        St.m=St.mvel;
    end

    velsurfsize=size(C_vel.zsurfvel,1);
    velbottomsize=size(C_vel.zbottomvel,1);

    backsurfsize=size(C_back.zsurfback,1);
    backbottomsize=size(C_back.zbottomback,1);

    %Control equal dimension surface

    if size(Dis.surf,1)==size(CssExtra.surf,1)
    elseif size(Dis.surf,1)>size(CssExtra.surf,1)
        %delete
        Dis.surf(1:size(Dis.surf,1)-size(CssExtra.surf,1),:)=[];
    elseif size(Dis.surf,1)<size(CssExtra.surf,1)
        CssExtra.surf(1:size(CssExtra.surf,1)-size(Dis.surf,1),:)=[];
    end


    if velsurfsize>backsurfsize
        surfsize=backsurfsize;
        cellsurf=C_vel.cellsurfVel;
    else
        surfsize=backsurfsize;
        cellsurf=C_back.cellsurfBack;
    end

    if velbottomsize<backbottomsize
        bottomsize=velbottomsize;
        cellbottom=C_vel.cellbottomVel;
    else
        bottomsize=backbottomsize;
        cellbottom=C_back.cellbottomBack;
    end

    for j=1:St.m
        if Cut_back.Backcut(1,j)<Cut_vel.Velcut(1,j)
            cut(1,j)=Cut_back.Backcut(1,j);
        else
            cut(1,j)=Cut_vel.Velcut(1,j);
        end

        if Cut_back.Backcut1(1,j)<Cut_vel.Velcut1(1,j)
            cut1(1,j)=Cut_back.Backcut1(1,j);
        else
            cut1(1,j)=Cut_vel.Velcut1(1,j);
        end

        if Cut_back.Backcut2(1,j)<Cut_vel.Velcut2(1,j)
            cut2(1,j)=Cut_back.Backcut2(1,j);
        else
            cut2(1,j)=Cut_vel.Velcut2(1,j);
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Transport Calculate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=1:St.m
         if cut(1,j)<3 %If have least 3 data per each ensemble                                          
             if cut(1,j)==0 
                NoMeas.Surf.coarse(1:surfsize,j)=0;
                NoMeas.Surf.fine(1:surfsize,j)=0;
             else                       
                for e=1:surfsize
                    %Surface interpolation% 
                      if j==1;
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(1,j))*Dis.surf(e,j);%
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);%
                      else
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(1,j))*Dis.surf(e,j);%
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);%
                      end
                end
             end

            if cut(1,j)==0
                NoMeas.Bottom.coarse(1:bottomsize,j)=0;
                NoMeas.Bottom.fine(1:bottomsize,j)=0;
            else
                for e=1:bottomsize
                   %Bottom interpolation% 
                    if j==1;
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(1,j))*Dis.bottom(e,j);%
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                    else
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(1,j))*Dis.bottom(e,j);%   
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                    end         
                end
            end

    elseif cut(1,j)>2%If have more than 3 data per each ensemble

            for e=1:surfsize%Surface
                 if e==1
               %Surface interpolation% 
                    if j==1;
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(e,j))*Dis.surf(e,j);
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);
                    else
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(e,j))*Dis.surf(e,j);
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);%
                    end
                elseif e>1
                    if j==1;
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(e,j))*Dis.surf(e,j);%
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);%
                    else
                        NoMeas.Surf.coarse(e,j)=(CssExtra.surf(e,j))*Dis.surf(e,j);%
                        NoMeas.Surf.fine(e,j)=(S.Csf/1000)*Dis.surf(e,j);%
                    end   
                end
            end

            for e=1:bottomsize%Bottom
                 if e==1
                    %Bottom interpolation%
                    if j==1;
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(e,j))*Dis.bottom(e,j);%   
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                    else
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(e,j))*Dis.bottom(e,j);%   
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                    end 
                elseif e>1
                   if j==1;
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(e,j))*Dis.bottom(e,j);%
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                   else
                        NoMeas.Bottom.coarse(e,j)=(CssExtra.bottom(e,j))*Dis.bottom(e,j);%   
                        NoMeas.Bottom.fine(e,j)=(S.Csf/1000)*Dis.bottom(e,j);%
                   end 
                end
            end
         end
    end
end

            