function [Graph]=ASET_GraphCal(S,Cal,R,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function create the graph and calculate the statistics parameters
% for the Calibration Module

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Graph 1 ST vs 10*log(Ms2)

for t=1:size(S.Dc,1)
    [minimo,indice]=min(abs(R{t}.Depth_m(:,1)-S.depthdata(t)));
    Graph.x1(t,1)=Cal.ST{t}(indice,1);
end

Graph.x1pred=linspace(min(Graph.x1),max(Graph.x1),length(Graph.x1));

for t=1:size(S.Dc,1)
     Graph.y1(t,1)=log10(S.Css(t)/1000.*P{t}.Ph.ks);
end

g1=fitlm(Graph.x1,Graph.y1);
Graph.coef1(1)=table2array(g1.Coefficients(2,1));
Graph.coef1(2)=table2array(g1.Coefficients(1,1));
Graph.r1=g1.Rsquared.Ordinary;
Graph.y1pred=Graph.coef1(1)*Graph.x1pred+Graph.coef1(2);

        
%%
% Graph2 log10(Ms2) vs SCB

Graph.x2=10*log10(S.Css/1000);
Graph.x2pred=linspace(min(Graph.x2),max(Graph.x2),length(Graph.x2));

for t=1:size(S.Dc,1)
    [minimo,indice]=min(abs(R{t}.Depth_m(:,1)-S.depthdata(t)));
    Graph.y2(t)=Cal.SCB{t}(indice,1);
end   
Graph.y2=Graph.y2';


g2=fitlm(Graph.x2,Graph.y2);
Graph.coef2(1)=table2array(g2.Coefficients(2,1));
Graph.coef2(2)=table2array(g2.Coefficients(1,1));
Graph.r2=g2.Rsquared.Ordinary;

Graph.y2pred=Graph.coef2(1)*Graph.x2pred+Graph.coef2(2);


