function  [clusters]=Corr_Clustering(A,optnumber,Labels,dis)


N = size(A,1);           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Dissim = dis;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering Analysis
Dist=squareform(Dissim);
clustTreeEuc = linkage(Dist,'weighted');        
clusters = cluster(clustTreeEuc,'MaxClust',optnumber);
%[h0,nodes0] = dendrogram(clustTreeEuc,0,'labels',Labels,'ColorThreshold','default');
%set(h0,'LineWidth',0.4)
%rotateticklabel(gca,90)




 
