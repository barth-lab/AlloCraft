rCut1 = 8;

 temp = calcdistancematrix(traj(1,:));
 temp(temp>=rCut1) = 0;
 [a,b] = find(temp);
imagesc(temp);
colorbar


         [pair, dist_list] = searchrange(traj(1,:), traj(1,:), rCut1);


%%
atom1 = find(selectname(mainEntry.pdb.chainid,'B') & selectname(mainEntry.pdb.name,'N*'));
atom2 = find(selectname(mainEntry.pdb.chainid,'A') & selectname(mainEntry.pdb.name,'CG') & selectid(mainEntry.pdb.resseq,83) );

temp = calcbond(traj,[atom1 atom2]);
axisLimits = axis;
line(axisLimits(1:2),[3.5 3.5],'Color','red','LineStyle','--')
line(axisLimits(1:2),[4.5 4.5],'Color','green','LineStyle','--')


%% 
rcut = 5.5;
frameEnd = 100;
[~, ~, chainBAtomFilter, ~, runContacts] = proteinContacts(pdb, traj{1}(1:frameEnd,:),  ligandChain.resIds, ligandChain.name, mainChain.name, rcut,0);

%%
rCut1 = 5;
rCut2 = 6;
[dist, ~, chainBAtomFilter, ~, runContacts2] = proteinContactsDualCutoff(pdb, traj{1}(1:frameEnd,:),  ligandChain.resIds, ligandChain.name, mainChain.name, rCut1,rCut2,0);


figure; imagesc(runContacts{1} - runContacts2{1}) 
colorbar
figure; imagesc(runContacts2{1})
colorbar
%%
thisResNdx = find(selectid(mainEntry.pdb.resseq,153) & selectname(mainEntry.pdb.name,'CD'));

for i = 1:length(thisResNdx)

    plot(reshape(dist{1}(3,thisResNdx(i),:),[1 2501]))
    hold on
end

%%
atom1 = find(selectname(mainEntry.pdb.chainid,'B') & selectname(mainEntry.pdb.name,'O2'));
atom2 =find(selectid(mainEntry.pdb.resseq,153) & selectname(mainEntry.pdb.name,'CD'));

temp = calcbond(traj{1}(1:frameEnd,:),[atom1 atom2]);
figure
plot(temp)
hold on
axisLimits = axis;
line(axisLimits(1:2),[rCut1 rCut1],'Color','red','LineStyle','--')
line(axisLimits(1:2),[rCut2 rCut2],'Color','green','LineStyle','--')
line(axisLimits(1:2),[rcut rcut],'Color','red','LineStyle','-')
