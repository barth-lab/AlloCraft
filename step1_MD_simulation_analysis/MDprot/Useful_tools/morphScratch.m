% Load test pdb:

%% PDB with exactly the same atoms (easy)
dirHere = 'D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\md2pathdev\pca';
[pdb, crd] = readpdb(fullfile(dirHere,'prot_highestDenpca_CA_d2-dop-WT.pdb'));

morphTraj = morph(crd(1,:),crd(2,:),'saveDir',dirHere,'saveName','morphTrajProtPCA.pdb','pdb',pdb);


%% PDB with different residues:
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\';
% dirHere = 'D:\Telework_library\dopamine_phase_3\14-d2_dop_F217M\md2pathdev\pca';
% dirHere = 'D:\Telework_library\dopamine_phase_3\16-d2_dop_I125N\md2pathdev\pca';
% dirRef = 'D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\md2pathdev\pca';
% mainPDBName = 'prot_pcaLigandReceptor_d2-dop-F217M.pdb';

% D2 DA F6.44M
config.Main.dir = 'D:\Telework_library\dopamine_phase_3\14-d2_dop_F217M\md2pathdev\pca';
config.Main.PDBName = 'prot_pcaLigandReceptor_d2-dop-F217M.pdb';
config.Main.frame = 1;
config.Main.name = 'F6.44M';

config.Ref.dir = 'D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\md2pathdev\pca';
config.Ref.PDBName = 'prot_highestDenpca_CA_d2-dop-WT.pdb';
config.Ref.frame = 1;
config.Ref.Name = 'WT';

% D2 BRC I4.46N
% config.Main.dir = 'D:\Telework_library\dopamine_phase_3\17-d2_brc_I125N\md2pathdev\pca_Test';
% config.Main.PDBName = 'Mutation environmenthighestDenpca_CA_d2-brc-I125N.pdb';
% config.Main.frame = 3;
% config.Main.name = 'I4.46N';

% D2 BRC F6.44M
% config.Main.dir = 'D:\Telework_library\dopamine_phase_3\15-d2_brc_F217M\md2pathdev\pca_Test';
% config.Main.PDBName = 'Mutation environmenthighestDenpca_CA_d2-brc-F217M.pdb';
% config.Main.frame = 1;
% config.Main.name = 'F6.44M';

% D2 BRC F6.44I
% config.Main.dir = 'D:\Telework_library\dopamine_phase_3\18b-d2_brc_F217I\md2pathdev\pca_Test';
% config.Main.PDBName = 'Mutation environmenthighestDenpca_CA_d2-brc-F217I.pdb';
% config.Main.frame = 5;
% config.Main.name = 'F6.44I';
% 
% config.Ref.dir = 'D:\Telework_library\dopamine_phase_3\3-d2_bromo_WT\md2pathdev\pca';
% config.Ref.PDBName = 'prot_highestDenpca_Distances_d2-bromo-WT.pdb';
% config.Ref.frame = 1;
% config.Ref.Name = 'WT';

% D1 DA F6.44M
% config.Main.dir = 'D:\Telework_library\gpcrdb_extra_simulations\d1_dpa\e-dpa-F228M_gp\md2pathdev\pca';
% config.Main.PDBName = 'prot_highestDenpca_CA_dpa-F228M-gp.pdb';
% config.Main.frame = 1;
% config.Main.name = 'F6.44M';

% D1 DA I4.46N
% config.Main.dir = 'D:\Telework_library\gpcrdb_extra_simulations\d1_dpa\d-dpa_I125N_gp\md2pathdev\pca';
% config.Main.PDBName = 'prot_highestDenpca_CA_dpa-I125N-gp.pdb';
% config.Main.frame = 2;
% config.Main.name = 'I4.46N';

% config.Ref.dir = 'D:\Telework_library\gpcrdb_extra_simulations\d1_dpa\a-d1r_dpa_gp\md2pathdev\pca';
% config.Ref.PDBName = 'prot_highestDenpca_Ligand_binding_d1r-dpa-gp_aligned.pdb';
% config.Ref.frame = 1;
% config.Ref.Name = 'WT';

% dirHere = 'D:\Telework_library\dopamine_phase_3\17-d2_brc_I125N\md2pathdev\pca';
% dirRef = 'D:\Telework_library\dopamine_phase_3\3-d2_bromo_WT\md2pathdev\pca';
% 
% mainPDBName = 'prot_pcaLigandReceptor_d2-dop-I125N.pdb';
% refPDBName = 'prot_highestDenpca_Distances_d2-bromo-WT.pdb';
% framePCAHere = 3;
% framePCARef = 1;

[pdbHere, crdHere] = readpdb(fullfile(config.Main.dir,config.Main.PDBName));
[pdbRef, crdRef] = readpdb(fullfile(config.Ref.dir,config.Ref.PDBName));


database = Database(databasePath);
% Chains: [receptor, G protein, ligand]
chains = 'ACB';


database.read(fullfile(config.Main.dir,config.Main.PDBName), chains, config.Main.name);
% Load the KLDiv reference structure:
database.read(fullfile(config.Ref.dir,config.Ref.PDBName), chains, config.Ref.Name );
% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

crdHere = crdHere(config.Main.frame,to3(selectname(pdbHere.chainid,chains(1))));
crdRef = crdRef(config.Ref.frame,to3(selectname(pdbRef.chainid,chains(1))));
%% Test functions

[topoAlignment,pdbAlch] = makeAlchTopo(database);

morphTrajAlch = morphAlchemical(pdbAlch,topoAlignment,crdRef(1,:),crdHere(1,:),'saveDir',config.Main.dir,'saveName','morphTrajProtPCAFrame1.pdb','nsteps',15,'smoothMerge',false,'interpolation','linear');



%% now in makeAlchTopo function
% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();
% Keep references to entries and chains used frequently
mainEntry = database.entries{1};
refEntry = database.entries{2}; % This entry spot usually reserved for active ref pdb

mainChain = mainEntry.chains{Chains.receptor};
refChain = refEntry.chains{Chains.receptor};

% Choose chain in PDB (should not be necessary!)
crdHere = crdHere(1,to3(selectname(pdbHere.chainid,chains(1))));
crdRef = crdRef(1,to3(selectname(pdbRef.chainid,chains(1))));
pdbHere = substruct(pdbHere,selectname(pdbHere.chainid,chains(1)));
pdbRef = substruct(pdbRef,selectname(pdbRef.chainid,chains(1)));

% Add mutations (from reference structure)
mutPos = find(database.residues{1}.Name(:,1) ~= database.residues{1}.Name(:,2));
% Remove any difference in length where alignment is not present, those are
% probably not mutations
mutPos(isspace(database.residues{1}.Name(mutPos,2)) | isspace(database.residues{1}.Name(mutPos,1))) = [];
mutRes = table2array(database.residues{1}(mutPos,1)); % Mutated residues
mutations = cellstr([database.residues{1}.Name(mutPos,2) num2str(mutRes) database.residues{1}.Name(mutPos,1)  ]);
mutResName = sprintf('%s2%s ',   database.residues{1}.Name(mutPos,2),database.residues{1}.Name(mutPos,1));


resIds = refChain.resIds; % Downside, only deals with current chain
resIdsPre = resIds(1:(find(resIds==mutPos)-1));
resIdsPost =resIds((find(resIds==mutPos)+1):end);
ndxPre = selectid(pdbRef.resseq,resIdsPre);
ndxPost = selectid(pdbRef.resseq,resIdsPost);

 % Should be common for both
pdbRefPre = substruct(pdbRef,ndxPre);
pdbRefPost = substruct(pdbRef,ndxPost); 

% Now get topology for Alchemical residue:
ndxResHere = selectid(pdbHere.resseq,mutPos);
ndxResRef = selectid(pdbRef.resseq,mutPos);

pdbResHere = substruct(pdbHere,ndxResHere);
pdbResRef = substruct(pdbRef,ndxResRef); 
pdbResAlch = addstruct(pdbResRef,pdbResHere);
% pdbResAlch.resname = repmat(mutResName,length(pdbResAlch.serial),1); % Rename alchemical residue

pdbAlch = addstruct(pdbRefPre,pdbResAlch);
pdbAlch = addstruct(pdbAlch,pdbRefPost);

% Get alignment:
mutRefName = pdbRef.resname(selectid(pdbRef.resseq,mutPos),:);
mutRefName = mutRefName(1,:);
MutHereName = pdbHere.resname(selectid(pdbHere.resseq,mutPos),:);
MutHereName = MutHereName(1,:);

% topoAlignment = true(length(pdbAlch.serial),3); % Columns Alch Ref Here
% topoAlignment(:,3) = not(selectname(pdbAlch.resname, mutRefName) & selectid(pdbAlch.resseq,mutPos));
% topoAlignment(:,2) = not(selectname(pdbAlch.resname, MutHereName) & selectid(pdbAlch.resseq,mutPos));

% another way to get alignments: 
% topoAlignment.Alch = false(length(pdbAlch.serial),3);
% topoAlignment.Alch(:,1)  = not(selectid(pdbAlch.resseq,mutPos)); % Columns notMut RefMut HereMut
% topoAlignment.Alch(:,2) = (selectname(pdbAlch.resname, mutRefName) & selectid(pdbAlch.resseq,mutPos));
% topoAlignment.Alch(:,3) = (selectname(pdbAlch.resname, MutHereName) & selectid(pdbAlch.resseq,mutPos));

clear topoAlignment
% With meaningful names 
topoAlignment.Alch.notMut  = not(selectid(pdbAlch.resseq,mutPos)); % Columns notMut RefMut HereMut
topoAlignment.Alch.RefMut = (selectname(pdbAlch.resname, mutRefName) & selectid(pdbAlch.resseq,mutPos));
topoAlignment.Alch.HereMut = (selectname(pdbAlch.resname, MutHereName) & selectid(pdbAlch.resseq,mutPos));
topoAlignment.Alch.BB = selectname(pdbAlch.name, 'CA', 'C', 'N', 'O');
% Ref topo
topoAlignment.Ref.notMut = not(selectid(pdbRef.resseq,mutPos));
topoAlignment.Ref.BB = selectname(pdbRef.name, 'CA', 'C', 'N', 'O');
% Here topo
topoAlignment.Here.notMut = not(selectid(pdbHere.resseq,mutPos));
topoAlignment.Here.BB = selectname(pdbHere.name, 'CA', 'C', 'N', 'O');

% Ok now unto the morphing:
% input 
% crdRef --> crdHere
NdofRef = size(crdRef,2);
NdofsHere = size(crdHere,2);
assert(NdofRef==3*sum(not(topoAlignment.Alch.RefMut)), ...
    sprintf("Ref. input coordinates and topology (topoAlignment) do not have the same size! Sizes are %i and %i",NdofRef,3*sum(not(topoAlignment.Alch.RefMut))))
assert(NdofsHere==3*sum(not(topoAlignment.Alch.HereMut)), ...
    sprintf("Input coordinates and topology (topoAlignment) do not have the same size! Sizes are %i and %i",NdofsHere,3*sum(not(topoAlignment.Alch.HereMut))))

%%

morphTrajAlch = morphAlchemical(pdbAlch,topoAlignment,crdRef(1,:),crdHere(1,:),'saveDir',dirHere,'saveName','morphTrajAlchTest5.pdb','nsteps',15,'smoothMerge',false);

%% Visualize to check:

crdHereCA = crdHere(1,to3(mainChain.getAtoms));
crdRefCA = crdRef(1,to3(refChain.getAtoms));
figure;
plot3(crdHereCA(1:3:end-2),crdHereCA(2:3:end-1),crdHereCA(3:3:end))
hold on
scatter3(crdHereCA(1:3:end-2),crdHereCA(2:3:end-1),crdHereCA(3:3:end))
% plot3(crdRefCA(1:3:end-2),crdRefCA(2:3:end-1),crdRefCA(3:3:end))
% hold on
% scatter3(crdRefCA(1:3:end-2),crdRefCA(2:3:end-1),crdRefCA(3:3:end))

axis equal

crdHereMut =  crdHere(1,to3(mainEntry.getAtoms("Chain",1,"Residues",125)));
crdRefMut =  crdRef(1,to3(refEntry.getAtoms("Chain",1,"Residues",125)));

hold on
scatter3(crdHereMut(1:3:end-2),crdHereMut(2:3:end-1),crdHereMut(3:3:end))
scatter3(crdRefMut(1:3:end-2),crdRefMut(2:3:end-1),crdRefMut(3:3:end))



%% 

M = getElementMasses(pdbAlch,topoAlignment.Alch.RefMut);

mRef = getElementMasses(pdbAlch,topoAlignment.Alch.RefMut);
mHere = getElementMasses(pdbAlch,topoAlignment.Alch.HereMut);

crd = crdRef(1,to3(not(topoAlignment.Ref.notMut)));

xyz = reshape(crd, 3, [])';

 sum(mRef.*xyz)/sum(mRef); % get COM


%% Morph functions

% Simple morph (same topology)
function morphTraj = morph(crd1,crd2,options)
    arguments
        crd1
        crd2
        options.nsteps = 10
        options.saveDir = ""
        options.saveName = []
        options.pdb = []
    end
    Ndofs = size(crd1,2);
    Ndofs2 = size(crd2,2);
    assert(Ndofs==Ndofs2, ...
        sprintf("Input coordinates do not have the same size! Sizes are %i and %i",Ndofs,Ndofs2))

    crdStep = (crd2 - crd1)/(options.nsteps-1);
    morphTraj = zeros(options.nsteps,Ndofs);

    for i = 1:options.nsteps
        morphTraj(i,:) = crd1+(i-1)*crdStep;
    end

    if ~isempty(options.saveName)
        writepdb(fullfile(options.saveDir,options.saveName), options.pdb, morphTraj);
%         writedcd(fullfile(options.saveDir,options.saveName), morphTraj)
    end
end




% Alchemical morph

function morphTraj = morphAlchemical(pdbAlch,topoAlignment,crdRef,crdHere,options)
    arguments
        pdbAlch
        topoAlignment
        crdRef
        crdHere
        options.nsteps = 10
        options.saveDir = ""
        options.saveName = []
        options.pdb = []
        options.smoothMerge = True % Show both amino acids toward the middle of the morph trajectory
        options.interpolation = 'spherical' % linear or spherical
    end
    NdofAlch = 3*length(pdbAlch.serial);
%     NdofRef = size(crdRef,2);
%     NdofsHere = size(crdHere,2);
%     assert(NdofRef==3*sum(not(topoAlignment(:,2))), ...
%         sprintf("Ref. input coordinates and topology (topoAlignment) do not have the same size! Sizes are %i and %i", ...
%         NdofRef,sum(not(topoAlignment(:,2)))))
%     assert(NdofsHere==3*sum(not(topoAlignment(:,3))), ...
%         sprintf("Input coordinates and topology (topoAlignment) do not have the same size! Sizes are %i and %i", ...
%         NdofsHere,sum(not(topoAlignment(:,3)))))
    
    if options.smoothMerge
     cond1 = options.nsteps/3;
     cond2 = 2*options.nsteps/3;
    else
     cond1 = options.nsteps/2;
     cond2 = options.nsteps/2;
    end
    
    morphTraj = zeros(options.nsteps,NdofAlch);
    % Morph the common part of the trajectory: (easy)
    crdStep = (crdHere(to3(topoAlignment.Here.notMut)) - crdRef(to3(topoAlignment.Ref.notMut)))/(options.nsteps-1);

    % Morphing of alchemical part will follow a vector from COM(ref) to
    % COM(here)
    if strcmp(options.interpolation,'linear')
        mRef = getElementMasses(pdbAlch,topoAlignment.Alch.RefMut);
        crdTemp = crdRef(1,to3(not(topoAlignment.Ref.notMut)));
        xyzTemp = reshape(crdTemp, 3, [])';
        comRef = sum(mRef.*xyzTemp)/sum(mRef); % get COM
    
        mHere = getElementMasses(pdbAlch,topoAlignment.Alch.HereMut);
        crdTemp = crdHere(1,to3(not(topoAlignment.Here.notMut)));
        xyzTemp = reshape(crdTemp, 3, [])';
        comHere = sum(mHere.*xyzTemp)/sum(mHere); % get COM
        % morphVector for linear interpolation
        morphVector = comHere - comRef;
        morphStep = morphVector/(options.nsteps-1);

    elseif strcmp(options.interpolation,'spherical')
        % Get SC COM
        mRef = getElementMasses(pdbAlch,topoAlignment.Alch.RefMut& ~topoAlignment.Alch.BB);
        crdTemp = crdRef(1,to3(not(topoAlignment.Ref.notMut)& not(topoAlignment.Ref.BB)));
        xyzTemp = reshape(crdTemp, 3, [])';
        comRef = sum(mRef.*xyzTemp)/sum(mRef); % get COM
    
        mHere = getElementMasses(pdbAlch,topoAlignment.Alch.HereMut& ~topoAlignment.Alch.BB);
        crdTemp = crdHere(1,to3(not(topoAlignment.Here.notMut)& not(topoAlignment.Here.BB)));
        xyzTemp = reshape(crdTemp, 3, [])';
        comHere = sum(mHere.*xyzTemp)/sum(mHere); % get COM

        % Try to do it with spherical interpolation: we need one more point
        % COM-BB for BOTH 
        mBBAlch = getElementMasses(pdbAlch,(topoAlignment.Alch.RefMut|topoAlignment.Alch.HereMut)&topoAlignment.Alch.BB);
        crdBBAlch = [crdRef(1,to3(not(topoAlignment.Ref.notMut)&topoAlignment.Ref.BB )) ... % Ref part
            crdHere(1,to3(not(topoAlignment.Here.notMut)&topoAlignment.Here.BB))]; % Here part
        xyzBBAlch = reshape(crdBBAlch, 3, [])';
        comBBAlch = sum(mBBAlch.*xyzBBAlch)/sum(mBBAlch); % get COM
        % now we have the three main COMs:
    
        % Do sth like this:
        % Move points to origin
        p1 = comRef - comBBAlch; 
        p2 = comHere - comBBAlch;
    
        % Transform to spherical
        [az1,elev1,r1] = cart2sph(p1(1),p1(2),p1(3));
        [az2,elev2,r2] = cart2sph(p2(1),p2(2),p2(3));
        % Interpolate in spherical space
        interAZ = linspace(az1,az2,options.nsteps);
        interELEV = linspace(elev1,elev2,options.nsteps);
        interR = linspace(r1,r2,options.nsteps);
    
        % Transform back int ocartesian space:
        [X,Y,Z] = sph2cart(interAZ,interELEV,interR);
        morphSphericalStep = [diff(X)' diff(Y)' diff(Z)'];
        % Translate points back to where they were (not needed if I only use
        % diffs)
    %     (X + dx,Y + dy,Z + dz);
    end
    
    for i = 1:options.nsteps
        morphTraj(i,to3(topoAlignment.Alch.notMut)) = crdRef(to3(topoAlignment.Ref.notMut))+(i-1)*crdStep;
       if strcmp(options.interpolation,'linear')
           if i < cond1 % Show ref and mut BB, send mut SC to SPACE
               % ref  
               morphTraj(i,to3(topoAlignment.Alch.RefMut)) = crdRef(to3(not(topoAlignment.Ref.notMut))) + ...
                    (i-1)*repmat(morphStep,1,  sum(topoAlignment.Alch.RefMut));
                
                % mut SC
                morphTraj(i,to3(~topoAlignment.Alch.BB&topoAlignment.Alch.HereMut)) = 100*ones(1, sum(to3(not(topoAlignment.Here.notMut)&~topoAlignment.Here.BB)));
                % mut BB
                morphTraj(i,to3(topoAlignment.Alch.BB & topoAlignment.Alch.HereMut)) = ...
                    crdHere(to3(not(topoAlignment.Here.notMut)&topoAlignment.Here.BB)) - ...
                    (options.nsteps - i)*repmat(morphStep,1,  sum(topoAlignment.Alch.BB &topoAlignment.Alch.HereMut));
    
           elseif i >= cond2 % Show mut and ref BB, send ref SC to SPACE 
                % mut
    %             morphTraj(i,to3(topoAlignment.Alch.RefMut)) = 100*ones(1, sum(to3(not(topoAlignment.Ref.notMut))));
                morphTraj(i,to3(topoAlignment.Alch.HereMut)) =  crdHere(to3(not(topoAlignment.Here.notMut))) - ...
                    (options.nsteps - i)*repmat(morphStep,1,  sum(topoAlignment.Alch.HereMut));
    
                % ref SC
                morphTraj(i,to3(~topoAlignment.Alch.BB&topoAlignment.Alch.RefMut)) = 100*ones(1, sum(to3(not(topoAlignment.Ref.notMut)&~topoAlignment.Ref.BB)));
                % ref BB
                morphTraj(i,to3(topoAlignment.Alch.BB & topoAlignment.Alch.RefMut)) = ...
                    crdRef(to3(not(topoAlignment.Ref.notMut)&topoAlignment.Ref.BB)) - ...
                    (options.nsteps - i)*repmat(morphStep,1,  sum(topoAlignment.Alch.BB &topoAlignment.Alch.RefMut));
           else  % Show both
                
                morphTraj(i,to3(topoAlignment.Alch.RefMut)) = crdRef(to3(not(topoAlignment.Ref.notMut))) + ...
                    (i-1)*repmat(morphStep,1,  sum(topoAlignment.Alch.RefMut));
                morphTraj(i,to3(topoAlignment.Alch.HereMut)) =  crdHere(to3(not(topoAlignment.Here.notMut))) - ...
                    (options.nsteps - i)*repmat(morphStep,1,  sum(topoAlignment.Alch.HereMut));
           end
       elseif strcmp(options.interpolation,'spherical')
            if i < cond1 % Show ref and mut BB, send mut SC to SPACE
               % ref  
               morphTraj(i,to3(topoAlignment.Alch.RefMut)) = crdRef(to3(not(topoAlignment.Ref.notMut))) + ...
                    repmat(sum(morphSphericalStep(1:i-1,:),1),1,  sum(topoAlignment.Alch.RefMut));
                
                % mut SC
                morphTraj(i,to3(~topoAlignment.Alch.BB&topoAlignment.Alch.HereMut)) = 100*ones(1, sum(to3(not(topoAlignment.Here.notMut)&~topoAlignment.Here.BB)));
                % mut BB
                morphTraj(i,to3(topoAlignment.Alch.BB & topoAlignment.Alch.HereMut)) = ...
                    crdHere(to3(not(topoAlignment.Here.notMut)&topoAlignment.Here.BB)) - ...
                    repmat(sum(morphSphericalStep(i:end,:),1),1,  sum(topoAlignment.Alch.BB &topoAlignment.Alch.HereMut));
    
           elseif i >= cond2 % Show mut and ref BB, send ref SC to SPACE 
                % mut
    %             morphTraj(i,to3(topoAlignment.Alch.RefMut)) = 100*ones(1, sum(to3(not(topoAlignment.Ref.notMut))));
                morphTraj(i,to3(topoAlignment.Alch.HereMut)) =  crdHere(to3(not(topoAlignment.Here.notMut))) - ...
                    repmat(sum(morphSphericalStep(i:end,:),1),1,  sum(topoAlignment.Alch.HereMut));
    
                % ref SC
                morphTraj(i,to3(~topoAlignment.Alch.BB&topoAlignment.Alch.RefMut)) = 100*ones(1, sum(to3(not(topoAlignment.Ref.notMut)&~topoAlignment.Ref.BB)));
                % ref BB
                morphTraj(i,to3(topoAlignment.Alch.BB & topoAlignment.Alch.RefMut)) = ...
                    crdRef(to3(not(topoAlignment.Ref.notMut)&topoAlignment.Ref.BB)) - ...
                    repmat(sum(morphSphericalStep(i:end,:),1),1,  sum(topoAlignment.Alch.BB &topoAlignment.Alch.RefMut));
           else  % Show both
                
                morphTraj(i,to3(topoAlignment.Alch.RefMut)) = crdRef(to3(not(topoAlignment.Ref.notMut))) + ...
                    repmat(sum(morphSphericalStep(1:i-1,:),1),1,  sum(topoAlignment.Alch.RefMut));
                morphTraj(i,to3(topoAlignment.Alch.HereMut)) =  crdHere(to3(not(topoAlignment.Here.notMut))) - ...
                    repmat(sum(morphSphericalStep(i:end,:),1),1,  sum(topoAlignment.Alch.HereMut));
           end

       end
    end


    if ~isempty(options.saveName)
        writepdb(fullfile(options.saveDir,options.saveName), pdbAlch, morphTraj);
%         writedcd(fullfile(options.saveDir,options.saveName), morphTraj)
    end
end

% Support functions:
function [topoAlignment,pdbAlch] = makeAlchTopo(database,options)
    arguments
        database % database object as defined in AlloDy
        options.chainName = database.entries{1}.chainNames(1);
    end

    database.align();
    % Align structures to the first database entry
    database.alignStructures();
    pdbHere = database.entries{1}.pdb;
    pdbRef = database.entries{2}.pdb;
    % Choose chain in PDB (should not be necessary!)
    pdbHere = substruct(pdbHere,selectname(pdbHere.chainid,options.chainName));
    pdbRef = substruct(pdbRef,selectname(pdbRef.chainid,options.chainName));

    % Keep references to entries and chains used frequently
    mainEntry = database.entries{1};
    refEntry = database.entries{2}; % This entry spot usually reserved for active ref pdb
    
    mainChain = mainEntry.chains{mainEntry.chainIndices(1)};
    refChain = refEntry.chains{refEntry.chainIndices(1)};
    
    % Add mutations (from reference structure)
    mutPos = find(database.residues{1}.Name(:,1) ~= database.residues{1}.Name(:,2));
    % Remove any difference in length where alignment is not present, those are
    % probably not mutations
    mutPos(isspace(database.residues{1}.Name(mutPos,2)) | isspace(database.residues{1}.Name(mutPos,1))) = [];
    mutRes = table2array(database.residues{1}(mutPos,1)); % Mutated residues
    mutations = cellstr([database.residues{1}.Name(mutPos,2) num2str(mutRes) database.residues{1}.Name(mutPos,1)  ]);
    mutResName = sprintf('%s2%s ',   database.residues{1}.Name(mutPos,2),database.residues{1}.Name(mutPos,1));
    
    
    resIds = refChain.resIds; % Downside, only deals with current chain
    resIdsPre = resIds(1:(find(resIds==mutPos)-1));
    resIdsPost =resIds((find(resIds==mutPos)+1):end);
    ndxPre = selectid(pdbRef.resseq,resIdsPre);
    ndxPost = selectid(pdbRef.resseq,resIdsPost);
    
     % Should be common for both
    pdbRefPre = substruct(pdbRef,ndxPre);
    pdbRefPost = substruct(pdbRef,ndxPost); 
    
    % Now get topology for Alchemical residue:
    ndxResHere = selectid(pdbHere.resseq,mutPos);
    ndxResRef = selectid(pdbRef.resseq,mutPos);
    
    pdbResHere = substruct(pdbHere,ndxResHere);
    pdbResRef = substruct(pdbRef,ndxResRef); 
    pdbResAlch = addstruct(pdbResRef,pdbResHere);
    % pdbResAlch.resname = repmat(mutResName,length(pdbResAlch.serial),1); % Rename alchemical residue
    
    pdbAlch = addstruct(pdbRefPre,pdbResAlch);
    pdbAlch = addstruct(pdbAlch,pdbRefPost);
    
    % Get alignment:
    mutRefName = pdbRef.resname(selectid(pdbRef.resseq,mutPos),:);
    mutRefName = mutRefName(1,:);
    MutHereName = pdbHere.resname(selectid(pdbHere.resseq,mutPos),:);
    MutHereName = MutHereName(1,:);
  
    % With meaningful names 
    topoAlignment.Alch.notMut  = not(selectid(pdbAlch.resseq,mutPos)); % Columns notMut RefMut HereMut
    topoAlignment.Alch.RefMut = (selectname(pdbAlch.resname, mutRefName) & selectid(pdbAlch.resseq,mutPos));
    topoAlignment.Alch.HereMut = (selectname(pdbAlch.resname, MutHereName) & selectid(pdbAlch.resseq,mutPos));
    topoAlignment.Alch.BB = selectname(pdbAlch.name, 'CA', 'C', 'N', 'O');
    % Ref topo
    topoAlignment.Ref.notMut = not(selectid(pdbRef.resseq,mutPos));
    topoAlignment.Ref.BB = selectname(pdbRef.name, 'CA', 'C', 'N', 'O');
    % Here topo
    topoAlignment.Here.notMut = not(selectid(pdbHere.resseq,mutPos));
    topoAlignment.Here.BB = selectname(pdbHere.name, 'CA', 'C', 'N', 'O');

    % Rename alchemical residue
    pdbAlch.resname(~topoAlignment.Alch.notMut,:) = repmat(mutResName,sum(~topoAlignment.Alch.notMut),1); 

end

function elementMass = getElementMasses(pdb,ndx)
       % Mass library according to charmm36 FF:
       m.H = 1.0080;
       m.C = 12.0110;
       m.N = 14.0070;
       m.O = 15.9990;
       m.S = 32.0600;
       
       queriedLength = length(find(ndx));
       elementMass = zeros(queriedLength,1);
       elementNames = pdb.element(ndx,:);
       elementNames = elementNames(~isspace(pdb.element(ndx,:)));

       for i = 1:queriedLength
           elementMass(i) = m.(elementNames(i));
       end
end
