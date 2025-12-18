function s3 = plotGPCRsnake(mainChain,options)
% plotGPCRsnake Makes a snake plot of a GPCR from a chain object
%   Detailed explanation goes here

arguments
    mainChain
    options.resText = true
    options.colorData
    options.colorDataName = 'User data'
    options.colorMap = 'default'
end
% Snake plot to represent GPCRs in 2D;

% Input: mainChain object
% 

% Make 7 helices to populate

helices = mainChain.helices;
helicalResidues = [];
for i = 1:length(helices)
    helicalResidues = [helicalResidues helices(i,1):helices(i,2)];
end
% Every helical turn is populated by alternating 3 and 4 residues
figure; 
xResScatter = [];
yResScatter = [];

% Draw the backbone
for helixHere = 1:7
    resHere = helices(helixHere,1):helices(helixHere,2); % Residues on this helix
    nTurns = ceil(length(resHere)/3.5);
    % 1 Helical turn for every 2*pi
    t = 0:pi/50:(nTurns*2-0.5)*pi;
    ct = 0.5*cos(t) + 1.5*helixHere;
    % st = sin(t);
    plot(ct,-t,'Color',[0.2 0.2 0.2])
    hold on

    % Populate the backbone
    resLeft = length(resHere);
    xResHelixScatter = [];
    yResHelixScatter = [];
    for thisTurn = 1:nTurns
        yStart = (thisTurn-1)*2*pi;
        yEnd = (thisTurn-1)*2*pi + pi;
        % we will populate the [0 pi] part of every turn
        if mod(thisTurn,2)==1 % Odd, three residues per turn
            yspace = linspace(yStart,yEnd,min(5,resLeft+2));
            xspace = 0.5*cos(yspace) + 1.5*helixHere;
            resLeft = resLeft - 3;
        else % even, four residues per turn
            yspace = linspace(yStart,yEnd,min(6,resLeft+2));
            xspace = 0.5*cos(yspace) + 1.5*helixHere;
            resLeft = resLeft - 4;
        end
        xResHelixScatter = [xResHelixScatter xspace(2:end-1)];
        yResHelixScatter = [yResHelixScatter yspace(2:end-1)];
    end
    if mod(helixHere,2)==1 % Odd helices, "top" to bottom
        xResScatter = [xResScatter xResHelixScatter];
        yResScatter = [yResScatter yResHelixScatter];
    else  % even helices, bottom to top
        xResScatter = [xResScatter fliplr(xResHelixScatter)];
        yResScatter = [yResScatter fliplr(yResHelixScatter)];
    end
end

if isfield(options,'colorData')
    s3 = scatter(xResScatter,-yResScatter,250,options.colorData,'filled','MarkerEdgeColor','k');
    resText = mainChain.formatResidues(helicalResidues,'BWonly',true);
    row = [dataTipTextRow('Residue',resText) dataTipTextRow(options.colorDataName,options.colorData)];
    s3.DataTipTemplate.DataTipRows(end+1:end+2) = row;
    colormap(options.colorMap)

    % colorbar
    if any(options.colorData<0) & any(options.colorData>0) % Make colormap symmetrical if it has both negative and positive elements
       caxis([-max(abs(caxis)) max(abs(caxis))]); 
    end
else
    s3 = scatter(xResScatter,-yResScatter,250,'filled','MarkerEdgeColor','k','MarkerFaceColor','auto');
    resText = mainChain.formatResidues(helicalResidues,'BWonly',true);
    row = dataTipTextRow('Residue',resText);
    s3.DataTipTemplate.DataTipRows(end+1) = row;
end

if options.resText
    temp  = char(resText);
    text(xResScatter,-yResScatter, temp(:,1),'HorizontalAlignment','center','FontName','Arial');   
end

axis off


end