function [my_pure_colormap] = colormapCustom(mapName,options)
%COLORMAPCUSTOM Colors Matlab plots with COOL colormaps, because the
%default Matlab ones are BOOOOOOOOOOOOOOOORING (even though one of them is
%LITERALLY called "cool")
% 
%% Usage:
% colormapCustom()
% [my_pure_colormap] = colormapCustom()
% [my_pure_colormap] = colormapCustom(colormapname, options)
%
%% Available colormaps: 
% 'default': legacy blue and red with white in the middle
%
%  'div': two colors with white in the middle, requires inputColor1 and
%  inputColor2. div works better with darker colors, so for example use: 
% [0 0.4 0.4] and [0.4 0 0.4] as input instead of [0 1 1] and [1 0 1]
%
% 'combine': Two colors with average color in the middle, requires 
%  inputColor1 and inputColor2. Unlike div, combine works better with
%  lighter colors
%
% 'seq': go from very light to  very dark passing through input color,
% requires inputColor1. Similarly, there are 'seqLight' and 'seqDark'
%
%% options:
%
% 'max_o': specifies limits of colormap in the following ways:
%
%   'div': controls how white the middle point is
%   'seq', 'seqLight' and 'seqDark': controls the light and dark limits
%
% 'inputColor1' and 'inputColor2': the color input to make the colormap
%
% 'isHex': is the input color a hex code? If yes, the function will
% transform the input colors to RGB before proceeding, defaults to false
%
% 'Nsteps': Number of steps in every direction, defaults to 128
%% Examples:
%
% colormapCustom('div','inputColor1', [0 0.4 0.4] ,'inputColor2', [0.8 0.4 0]);
%
% colormapCustom('combine','inputColor1', [0 0.8 0.8] ,'inputColor2', [0.8 0 0.8])
%
% colormapCustom('seq','inputColor1', [0 0.65 0.65])
%
% Mahdi Hijazi (2023)



    arguments
        mapName = 'default'
        options.max_o = 0.9 % Limits of color intensity (between 0 and 1)
        options.inputColor1 
        options.inputColor2 
        options.isHex = false % Is the input in Hex or RGB?
        options.Nsteps = 128 % How many steps to consider?
    end

if options.isHex % Transform from Hex to [0 1] range RGB
    colorHex(1) = options.inputColor1;
    if isfield(options,'inputColor2')
        colorHex(2) = options.inputColor2;
    end
    for i = 1:length(colorHex)
        hexHere = colorHex(i);
        if isstring(hexHere)
            hexHere = char(colorHex(i));
        end
        if strcmpi(hexHere(1,1),'#')
            hexHere(:,1) = [];
        end
        if i == 1
            options.inputColor1 = reshape(sscanf(hexHere.','%2x'),3,[]).'/255;
        elseif i == 2
            options.inputColor2 = reshape(sscanf(hexHere.','%2x'),3,[]).'/255;
        end
    end
end

max_o = options.max_o;
Nsteps = options.Nsteps;
small_step = max_o/Nsteps;
if  strcmp(mapName,'default') % Default: From pure blue to pure red  
    pure_blue_map = zeros(Nsteps,3);
    pure_red_map = zeros(Nsteps,3);
    
    for i=1:Nsteps
        pure_blue_map(i,:) = [0 0 1] + i.*[small_step small_step 0];
        pure_red_map(i,:) =  [1 max_o max_o] - i.*[0 small_step small_step];
    end
    my_pure_colormap = [pure_blue_map;1 1 1; pure_red_map];
    my_pure_colormap(my_pure_colormap<0) = 0;
    my_pure_colormap(my_pure_colormap>1) = 1;

elseif strcmp(mapName,'div') % Two colors with white in the middle

    for i =1:3
        x1(:,i) = linspace(options.inputColor1(i), max_o,Nsteps); % Go from color1 to white
        x2(:,i) = linspace(max_o, options.inputColor2(i),Nsteps); % Go from white to color2
    end
    my_pure_colormap = [x1; x2];

elseif strcmp(mapName,'combine') % Two colors with average color in the middle

    colorcombine =  (options.inputColor1 + options.inputColor2)/2;
    for i =1:3
        x1(:,i) = linspace(options.inputColor1(i), colorcombine(i),Nsteps); % Go from color1 to white
        x2(:,i) = linspace(colorcombine(i), options.inputColor2(i),Nsteps); % Go from white to color2
    end
    my_pure_colormap = [x1; x2];


elseif strcmp(mapName,'seqDark') % Going from input color to fully intense
    my_pure_colormap = zeros(2*Nsteps,3);
    for i=1:2*Nsteps
        my_pure_colormap(i,:) = options.inputColor1 - 0.5*i.*[small_step small_step small_step]; % Intensifying color 1
    end
    my_pure_colormap(my_pure_colormap<0) = 0;
    my_pure_colormap(my_pure_colormap>1) = 1;

elseif strcmp(mapName,'seqLight') % Going from very light color to input
    my_pure_colormap = zeros(2*Nsteps,3);
    for i=1:2*Nsteps
        my_pure_colormap(i,:) = options.inputColor1 + 0.5*i.*[small_step small_step small_step]; % Intensifying color 1
    end
    my_pure_colormap(my_pure_colormap<0) = 0;
    my_pure_colormap(my_pure_colormap>1) = 1;
    my_pure_colormap = flipud(my_pure_colormap);

elseif strcmp(mapName,'seq') % Going from very "light" color to fully intense going through input
    
    for i =1:3
        limit1 = max_o;
        if options.inputColor1(i) > limit1
            limit1 =  options.inputColor1(i);
        end
        limit2 = 1 - max_o;
        if options.inputColor1(i) < limit2
            limit2 =  options.inputColor1(i);
        end
        x1(:,i) = linspace(options.inputColor1(i), (limit1 - 0.25*(max(options.inputColor1) - options.inputColor1(i))),Nsteps); % Lighten it up
        x2(:,i) = linspace(options.inputColor1(i), (limit2 + 0.25*(options.inputColor1(i) - min(options.inputColor1) )),Nsteps); % Darken it up
%         x2(:,i) = linspace(options.inputColor1(i), limit2,Nsteps); % Darken it up
    end
    my_pure_colormap = [flipud(x1); options.inputColor1; x2];
end


colormap(my_pure_colormap);
end

