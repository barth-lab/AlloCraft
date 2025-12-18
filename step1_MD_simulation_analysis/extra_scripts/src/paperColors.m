% Colors to be used in plots for the paper
paperColorPalette.Hex = ["#B55475" "#ECD3DC" "#344E6D" "#CCD8E7" "#7D9263" "#DBE1D3"  ...
    "#7153A1" "#D7CEE6" "#B79315" "#F5E4A9"];
paperColorPalette.Name = ["LavenderDark" "LavenderLight" "BluegreyDark" "BluegreyLight" "OliveDark" "OliveLight"  ...
   "PurpleDark" "PurpleLight"  "GoldDark" "GoldLight"];
paperColorPalette.RGB = [    0.7098    0.3294    0.4588
    0.9255    0.8275    0.8627
    0.2039    0.3059    0.4275
    0.8000    0.8471    0.9059
    0.4902    0.5725    0.3882
    0.8588    0.8824    0.8275
    0.4431    0.3255    0.6314
    0.8431    0.8078    0.9020
    0.7176    0.5765    0.0824
    0.9608    0.8941    0.6627];

myColorLavenderGrey = colormapCustom('div','inputColor1', paperColorPalette.Hex(1) ,'inputColor2', paperColorPalette.Hex(3),'isHex',true);
myColorOlivePurple = colormapCustom('div','inputColor1', paperColorPalette.Hex(5) ,'inputColor2', paperColorPalette.Hex(7),'isHex',true);
