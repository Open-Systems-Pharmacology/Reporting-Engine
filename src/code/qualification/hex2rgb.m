function rgb = hex2rgb(hexColorString)

% hexColorString is always "#RRGGBB"
r_hex = hexColorString(2:3);
g_hex = hexColorString(4:5);
b_hex = hexColorString(6:7);

rgb = [colorIntensity(r_hex) colorIntensity(g_hex) colorIntensity(b_hex)];

function color = colorIntensity(hexString)

%returns color intensity between 0 and 1 from the hex string color value (between 0 and FF)
color = hex2dec(hexString) / 255.0;