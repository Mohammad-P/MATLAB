chipNames = {
'Dark Skin';
'Light Skin';
'Blue Sky';
'Foliage';
'Blue Flower';
'Bluish Green';
'Orange';
'Purple Red';
'Moderate Red';
'Purple';
'Yellow Green';
'Orange Yellow';
'Blue';
'Green';
'Red';
'Yellow';
'Magenta';
'Cyan';
'White';
'Neutral 8';
'Neutral 65';
'Neutral 5';
'Neutral 35';
'Black'};
sRGB_Values = [...
115,82,68;
194,150,130;
98,122,157;
87,108,67;
133,128,177;
103,189,170;
214,126,44;
80,91,166;
193,90,99;
94,60,108;
157,188,64;
224,163,46;
56,61,150;
70,148,73;
175,54,60;
231,199,31;
187,86,149;
8,133,161;
243,243,242;
200,200,200;
160,160,160;
122,122,121;
85,85,85;
52,52,52];

% Chip width is pretty close to 14.0% on an actual, real x-rite Color Checker chart.
chipWidth = int32((imageWidth - 7 * gridWidth)/ 6);
% Image height is 4 chip widths + 5 grid widths.
imageHeight = 4 * (chipWidth + gridWidth) + gridWidth;
cols = (gridWidth + 1) : (chipWidth + gridWidth) : imageWidth;
rows = (gridWidth + 1) : (chipWidth + gridWidth) : imageHeight;
% width = max(cols); % Update this to get it accurate.
rgbImage = zeros(imageHeight, imageWidth, 3, 'uint8');
% Assign the chip colors to the proper locations in the RGB image.
chipNumber = 1;
for row = 1 : 4
	row1 = rows(row);
	row2 = row1 + chipWidth;
	for col = 1 : 6
		col1 = cols(col);
		col2 = col1 + chipWidth;
		% Assign red, green, and blue chips.
		rgbImage(row1:row2, col1:col2, 1) = sRGB_Values(chipNumber, 1);
		rgbImage(row1:row2, col1:col2, 2) = sRGB_Values(chipNumber, 2);
		rgbImage(row1:row2, col1:col2, 3) = sRGB_Values(chipNumber, 3);
		chipNumber = chipNumber + 1;
	end
end
imshow(rgbImage);
axis on;