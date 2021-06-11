function [center] = extract_center_from_img(img)

figure
imagesc(img)
hold on
[x,y] = getpts();

center = [y,x];

end