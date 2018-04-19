
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
t=80
file=['temp' num2str(t) '.txt']
data=importdata(file);
x=data(1:end,1);
x=x*10^10;
y=data(1:end,2);
y=y*10^10;
% Make 99,000 spheres in a cube 100 by 100 by 100
bin=[0:0.5:100];
allradial=zeros(20*20,length(bin)-1);
for i=1:1:20*20
xCentroid = x(i);
yCentroid = y(i);
% Find the distances of all the particles from the centroid
distances = sort(sqrt((x - xCentroid).^2+(y - yCentroid).^2));
% Get a distribution (count) of how many particles are at each distance
re=sort(distances);
h=histogram(re,bin);
count=h.Values;
tick=bin(2:end);
radial=count./tick;
for j=1:1:length(radial)
    allradial(i,j)=radial(j);
end
end
allradial=sum(allradial);
h=plot(tick(2:end),allradial(2:end));
ylabel("G(r)/A^-1");
xlabel("r/A");
saveas(h,['r' num2str(t) '.png'],'png')