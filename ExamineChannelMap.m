% Make a grid the size of one of the channels:
meshspace = 20;
x = 1:meshspace:512/2;
y = ones(length(x),1)*[1:meshspace:512];
% Turn the x,y coordinates of this mesh into a 2xnumpoints matrix that can
% be multiplied by a channel map:
points(1,:) = repmat(x,1,size(y,1)*size(y,2)/length(x));
points(2,:) = reshape(y,1,size(y,1)*size(y,2));

% Plot the starting mesh to check:
% figure('Position',[200 200 325 625])
% %plot(x,y,'.b')
% plot(points(1,:),points(2,:),'.b')
% xlim([-10 512/2+10])
% ylim([-10 512+10])

% Create the transformed mesh from my map (load a ChannelMapping.mat file
% first)
newSteph = A*points+repmat(b,1,size(points,2));
newMatlab = Amatlab*points+repmat(bmatlab,1,size(points,2));

% Plot and compare
figure('Position',[200 200 325 625])
plot([newMatlab(1,:);points(1,:)],[newMatlab(2,:);points(2,:)],'-b')
hold on
plot(newMatlab(1,:),newMatlab(2,:),'xr')
xlim([-10 512/2+10])
ylim([-10 512+10])
title('FitGeoTrans')

figure('Position',[400 200 325 625])
plot([newSteph(1,:);points(1,:)],[newSteph(2,:);points(2,:)],'-b')
hold on
plot(newSteph(1,:),newSteph(2,:),'xr')
xlim([-10 512/2+10])
ylim([-10 512+10])
title('CalcChannelMapping')