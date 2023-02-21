% test parametric bar predictions for Duck.  This is a 2D version of previous 1D code.  
% This code is intended to be an example illustration of how to use the 2D
% parametric beach code.
% Holman, May 2016

clear
load testCase2D     % load support data for 09/16/09 case from HLEV16 paper
                    % This includes a CRAB survey in vector and gridded
                    % form, a rectified time exposure image plus the
                    % manually digitized shoreline, sand bar and deep input
                    % data.

% create an output grid (optional)
xg = [100: 5: max(beach.deep.x)];
yg = [0: 10: 1000];

hSea = 4.5;         % offshore limit of bar envelope.  See 1D paper
[xp,yp,hp] = make2DBeachEngine( ...
        beach.shore,beach.bar,beach.deep, hSea, xg, yg);
                   
figure(1); clf; imagesc(beach.timex.y,beach.timex.x,beach.timex.I);
axis image; title(beach.txName); grid on
hold on; plot(beach.shore.y, beach.shore.x, 'w+-');
plot(beach.bar.y, beach.bar.x, 'w+-')
plot(beach.survey.y,beach.survey.x, '.')
xlabel('y (m)'); ylabel('x (m)')
    

figure(2); clf
subplot(121); imagesc(xp,yp,hp); 
axis image; set(gca,'ydir', 'nor')
cmap = flipud(colormap(jet));
colormap(cmap); colorbar; grid on; xlim([100 700])
xlabel('x (m)'); ylabel('y (m)'); caxis([0 7])
title(['Parametric beach, ' beach.bathyName(1:8)])
subplot(122); imagesc(beach.CRAB.x, beach.CRAB.y, -beach.CRAB.z)
axis image; set(gca,'ydir', 'nor')
colorbar; grid on; xlim([100 700])
xlabel('x (m)'); ylabel('y (m)'); caxis([0 7])
title(['CRAB survey, ' beach.bathyName(1:8)])
   
    