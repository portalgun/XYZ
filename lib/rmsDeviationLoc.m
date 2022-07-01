function [Irms,Ibar] = rmsDeviationLoc(I,W,bPLOT,Ibar)

% function [Irms,Ibar] = rmsDeviationLoc(I,W,bPLOT,Ibar)
%
%   example call:  W = gaussKernel2D(3,[],1);
%                  Irms = rmsDeviationLoc(I,W,1);
%
% RMS deviation (local) at each point in image
%
% I:       luminance image
% W:       kernel specifying local area over
%          which to compute contrast
%          e.g. W = cosWindow([17 17]);
%               W = gaussKernel2D(4,21,1);
% bPLOT:   plot or not
%          1 -> plot
%          0 -> not
% Ibar:    local mean matched to window (optional)
%          passing Ibar speeds up computation
%          NOTE! Ibar MUST be computed correctly!!!
% %%%%%%%%%%%%%%%%%%
% Irms:    deviation image, where every point in image
%          gives the local RMS deviation

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end

% ENSURE W SUMS TO 1.0
W = W./sum(W(:));

if ~exist('Ibar','var') || isempty(Ibar) || size(Ibar,1) ~= size(I,1) || size(Ibar,2) ~= size(I,2)
    % LOCAL MEAN LUMINANCE
    Ibar   = conv2(I,W,'same');
end

% LOCAL DEVIATION AT EACH POINT IN IMAGE
Irms    = sqrt( ( conv2( I.^2,W,'same')-Ibar.^2 ) );

% PLOT STUFF
if bPLOT == 1
    % CROP PATCH FOR QUALITY CONTROL
    [P,PctrRC] = cropImageCtr(I,[],size(W'));
    % RMS DEVIATION OF INDIVIDUAL PATCH
    Prms = rmsDeviation(P-intensity(P,W),W);

    % PLOT IMAGE OF RMS DEVIATIONS
    figure('position',[1000 661 759 677]);
    subplot(3,1,[1:2]);
    imagesc(log(Irms)); hold on;
%     imagesc(I.^.4); hold on; % TO CHECK REGISTRATION
    axis image; colormap gray;
    Fig.format([],[],['RMSdev(PctrRC)=' num2str(Irms(PctrRC(1),PctrRC(2)),'%.2f')]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    plotSquare(fliplr(PctrRC),fliplr(size(W)),'y',1);

    % PLOT CROPPED PATCH
    subplot(3,3,[7]);
    imagesc(P.^.4);
    Fig.format([],[],['RMSdev=' num2str(Prms,'%.2f')]);
    caxis(minmax(I.^.4)); axis image;
    set(gca,'xtick',[]); set(gca,'ytick',[]);

    % HISTOGRAM OF RMS DEVIATION VALUES
    subplot(3,3,[8]);
    hist(Irms(:),max([21 ceil(.0002.*numel(Irms))]));
    Fig.format(['RMS'],['Histogram']);

    % HISTOGRAM OF LOGGED RMS DEVIATION VALUES
	subplot(3,3,[9]);
    hist(log(Irms(:)),max([21 ceil(.0002.*numel(Irms))]));
    Fig.format(['log(RMS)'],['Histogram']);
end
