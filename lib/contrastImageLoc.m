function Irms = contrastImageLoc(Ipht,W,bPLOT)

% function contrastImageLoc(Ipht,W,bPLOT)
% 
%   example call: Irms = contrastImageLoc(Lpht,cosWindow([16 16],1),1);
% 
% RMS contrast (local) at each point in image
%
% Ipht:    luminance image
% W:       kernel specifying local area over 
%          which to compute contrast
%          e.g. W = cosWindow([17 17]);
%               W = gaussKernel2D(4,21,1);
% bPLOT:   plot or not
%          1 -> plot
%          0 -> not
% %%%%%%%%%%%%%%%%%%
% Irms:    contrast image, where every point in image 
%          gives the local RMS contrast

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end

% ENSURE W SUMS TO 1.0
W = W./sum(W(:));
% LOCAL MEAN LUMINANCE
DC   = conv2(Ipht,W,'same');
% LOCAL CONTRAST AT EACH POINT IN IMAGE
Irms    = sqrt( conv2( ((Ipht - DC)./DC).^2,W,'same') );
% SCREEN OUT BAD POINTS
Irms(isnan(Irms)) = 0;

% PLOT STUFF
if bPLOT == 1
    % PLOT IMAGE OF RMS CONTRASTS
    figure('position',[1000 661 759 677]);
    subplot(3,1,[1:2]);
    imagesc(log(Irms)); 
    axis image; colormap gray;
    formatFigure([],[],'log(RMS)');
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    
    % PLOT HISTOGRAM OF RMS CONTRASTS
    subplot(3,2,[5]);
    hist(Irms(:),max([21 ceil(.0002.*numel(Irms))])); 
    formatFigure(['RMS'],['Histogram']); 
    
    % PLOT HISTOGRAM OF LOG RMS CONTRASTS
	subplot(3,2,[6]);
    hist(log(Irms(:)),max([21 ceil(.0002.*numel(Irms))])); 
    formatFigure(['log(RMS)'],['Histogram']); 

    
end