function [IvrgDeg] = vergenceFromRangeXYZVec(LorRorC,IPDm,Ixyz,dim)

% function [IvrgDeg] = vergenceFromRangeXYZ(LorRorC,IPDm,Ixyz)
%
%   example call: vergenceFromRangeXYZ('L',0.065,Lxyz)
%
% vergence angle in the epipolar plane in deg given a range xyz image and an IPDm
%
% NOTE! if LorRorC -> 'L' Ixyz is expressed in a coordinate system centered at the LE nodal point
%       if LorRorC -> 'R' Ixyz is expressed in a coordinate system centered at the RE nodal point
%       if LorRorC -> 'C' Ixyz is expressed in a coordinate system centered at the CE nodal point
%
% LorRorC:      left,right or cyclopean co-ordinate system to use
%            'L' -> left      image used   as reference
%            'R' -> right     image used   as reference
%          ? 'C' -> cyclopean image ?
% IPDm:      inter-ocular distance to simulate
% Ixyz:      range data in cartesian coordinates [ r x c x 3 ]
%            Ixyz(:,:,1) -> x-values
%            Ixyz(:,:,2) -> y-values
%            Ixyz(:,:,3) -> z-values
%            NOTE: If size(Ixyz) = [1 3], automatically resized to [1 1 3]
% dim:       1   -> x vergence
% dim:       2   -> y vergence
% dim:       1:2 -> x & y vergence (default)
%%%%%%%%%%%%%%%%%%%%%%%
% IvrgDeg:   vergence angle to point

    N=size(Ixyz,1);
    if length(size(Ixyz))==2 && size(Ixyz,2) == 3, Ixyz = reshape(Ixyz,[N 1 3]); end
    if size(Ixyz,3) ~= 3, error(['vergenceFromRangeXYZ: WARNING! size(Ixyz,3) = ' num2str(size(Ixyz,4)) ' instead of 3']); end
    if ~exist('dim','var') || isempty(dim)
        dim=[1:2];
    end

    sz=[1,1,numel(dim)+1];
    switch LorRorC
    case 'L'
        LExyz = [   0 0 0];
        RExyz = [IPDm 0 0];
    case 'R'
        LExyz = [-IPDm 0 0];
        RExyz = [    0 0 0];

    case 'C'
        LExyz = [-IPDm/2 0 0];
        RExyz = [ IPDm/2 0 0];
    otherwise
        error(['vergenceAngleFromRangeXYZ: WARNING! unhandled LorRorC string: ' LorRorC]);
    end
    if numel(dim)==1 & dim==2
        C=0;
    else
        C=IPDm;
    end
    LExyz=reshape(LExyz([dim 3]),sz);
    RExyz=reshape(RExyz([dim 3]),sz);
    A=sqrt(abs(sum(bsxfun(@minus,Ixyz(:,:,[dim 3]),RExyz).^2,3)));
    B=sqrt(abs(sum(bsxfun(@minus,Ixyz(:,:,[dim 3]),LExyz).^2,3)));
    IvrgDeg = SSSdeg(A,B,C);
    negInd=Ixyz(:,:,3)<=0;
    IvrgDeg(negInd)=IvrgDeg(negInd)*-1;

end
function thetaDeg=SSSdeg(a,b,c)
    thetaDeg = acosd( (a.^2 + b.^2 - c.^2)./(2.*a.*b));
end


    % LangDeg = lawofcosinesSSSd(sqrt(sum(Lxyz.^2,3)), sqrt(sum(bsxfun(@minus,Lxyz,reshape([ IPDm 0 0],[1 1 3])).^2,3)),IPDm  );
    % RangDeg = lawofcosinesSSSd(sqrt(sum(Rxyz.^2,3)), sqrt(sum(bsxfun(@minus,Rxyz,reshape([-IPDm 0 0],[1 1 3])).^2,3)),IPDm  );

%clear
% syms V A B Ixyz RExyz LExyz IPDm d(X,Y)
% d(X,Y) =  sqrt(sum((X-Y).^2,3))
% eqn=V==acosd(d(Ixyz,LExyz).^2 + d(Ixyz,RExyz).^2 - IPDm.^2)./(2.*d(Ixyz,LExyz).*(d(Ixyz,RExyz)))
%
% eqn=V==acosd(A^2 + B^2 - IPDm.^2)./(2.*A*B)
%
% eqns=[eqn,a,b]
% solve(eqns,V)
%
%%
