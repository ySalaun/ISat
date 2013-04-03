%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 		ISat Project		%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD PARAMETERS

% load left and right pictures
ILeft	=	imread('im1.png');
IRight	=	imread('im2.png');

% load disparity map computed with IPOL article
i_dMap	=	imread('disp.png');

% noise standard deviation (to fill)
sigma	=	50;

% %% COMPUTE ERROR MAP
% 
% % window parameters
% wi      =   10;
% wj      =   10;
% phi     =   ones(2*wi+1,2*wj+1);
% phi     =   phi./sum(sum(phi));
% 
% % error map parameters
ErrMap	=	imread('disp.png');
s       =   size(ErrMap);
nrow    =   s(1);
ncol    =   s(2);
% 
% horizontal gradient picture
gradI   =   double(zeros(nrow,ncol-1));
for i=1:3
    gradI   =   gradI + double((ILeft(:,2:ncol,3)-ILeft(:,1:ncol-1,3)).^2);
end
gradI = [zeros(nrow,1) gradI];
% 
% % compute error for each pixel
% for i=1:nrow
%     disp(i)
%     for j=1:ncol
%         disparity       =   ErrMap(i,j,3);
%         error           =   (disparity ~= ErrMap(i,j,1));
%         m               =   255;
%         M               =   0;     
%         num             =   0;
%         den             =   0;
%         for ii=-wi:wi
%             for jj=-wj:wj
%                 if i+ii > 0 && i+ii <= nrow && j+jj > 0 && j+jj <= ncol && ErrMap(i,j,1) ==  ErrMap(i,j,3)
%                     m       =   min(m, DispMap(i+ii,j+jj,3));
%                     M       =   max(M, DispMap(i+ii,j+jj,3));
%                     uprim   =   gradI(i+ii,j+jj);
%                     num     =   1;
%                     den     =   den + uprim^2*phi(wi+ii+1,wj+jj+1);
%                 end
%             end
%         end
%         Aerr            =   M-m;
%         Nerr            =   sigma^2/den;
%         err             =   Aerr + Nerr;
%         ErrMap(i,j,1)   =   err;
%         ErrMap(i,j,2)   =   err;
%         ErrMap(i,j,3)   =   err;
%     end
% end
% newimage(ErrMap);

%% GRAPH CUT

% graph cut methods (cpp programm)
label=0:100;
mex GC/ISat.cpp
[o_dMap, varMap] = ISat(double(ILeft), double(IRight), double(i_dMap), double(gradI), 5, 23.7, 10);