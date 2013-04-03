% NEWIMAGE: Affichage d'une image dans une nouvelle figure
%           Affiche une image sans axe après rescaling des valeurs
%
% Usage: newimage(I,title,range);
%

function newimage(I,title,range);

try
    % Vérification des arguments
    if(~isreal(I))
        disp('L''image est complexe, seul le module est affiché...');
        I=abs(I);
    end
    I=double(I);
    
    if nargin>=3 && isvector(range) && isreal(range) && length(range)==2 && range(1)<range(2) % Intervalle donné?
        I=(I<range(2)).*(I>range(1)).*((I-range(1))/(range(2)-range(1)))+(I>=range(2));
    else
        % Normalisation de l'image
        D=max(I(:))-min(I(:));
        range(1)=min(I(:));
        range(2)=max(I(:));
        if D~=0
            I=(I-min(I(:)))/(max(I(:))-min(I(:)));
        else
            I=(I-min(I(:))+1)/(max(I(:))-min(I(:))+1);
        end
    end

    % Affichage de l'image
    [m,n,c] = size(I);
    f=min([max([m,n]),512])/max([m,n]);
    if(f<1)
        txt=sprintf('Affichage à un facteur %f',f);
        disp(txt);
        I=imresize(I,f);
    end
    [m,n,c] = size(I);
    figure('Units','pixels','Position',[100 100 round(n) round(m)])
    if(c<3)
        frame=zeros(m,n,3);
        frame(:,:,1)=I(:,:);
        frame(:,:,2)=I(:,:);
        frame(:,:,3)=I(:,:);
        image(frame);
    else
        image(I);
    end
    set(gca,'Visible','off');
    colormap(gray);
    %axis equal;
    
    set(gca,'Position',[0 0 1 1]);
    
    if nargin>1
        if ischar(title)
            set(gcf,'NumberTitle','off');
            title=sprintf('%s [%.2f ; %.2f]',title,range(1),range(2));
            set(gcf,'Name',title);
        end
    else
        set(gcf,'NumberTitle','off');
        title=sprintf('Image [%.2f ; %.2f]',range(1),range(2));
        set(gcf,'Name',title);
    end
    


catch
    disp('Echec de l''affichage de l''image');
end
