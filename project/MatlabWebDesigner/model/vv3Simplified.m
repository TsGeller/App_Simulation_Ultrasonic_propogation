function coordinates = vv3Simplified(mri_inter)
try
x = 0;
y = 0;
z = 0;
coordinates =[x,y,z];
buttonPressed = false;
buttonexitpressed= false;
[~, ~, numSlices] = size(mri_inter);
position_plane = fix(numSlices/2);
hFigure = figure;
set(hFigure, 'Name', 'Choice point', 'color', 'w', 'WindowStyle', 'modal');
createFigureWithButton(hFigure)
choiceInput(hFigure)
colormap gray;
axis equal;
axis tight;
set(gcf, 'WindowScrollWheelFcn', @scrollCallback);
scrollDirection = 0;
while ~buttonexitpressed 
    if ~ishandle(hFigure)    
    disp('La figure n est plus affichée.');
    break;
    end
while ~buttonPressed &&~buttonexitpressed
    if ~ishandle(hFigure)    
    disp('La figure n est plus affichée.');
    break;
    end
z= stateLoop();
end
if buttonexitpressed
    break;
end
if ~ishandle(hFigure)    
    disp('La figure n est plus affichée.');
    break;
end
x=0;
y=0;
target_plane = squeeze(mri_inter(:,:,position_plane));
imagesc(target_plane);
showCross(x,y);
title('MRI picture ')
colormap turbo;
[x,y] = ginput(1);
coordinates =[fix(x),fix(y),z];
showCross(x,y);
buttonPressed= false;
end
catch exception
    disp(exception)
end

return;

    function z = stateLoop()        
            set(gcf, 'WindowScrollWheelFcn', @scrollCallback);
        if(position_plane + scrollDirection<numSlices &&position_plane + scrollDirection >0 )
            position_plane = position_plane + scrollDirection;
            z=position_plane;
        end
            target_plane = squeeze(mri_inter(:,:,position_plane));
            imagesc(target_plane);
            colormap gray;
            title('MRI picture')
            scrollDirection = 0; 
            showCross(x,y);
            drawnow;
        end
    
%%Creation d'un bouton pour arreter la navigation dans la matrice et
%%choisir un point
    function choiceInput(hFigure)       
        hButton = uicontrol('Style', 'pushbutton', 'String', 'Choice Target', ...
            'Position', [0, 5, 100, 20], 'Callback', @getInput);
        % Définissez la fonction de rappel (callback) du bouton
        function getInput(~, ~)
            x=0;
            y=0;
            showCross(x,y);
            buttonPressed = true;
            
        end
    end
%%Creation d'un bouton pour entré les valeur choisi par le curseur de
%%l'utilisateur.
    function createFigureWithButton(hFigure)
        hButton = uicontrol('Style', 'pushbutton', 'String', 'save point', ...
            'Position', [110, 5, 100, 20], 'Callback',@buttonCallback);
        % Définissez la fonction de rappel (callback) du bouton
        function buttonCallback(~, ~)
            buttonexitpressed = true;            
            disp('Exit Cliqued');
            close(hFigure);
        end
    end
%%Affiche la  derniere cross sur l'image IRM.
    function showCross(x,y)
        if(x ~=0 && y~=0)
         % Charger à nouveau l'image (pour vous assurer qu'elle est à jour)
        target_plane = squeeze(mri_inter(:,:,position_plane));
        
        % Afficher l'image
        imagesc(target_plane);
        title('MRI picture');

        % Dessiner la croix en utilisant des lignes fines
        hold on;
        line([x, x], [y - 5, y + 5], 'Color', 'red', 'LineWidth', 1);
        line([x - 5, x + 5], [y, y], 'Color', 'red', 'LineWidth', 1);
        hold off;

        drawnow;
        end
    end
%Evenement à un scroll
    function scrollCallback(~, event)
            scrollDirection = event.VerticalScrollCount;
    end
    

 end
