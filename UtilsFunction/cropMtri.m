function dimensions = cropMtri(Mtri,target)
    x = size(Mtri, 1);
    y = size(Mtri, 2);
    z = size(Mtri, 3);
    
    % Créez une nouvelle figure interactive avec App Designer
    fig = uifigure('Name', 'Crop 3D Matrix', 'Position', [100, 100, 800, 400]);

    % Créez des uisliders pour changer les tranches
    sliderxmin = uislider(fig, 'Position', [85, 40, 150, 3]);
    sliderxmax = uislider(fig, 'Position', [85, 90, 150, 3]);
    sliderymin = uislider(fig, 'Position', [360, 40, 150, 3]);
    sliderymax = uislider(fig, 'Position', [360, 90, 150, 3]);
    sliderzmin = uislider(fig, 'Position', [635, 40, 150, 3]);
    sliderzmax = uislider(fig, 'Position', [635, 90, 150, 3]);
    %for x
    sliderxmin.Limits = [1, fix(y/2)];
    sliderxmin.Value = 1;
    sliderxmax.Limits = [fix(y/2), y];
    sliderxmax.Value = y;
    %for y
    sliderymin.Limits = [1, fix(x/2)];
    sliderymin.Value = 1;
    sliderymax.Limits = [fix(x/2), x];
    sliderymax.Value = x;
    %for z
    sliderzmin.Limits = [1, fix(z/2)];
    sliderzmin.Value = 1;
    sliderzmax.Limits = [fix(z/2), z];
    sliderzmax.Value = z;
    % Créez des axes pour afficher les images
    ax1 = uiaxes(fig, 'Position', [50, 100, 250, 250]);
    ax2 = uiaxes(fig, 'Position', [320, 100, 250, 250]);
    ax3 = uiaxes(fig, 'Position', [590, 100, 250, 250]);

    % Affichez les images initiales
    
    % Ajoutez des listeners aux sliders pour mettre à jour les images
    addlistener(sliderxmin, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    addlistener(sliderxmax, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    addlistener(sliderymin, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    addlistener(sliderymax, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    addlistener(sliderzmin, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    addlistener(sliderzmax, 'ValueChanged', @(src, event) updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z));
    
    while ishandle(fig)
    dimensions =updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z);
    pause(0.5);
    end
    %croppedMatrix= cutMatrix(dimensions(3),dimensions(4),dimensions(1),dimensions(2),dimensions(5),dimensions(6),Mtri);
end

% Fonction pour mettre à jour les images en fonction de la position du slider
function dimensions =updateImages(target,sliderxmin, sliderxmax,sliderymin, sliderymax,sliderzmin, sliderzmax, Mtri, ax1, ax2, ax3, x, y, z)
    xmin_position = round(sliderxmin.Value);
    xmax_position = round(sliderxmax.Value);

    ymin_position = round(sliderymin.Value);
    ymax_position = round(sliderymax.Value);

    zmin_position = round(sliderzmin.Value);
    zmax_position = round(sliderzmax.Value);
    
    % Mise à jour du premier sous-tracé (sagittal)
    sagittal = squeeze(Mtri(:, :, target(3))); % affiche le plan du milieu de la matrice en sagittal
    imagesc(ax1, sagittal);
    title(ax1, 'Sagittal');
    xlabel(ax1, 'x');
    ylabel(ax1, 'y');
    axis(ax1, 'equal');
    axis(ax1, 'tight');
    colormap(ax1, 'gray');
    hold(ax1, 'on');
    line([0, y], [ymin_position, ymin_position], 'Color', 'green', 'LineWidth', 1, 'Parent', ax1, 'LineStyle', '--');
    line([0, y], [ymax_position, ymax_position], 'Color', 'green', 'LineWidth', 1, 'Parent', ax1, 'LineStyle', '--');
    plot(ax1, [xmin_position, xmin_position], [1, size(sagittal, 1)], 'r', 'LineWidth', 1);
    plot(ax1, [xmax_position, xmax_position], [1, size(sagittal, 1)], 'r', 'LineWidth', 1);
    hold(ax1, 'off');

    % Mise à jour du deuxième sous-tracé (frontal)
    frontal = squeeze(Mtri(:, target(2), :)); % affiche le plan du milieu de la matrice en frontal
    frontal = rot90(frontal, 1); 
    frontal = flipud(frontal);
    imagesc(ax2, frontal);
    title(ax2, 'Frontal');
    xlabel(ax2, 'y');
    ylabel(ax2, 'z');
    axis(ax2, 'equal');
    axis(ax2, 'tight');
    colormap(ax2, 'gray');
    hold(ax2, 'on');
    line([0, x], [zmin_position, zmin_position], 'Color', 'blue', 'LineWidth', 1, 'Parent', ax2, 'LineStyle', '--');
    line([0, x], [zmax_position, zmax_position], 'Color', 'blue', 'LineWidth', 1, 'Parent', ax2, 'LineStyle', '--');
    plot(ax2, [ymin_position, ymin_position], [1, size(frontal, 1)], 'g', 'LineWidth', 1);
    plot(ax2, [ymax_position, ymax_position], [1, size(frontal, 1)], 'g', 'LineWidth', 1);
    hold(ax2, 'off');

    % Mise à jour du troisième sous-tracé (transversal)
    transversal = squeeze(Mtri(target(1), :, :)); % affiche le plan du milieu de la matrice en transversal
    imagesc(ax3, transversal);
    title(ax3, 'Transversal');
    xlabel(ax3, 'z');
    ylabel(ax3, 'x');
    axis(ax3, 'equal');
    axis(ax3, 'tight');
    colormap(ax3, 'gray');    
    hold(ax3, 'on');
    line([0, z], [xmin_position, xmin_position], 'Color', 'red', 'LineWidth', 1, 'Parent', ax3, 'LineStyle', '--');
    line([0, z], [xmax_position, xmax_position], 'Color', 'red', 'LineWidth', 1, 'Parent', ax3, 'LineStyle', '--');
    plot(ax3, [zmin_position, zmin_position], [1, size(transversal, 1)], 'b', 'LineWidth', 1);
    plot(ax3, [zmax_position, zmax_position], [1, size(transversal, 1)], 'b', 'LineWidth', 1);   
    hold(ax3, 'off');
    dimensions = [xmin_position,xmax_position,ymin_position,ymax_position,zmin_position,zmax_position];
end

