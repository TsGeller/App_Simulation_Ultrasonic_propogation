%% Le programme principal de simulation de champ sonore ultrasonore transcrânien comprend principalement : la sélection du point cible et du point incident, l'acquisition des paramètres acoustiques du crâne, la simulation K_WAVE
%%%Version : 2022.9.21, modifiée par : Li Dapeng, utilisant le programme : AB data, macaque E (enregistré)

%% initialisation
clc;
clear;
% close all;

%% Ouvrir l'image CT (CT a été enregistré)
[filename,filepath1]=uigetfile('E:\data\配准后数据\E\CT\.dcm'); 
%Sélectionnez simplement n'importe quel fichier DIOCM pour obtenir les informations de l'image CT
filepath2=[filepath1 filename];
info1 = dicominfo(filepath2); % Lire les informations des données DICOM
[V1,spatial1,dim1] = dicomreadVolume(filepath1); 
%Lire les données, V1 est une matrice tridimensionnelle, représentant les données brutes CT (par rapport à l'enregistrement IRM)
%spatial:Résolution spatiale et coordonnées, dim : dimension:dimension
V1 = squeeze(V1);%Éliminer les cotes de longueur 1
D1=single (V1);%convertir en simple précision
clear V1;
%%Échelle de gris à la valeur CT
CT = D1.* info1.RescaleSlope + info1.RescaleIntercept;

%% Ouvrir l'image IRM
% filename2 est le nom du fichier DICOM que vous souhaitez utiliser
[filename,filepath1]=uigetfile('E:\data\配准后数据\E\MRI\.dcm'); 
filepath2=[filepath1 filename];
info2 = dicominfo(filepath2);
[V2,spatial2,dim2] = dicomreadVolume(filepath1);
%V:Données brutes (4D),spatial:Résolution spatiale et coordonnées, dim : dimension:dimension
V2 = squeeze(V2);%Éliminer les cotes de longueur 1
D2=single (V2);%convertir en simple précision
clear V2;
MRI=D2;


%% Acquérir la résolution des images CT et IRM
MRI_space=[];%Résolution tridimensionnelle de l'IRM (mm)
MRI_space(1:2)=info2.PixelSpacing;
MRI_space(3)=info2.SliceThickness;

CT_space=MRI_space;%La résolution tridimensionnelle du CT (mm), du CT et de l'IRM a été enregistrée et a la même résolution spatiale

%% Définissez les limites supérieures et inférieures et les couches de données CT du crâne (eau, tissus mous et crâne de singe), différents seuils de données sont différents
CT_th1=0;%Réglez le seuil CT de l'eau sur 0 HU, et tous ceux en dessous de 0 sont définis comme de l'eau
CT_th2=1800;%Réglez la limite inférieure du seuil CT du crâne à 1800HU, et 0-1800HU est considéré comme un tissu mou
CT_th3=7000;%Définissez la limite supérieure du seuil CT du crâne sur 7000HU, 1800-7000HU est considéré comme le crâne et le seuil CT supérieur à 7000 est défini sur 2400HU
CT(CT<CT_th1)=CT_th1;
CT(CT<CT_th2&CT>CT_th1)=CT_th2;
CT(CT>CT_th3)=CT_th3;

%% Sélectionnez le même emplacement (centre du plan horizontal) images IRM et CT pour voir l'effet d'enregistrement (facultatif)
axial_plane_center1 = round(size(CT,3)/2);
axial_plane_center2 = round(size(MRI,3)/2);
CT_im=CT(:,:,axial_plane_center1);
MRI_im=MRI(:,:,axial_plane_center2);
CTandMRI=CT_im*0.4+MRI_im;%image de superposition
figure()
subplot(131)
imagesc(CT_im)
axis equal
axis tight
colormap gray;
colorbar;
title('CT image')
subplot(132)
imagesc(MRI_im)
axis equal
axis tight
colormap gray;
colorbar;
title('MRI image')
subplot(133)
imagesc(CTandMRI)
axis equal
axis tight
colormap gray;
colorbar;
title('CT and MRI iamge')


%% Calculer le pas spatial de l'onde k de la longueur d'onde ultrasonore (la simulation tridimensionnelle contient au moins 3 longueurs d'onde par pas)
wave_speed=1500;%vitesse des vagues dans l'eau m/s
wave_frequency=0.5; %la fréquence MHz
wave_length=wave_speed/(wave_frequency*1e6)*1000;%longueur d'onde mm
ppw=5;% Combien de longueurs d'onde par point de grille（k-wave）
kwave_step=wave_length/ppw;% La taille de pas de chaque point de la grille dans k-wave



%% Interpolation des données brutes CT et IRM
[m,n,k]=size(CT);
[x,y,z] = meshgrid(1:n,1:m,1:k);
%%%Remarquer：Vq = interp3(V,Xq,Yq,Zq) Supposons une grille par défaut de points d'échantillonnage. Zone de couverture des points de grille par défaut X=1:n、Y=1:m 和 Z=1:p，其中 [m,n,p] = size(V)
[xq,yq,zq] = meshgrid(1:kwave_step/CT_space(1):m,1:kwave_step/CT_space(2):n,1:kwave_step/CT_space(3):k);
CT_inter = interp3(x,y,z,CT,xq,yq,zq);%Données CT interpolées
MRI_inter = interp3(x,y,z,MRI,xq,yq,zq);%Données IRM interpolées

%% Définissez les limites supérieure et inférieure et la stratification des données CT du crâne (eau, tissus mous et crâne de singe), et devez stratifier à nouveau après l'interpolation
CT_th1=0;%Réglez le seuil CT de l'eau sur 0 HU et définissez tout ce qui est en dessous de 0 comme de l'eau
CT_th2=1800;%Réglez la limite inférieure du seuil CT du crâne à 1800HU, et 0-1800HU est considéré comme un tissu mou
CT_th3=7000;%Définissez la limite supérieure du seuil CT du crâne sur 7000HU, 1800-7000HU est considéré comme le crâne et le seuil CT supérieur à 7000 est défini sur 2400HU
CT_inter(CT_inter<CT_th1)=CT_th1;
CT_inter(CT_inter<CT_th2&CT_inter>CT_th1)=CT_th2;
CT_inter(CT_inter>CT_th3)=CT_th3;

%% Recadrez ou amplifiez les données brutes pour réserver une place au transducteur
%%Nécessité d'introduire le transducteur pour connaître la taille d'amplification
% %augmentation des données
[m,n,k]=size(MRI_inter);%taille actuelle des données
add_number=30;%nombre de calques à ajouter
%Augmentation des données IRM
MRI_add=zeros(m+add_number,n,k);
MRI_add(add_number+1:end,:,:)=MRI_inter;
%Augmentation des données CT
CT_add=zeros(m+add_number,n,k);
CT_add(add_number+1:end,:,:)=CT_inter;
%Observer la taille de l'image après amplification
% VolumeViewer3D(MRI_add)%boîte à outils
MRI_inter=[];
MRI_inter=MRI_add;
CT_inter=[];
CT_inter=CT_add;

%% Trouvez le plan où se trouve le point cible et trouvez le reste des coordonnées en fonction de l'image du plan
%%value = VV3DModified(MRI_inter);
%La boîte à outils est requise et Coordinate correspond respectivement aux lignes, aux colonnes et aux pages de la matrice
%Trouver le plan où se trouve la cible, qu'elle soit sagittale, transversale ou coronale
val = vv3Simplified(MRI_inter);
target_plane_position=val(3);%Trouver le plan où se trouve l'ACC macaque k: 99, 140, 110, macaque E: 98,120, 120


%% Dessinez le plan où se trouve la cible
target_plane=squeeze(MRI_inter(:,:,target_plane_position));
figure
imagesc(target_plane)
axis equal
axis tight
colormap gray;

%% Déterminer les coordonnées restantes de la cible
ACC_target_xy=[val(1),val(2)];%cible
figure
imagesc(target_plane)
axis equal
axis tight
colormap gray;
hold on
plot((ACC_target_xy(1)),(ACC_target_xy(2)),'r+');

%% Calculer le contour du crâne comme sa direction tangente
target_plane_CT=squeeze(CT_inter(:,:,target_plane_position));

figure
imagesc(target_plane_CT)
axis equal
axis tight
colormap gray;


[m n]=find(target_plane_CT>1800);
[k1,av] = convhull(n,m);
x_convhull=n(k1); 
y_convhull=m(k1);
figure()
plot(n,m,'bx')
hold on
plot(x_convhull,y_convhull,'r')
axis equal
axis tight
%% Déterminer la position du transducteur et calculer l'angle d'incidence et la trajectoire
focol_distance=63.6;%La distance focale réelle du transducteur focalisé (calculée à partir du sommet du transducteur)
%%Transformez le transducteur dans le plan sagittal pour maintenir la mise au point
%%Générer un cercle avec le point cible comme centre et la distance focale comme rayon
theta = 1:0.5:360;
x0 = (ACC_target_xy(1));
y0 = (ACC_target_xy(2));
r = focol_distance/kwave_step;%combien de grilles le rayon occupe
x_circle= x0 + r*cosd(theta);
% x_circle=round(x_circle);
y_circle = y0 + r*sind(theta);
% y_circle=round(y_circle);
%%point réservé à l'extérieur du crâne
in=[];%La caractéristique de la position enregistrée, 0 signifie à l'extérieur du polygone, 1 signifie à l'intérieur et à l'extérieur du polygone
for i=1:1:length(x_circle)
    in(i)=inpolygon(x_circle(i),y_circle(i),x_convhull,y_convhull);
end
x_circle1=x_circle(in==0);%obtenir le point intérieur
y_circle1=y_circle(in==0);%obtenir le point intérieur

%% Limiter les points à sélectionner et exclure les points déraisonnables
x_circle1_cut=[];
y_circle1_cut=[];
i=find((x_circle1>50 & x_circle1<190 & y_circle1<ACC_target_xy(2)));
%i=find((x_circle1>50 & x_circle1<190));
x_circle1_cut= x_circle1(i);%obtenir le point intérieur
y_circle1_cut=y_circle1(i);%obtenir le point intérieur


%% dessin
figure()
imagesc(target_plane_CT)
xlabel('z distance ')
ylabel('y distance ')
axis equal
axis tight
colormap gray;
title('CT Sagittal plane and target')
hold on
plot((ACC_target_xy(1)),(ACC_target_xy(2)),'r+')
hold on
plot(x_circle1_cut,y_circle1_cut,'b+')
hold on
plot(x_convhull,y_convhull,'r')

%% Obtenir les informations sur le crâne, les coordonnées x, y et la valeur HU du faisceau sonore à travers le crâne (modifié)
pos_final_all={};%L'ensemble des coordonnées de tous les points passant par le crâne
skull_profile_all={};%Une collection de pixels avec tous les points passant par le crâne
angle_inside_all=[];%Ensemble des angles de tous les points passant par le crâne
x2=x0;y2=y0;%coordonnées cibles
for i1=1:length(x_circle1_cut)
% x1=199;y1=121;
x1=x_circle1_cut(i1);%Coordonnées du sommet du transducteur
y1=y_circle1_cut(i1);

%%Calculer les coordonnées en pixels de l'image traversées par la ligne
pos1=[];
pos2=[];
pos3=[];
pos_final=[];
pos1 = [x1,y1;x2,y2];%Obtenir les coordonnées du point final de la ligne
dx=pos1(2,1)-pos1(1,1);%La première colonne est l'abscisse
dz=pos1(2,2)-pos1(1,2);%La deuxième colonne est la coordonnée y
if dz==0
    pos2(:,1)=min(pos1(1,1),pos1(2,1)):max(pos1(1,1),pos1(2,1));
    pos2(:,2)=pos1(2,2);
    pos_final=round(pos2);
elseif dx==0
    pos2(:,2)=min(pos1(1,2),pos1(2,2)):max(pos1(1,2),pos1(2,2));
    pos2(:,1)=pos1(1,1);
    pos_final=round(pos2);
else
    k=dz/dx;
    b=pos1(1,2)-pos1(1,1)*k;
    %générer une ligne：(1)z=kx+b;(2)x=(z-b)/k;
%Calculer avec l'abscisse
x11=min(pos1(2,1),pos1(1,1));x22=max(pos1(2,1),pos1(1,1));
    pos2(:,1)=x11:0.25:x22;
    pos2(:,2)=pos2(:,1).*k+b;
    pos2=round(pos2);
    pos2=unique(pos2,'rows','stable');
%Calculé avec les coordonnées z
z11=min(pos1(1,2),pos1(2,2));z22=max(pos1(1,2),pos1(2,2));
    pos3(:,2)=z11:0.25:z22;
    pos3(:,1)=(pos3(:,2)-b)./k;
    pos3=round(pos3);
    pos3=unique(pos3,'rows','stable');
    if size(pos2,1)>=size(pos3,1)
    pos_final=pos2;
    else
      pos_final=pos3;  
    end
end
clear pos1 pos2 pos3 
if pos_final(end,:)~=[x2 y2]
pos_final=flip(pos_final);
end
pos_final_all{i1}=pos_final;
%%Obtenez le contour du crâne en fonction des pixels passés par la ligne droite et trouvez le point d'entrée du crâne
skull_profile=[];
for i2=1:1:length(pos_final)
    skull_profile(i2)=target_plane_CT(pos_final(i2,2),pos_final(i2,1));
end
skull_profile_all{i1}=skull_profile;
%%Calculer quelques informations crâniennes (épaisseur, porosité, variation de vitesse, etc.)
a=find(skull_profile~=0);%Trouver tous les points du crâne qui ne sont pas 0
x_inside=pos_final(a(1),1);
y_inside=pos_final(a(1),2);%Trouver les coordonnées du point incident

%%Calculer l'angle entre deux droites
%%Trouver les coordonnées des extrémités tangentes correspondantes des deux côtés du point incident
for i3=1:1:length(x_convhull)-1
    if (x_inside-x_convhull(i3))*(x_inside-x_convhull(i3+1))<=0 && (y_inside-y_convhull(i3))*(y_inside-y_convhull(i3+1))<=0
        x3=x_convhull(i3);
        y3=y_convhull(i3);
        x4=x_convhull(i3+1);
        y4=y_convhull(i3+1);
    elseif (x_inside-x_convhull(i3))*(x_inside-x_convhull(i3+1))<=0 && y_convhull(i3)==y_convhull(i3+1)
        x3=x_convhull(i3);
        y3=y_convhull(i3);
        x4=x_convhull(i3+1);
        y4=y_convhull(i3+1);
    elseif (y_inside-y_convhull(i3))*(y_inside-y_convhull(i3+1))<=0 && x_convhull(i3)==x_convhull(i3+1)
        x3=x_convhull(i3);
        y3=y_convhull(i3);
        x4=x_convhull(i3+1);
        y4=y_convhull(i3+1);
        continue
    end
end

%%Calculez l'angle en utilisant la formule du produit scalaire des vecteurs a.b=|a|*|b|*con(l);
v1 = [x1,y1] - [x2,y2];
v2 = [x4,y4] - [x3,y3];
angle_inside=acos(dot(v1,v2)/(norm(v1)*norm(v2)))*(180/pi);
%Angle de rotation radian* (180/pi)
angle_inside_all(i1)=angle_inside;
end

%% Dessinez un histogramme de la distribution des angles
figure()
histogram(angle_inside_all)
xlabel('Angle range  ')
ylabel('Number ')
title('Potential ultrasound angle of incidence')
%% Trouvez la position du sommet et l'angle de l'angle d'incidence le plus proche de l'ange ° (x est l'angle requis, tel que 90 degrés)
angel=90;
x_angel=[];
y_angel=[];
angle_angel=[];
j=1;
diff_angel=abs(angle_inside_all-angel);
k=find(diff_angel==min(diff_angel));
x_angel=x_circle1_cut(k);
y_angel=y_circle1_cut(k);
angle_angel=angle_inside_all(k);
%%画图
figure()
imagesc(target_plane_CT)
xlabel('z distance ')
ylabel('y distance ')
axis equal
axis tight
colormap gray;
title(['CT iamge and target，','The angle of ultrasound incidence is ',num2str(angle_angel),'°'])
hold on
plot((ACC_target_xy(1)),(ACC_target_xy(2)),'r+')
hold on
plot(x_circle1_cut,y_circle1_cut,'b+')
hold on
plot(x_convhull,y_convhull,'r')
hold on
plot([x_angel(1) x0],[y_angel(1),y0],'b')
% view(-90,90);%Champ de vision de l'image pivoté de 90°

% %% Trouver le profil transcrânien pour un angle d'incidence angel° (l'incidence oblique n'est pas précise)
skull_profile_1=skull_profile_all{k};
skull_profile_x=0:kwave_step:(length(skull_profile_1)-1)*kwave_step;%mm
skull_profile_2=skull_profile_1(skull_profile_1>1800);
Skull_thickness=(length(skull_profile_2)-1)*kwave_step;
figure()
plot(skull_profile_x,skull_profile_1)
xlabel('distance mm')
ylabel('HU ')
title(['Skull profile on ultrasound path,Skull thickness is ',num2str(Skull_thickness),' mm'])


%% Tracez une ligne droite reliant le point cible et le point incident, et déterminez l'angle de rotation (en rendant le faisceau sonore perpendiculaire à l'image pour une observation facile)
transducer_xy=[x_angel(1),y_angel(1)];%cibler
v1=ACC_target_xy;%cibler
v2=transducer_xy;%sommet du transducteur
v3=[];%pointe verticale
v3(1)=ACC_target_xy(1);v3(2)=transducer_xy(2);
distance_tran_tar=round(norm(v2-v1));%Distance entre le sommet du transducteur et le point cible (maille)
%Calculer l'angle entre la direction du faisceau sonore ultrasonique et la normale
rotate_angle=acos((norm(v3-v1)/norm(v2-v1)))*(180/pi);
figure
imagesc(target_plane_CT)
axis equal
axis tight
colormap gray;
hold on
plot([v1(1),v2(1)],[v1(2),v2(2)],'b');
hold on
plot([v1(1),v3(1)],[v1(2),v3(2)],'r');
hold on
plot(x_circle1_cut,y_circle1_cut,'b+')
hold on
plot(x_convhull,y_convhull,'r')
title(['MRI image include target，','rotate angle is ',num2str(rotate_angle),'°'])

%% Déterminer les coordonnées du point cible et les coordonnées du point incident ultrasonore à partir de l'image planaire contenant le point cible, voici un exemple
transducer_position=[v2(2),v2(1),target_plane_position];
target_position=[v1(2),v1(1),target_plane_position];
transform_matrx=MRI_inter*0;%Utilisé pour enregistrer les sommets et les cibles du transducteur
transform_matrx(target_position(1)-2:target_position(1)+2,target_position(2)-2:target_position(2)+2,target_position(3)-2:target_position(3)+2)=200;

% VolumeViewer3D(transform_matrx)

%% Faites pivoter les données en fonction de l'angle calculé (ici, faites pivoter autour de l'axe z, volumeViewer affiche la zone cible dans la tranche xy)

MRI_rotate=imrotate3(MRI_inter,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);
CT_rotate=imrotate3(CT_inter,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);
transform_matrx_rotate=imrotate3(transform_matrx,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);

%% Observez les données tournées
% VolumeViewer3D(MRI_rotate);%La boîte à outils est requise et Coordinate correspond respectivement aux lignes, aux colonnes et aux pages de la matrice

%% Découpez les données (y compris CT, IRM et matrice de position) pour réduire la quantité de calculs
MRI_final=[];
CT_final=[];
transform_matrx_final=[];
row_set=[30,310];
column_set=[80,250];
page_set=[30,170];
MRI_final=MRI_rotate(row_set(1):row_set(2),column_set(1):column_set(2),page_set(1):page_set(2));
CT_final=CT_rotate(row_set(1):row_set(2),column_set(1):column_set(2),page_set(1):page_set(2));
transform_matrx_final=transform_matrx_rotate(row_set(1):row_set(2),column_set(1):column_set(2),page_set(1):page_set(2));
% VolumeViewer3D(MRI_final);
% VolumeViewer3D(CT_final);
% VolumeViewer3D(transform_matrx_final);

%% Trouver les sommets et les points cibles des transducteurs transformés
Lax=find(transform_matrx_final==200);%
s=size(transform_matrx_final);%Calculer la taille d'un tableau tridimensionnel
[m,n,k]=ind2sub(s,Lax(round(length(Lax)/2)));%Convertir un indice unique maximum en indice multidimensionnel tridimensionnel
targrt_position_new=[m,n,k];
transducer_position_new=[m-distance_tran_tar,n,k];


%% Observez la position relative de la cible de stimulus transformée et du point incident du transducteur dans le nouveau plan
v11=[targrt_position_new(2),targrt_position_new(1)];%cible
v22=[transducer_position_new(2),transducer_position_new(1)];%sommet du transducteur
figure
imagesc(squeeze(MRI_final(:,:,targrt_position_new(3))))
axis equal
axis tight
colormap gray;
hold on
plot([v11(1),v22(1)],[v11(2),v22(2)],'r');
title(['MRI image'])

CT_final_Sagittla=squeeze(CT_final(:,:,targrt_position_new(3)));
figure
imagesc(CT_final_Sagittla)
axis equal
axis tight
colormap gray;
hold on
plot([v11(1),v22(1)],[v11(2),v22(2)],'r');
title(['CT image in Sagittla Plan'])

%% Calculer l'épaisseur du crâne et la distribution des valeurs CT à travers lesquelles passe le faisceau sonore
Transcranial_HU=CT_final_Sagittla(v22(2):v11(2),v11(1));
position_skll=(0:abs(v11(2)-v22(2))).*kwave_step;

skull_index=find(Transcranial_HU>1800);%Trouver le point de coupure du crâne
skull_length=abs((skull_index(1)-skull_index(end)))*kwave_step;

figure()
plot(position_skll,Transcranial_HU)
xlabel('Ditance (mm)','FontSize',16);
ylabel('HU','FontSize',16);
title(['The thickness of the skull crossed by the ultrasound is ',num2str(skull_length),'mm'],'FontSize',16)


%% Définissez les limites supérieure et inférieure et la stratification (eau, tissus mous et crâne) des données CT du crâne, et devez stratifier à nouveau après l'interpolation
CT_th1=0;%Réglez le seuil CT de l'eau sur 0 HU et définissez tout ce qui est en dessous de 0 comme de l'eau
CT_th2=1800;%Réglez la limite inférieure du seuil CT du crâne à 150HU, et 0-150HU est considéré comme un tissu mou
CT_th3=7000;%Définissez la limite supérieure du seuil CT du crâne sur 2400HU, 150-2400HU est considéré comme le crâne et le seuil CT supérieur à 2400 est défini sur 2400HU
CT_final(CT_final<CT_th1)=CT_th1;
CT_final(CT_final<CT_th2&CT_final>CT_th1)=1800;
CT_final(CT_final>CT_th3)=CT_th3;

%% définir les constantes acoustiques
HU_max=max(CT_final(:));
HU_water=0;
swater=wave_speed;%Règle la vitesse du son dans l'eau, correspondant au minimum HU
dwater=1000;%densité de l'eau
stissue=1560;%vitesse du son des tissus mous
dtissue=1030;%densité des tissus mous
sbone=3100;%Supposée être la vitesse du son dans l'os cortical, correspondant au maximum HU
dbone=1900;%Régler sur la densité dans l'os cortical, correspondant au plus grand HU
f=wave_frequency;%Tester la fréquence ultrasonique pour ce coefficient d'atténuation (MHz)
b=1.1;%Coefficient d'absorption a=a0*f^y ; a0 doit être transformé en [dB/(MHz^y cm)]
amin=0.2;%Attenuation  [dB/(MHz^y cm)]
amax=8;%%Attenuation  [dB/(MHz^y cm)]
atissue=0.6;%Attenuation  [dB/(MHz^y cm)]
%alpha_dB = 20*log10 （exp（alpha_Neper）） 
%alpha_dB= alpha_Neper * 20*log10 （exp） = 8.686 * alpha_Neper

%% Calcul des constantes de champ sonore, y compris la vitesse, la densité et l'atténuation du son
axial_plane_center = round(size(CT_final,3)/2);
% figure()
% imagesc(CT_final(:,:,axial_plane_center))
% axis equal
% axis tight
% colormap;
% h=colorbar;
% set(get(h,'title'),'string','HU');

%Utilisez HU pour former une relation linéaire avec la densité et la vitesse pour résoudre
density=(CT_final-HU_water).*((dbone-dwater)/(HU_max-HU_water))+dwater;%Calculer la densité
speed=(CT_final-HU_water).*((sbone-swater)/(HU_max-HU_water))+swater;%Calculer la vitesse
k=(dbone-density)./(dbone-dwater);%Calculer la porosité
atta=amin+(amax-amin).*(k.^0.5);%Calculer le facteur d'atténuation
density(CT_final==1800)=dtissue;%densité des tissus mous
speed(CT_final==1800)=stissue;%vitesse des tissus mous
atta(CT_final==1800)=atissue;%coefficient d'atténuation des tissus mous
atta(CT_final==0)=amin;%Régler l'atténuation de l'eau sur amin

%% Afficher les résultats de calcul des paramètres acoustiques (facultatif)
axial_plane_center = round(size(CT_final,3)/2);
HU_im=CT_final(:,:,axial_plane_center);
density_im=density(:,:,axial_plane_center);
speed_im=speed(:,:,axial_plane_center);
atta_im=atta(:,:,axial_plane_center);
figure()
subplot(141)
imagesc(HU_im)
axis equal
axis tight
colormap;
h=colorbar;
set(get(h,'title'),'string','HU');
title('CT image')
subplot(142)
imagesc(density_im)
axis equal
axis tight
colormap;
h=colorbar;
set(get(h,'title'),'string','kg/m^3');
title('Density')
subplot(143)
imagesc(speed_im)
axis equal
axis tight
colormap;
h=colorbar;
set(get(h,'title'),'string','m/s');
title('Speed')
subplot(144)
imagesc(atta_im)
axis equal
axis tight
colormap;
h=colorbar;
set(get(h,'title'),'string','dBMHz^−^1^.^1cm^−^1');
title('Attenuation')


%% Définir les conditions de simulation de l'onde k %%

% Définir la grille de calcul
[Nx,Ny,Nz]=size(CT_final);
dx = kwave_step/1e3;    	% La taille du pas dans la direction x
dy = dx ;            % La taille du pas dans la direction y
dz=dx ; % La taille du pas dans la direction z
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%Définir les propriétés du support
medium.sound_speed = swater * ones(Nx, Ny, Nz);	% [m/s]
medium.density = dwater * ones(Nx, Ny, Nz);       % [kg/m^3]
medium.alpha_coeff = 0 * ones(Nx, Ny, Nz);  % [dB/(MHz^y cm)]
medium.alpha_power=b;

%définir la couche d'atténuation
medium.sound_speed = speed;	% [m/s]
medium.density = density;       % [kg/m^3]
medium.alpha_coeff=atta; % [dB/(MHz^y cm)]

source_f0       = wave_frequency*1e6;      % source frequency [Hz]
%Définir la grille horaire
% cfl             = 0.5;      % CFL Numéro
cfl             = 0.2;      % CFL Numéro
PPP = round(ppw / cfl);%Combien de cycles par point
dt = 1 / (PPP * source_f0);
%(1)Définir automatiquement l'heure
% kgrid.makeTime(medium.sound_speed);
%(2)Définir manuellement l'heure
% cfl             = 0.5;      % CFL number
max_distance=sqrt(Nx^2+Ny^2+Nz^2)*dx;%La distance maximale de propagation des ultrasons
t_end           = max_distance/1500;    % Temps de calcul total (doit être supérieur au temps de stabilisation du système)
% PPP = round(ppw / cfl);%Calculer combien de points par période
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

%définir la source
% Définir les propriétés du transducteur
source_f0       = wave_frequency*1e6;      % source frequency [Hz]
source_roc      = 66e-3;    % Rayon de courbure du transducteur [m]
source_diameter = 64e-3;    % Diamètre effectif du transducteur [m]
source_amp      = 1e6;      % La pression acoustique de surface du transducteur [Pa]

%Définit le type de transducteur, ici un bol de focalisation
% Définir la position spatiale et la rotation du bol
bowl_pos = [kgrid.x_vec(transducer_position_new(1)), kgrid.y_vec(transducer_position_new(2)), kgrid.z_vec(transducer_position_new(3))];
focus_pos = [kgrid.x_vec(targrt_position_new(1)), kgrid.y_vec(targrt_position_new(2)), kgrid.z_vec(targrt_position_new(3))];

% définir un vide kWaveArray
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 10);
% Augmenter le transducteur en forme de bol
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);

%% définition binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);
% Définir la série temporelle du signal
source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);
% Définir la source du signal, la redistribution du signal, cette étape demande beaucoup de temps et une grande mémoire
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
% VolumeViewer3D(int8(source.p_mask));
% 
% VolumeViewer3D(medium.sound_speed);
% VolumeViewer3D(medium.density);
% VolumeViewer3D(medium.alpha_coeff);
% définir le capteur
sensor.mask = [1,1,1, Nx, Ny,Nz].';%La grille à enregistrer est la grille entière 
sensor.record = {'p_max'};%Le paramètre à enregistrer est la valeur maximale de la pression acoustique
% sensor.record = {'p_min'};%Le paramètre à enregistrer est la valeur minimale de la pression acoustique
% N'enregistrer que des valeurs stables
record_periods  = 3;        % nombre de cycles à enregistrer
sensor.record_start_index = kgrid.Nt - record_periods * PPP + 1;

%% Définir les conditions de sortie de la simulation
% display_mask = source.p_mask;%Afficher la grille des transducteurs
display_mask = CT_final*0;%montrer le crâne
display_mask(CT_final>1800)=1;
% VolumeViewer3D(display_mask);
input_args = {'DisplayMask', display_mask, 'PlotLayout', true, 'PMLInside', false, 'PlotPML', false,...
    'DataCast', 'single', 'DataRecast', true,'PlotScale', ...
    [-1, 1] * source_amp};%cpu opération
% input_args = {'DisplayMask', display_mask, 'RecordMovie',true,'MovieName','1','PlotLayout', true,  'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10},...
%     'PMLInside', false, 'PlotPML', false,...
%     'DataCast', 'gpuArray-single', 'DataRecast', true,'PlotScale', ...
%     [-1, 1] * source_amp};%gpu opération,et enregistrer une vidéo
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});%Démarrer la simulation d'onde k 3D

%% Évaluation des résultats de simulation %%
%% Extraire l'amplitude du signal du capteur
% amp = extractAmpPhase(sensor_data.p_max, 1/kgrid.dt, source_f0, ...
%     'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);%La deuxième dimension représente le temps
% % Matrice de reconstruction
% amp = reshape(amp, Nx, Ny,Nz);
 amp=sensor_data.p_max;


%% Obtenir les coordonnées 3 axes
x_vec = kgrid.x_vec;
y_vec = kgrid.y_vec;
z_vec = kgrid.z_vec;

%% Trouver la mise au point, calculer la profondeur de mise au point
% %Exclure les points de valeur maximale sur le support
Estimated_focal_length=10/1000;%distance de mise au point visuelle
start_n=round(Estimated_focal_length/dx);
amp_exculude=amp;
amp_exculude(1:start_n,:,:)=0;

amp_max=max(amp_exculude(:));%Calculer la valeur maximale d'un tableau tridimensionnel,Pa
s=size(amp);%Calculer la taille d'un tableau tridimensionnel
Lax=find(amp==amp_max);%Calculer l'indice unique de la position maximale
[m,n,k]=ind2sub(s,Lax);%Convertir un indice unique maximum en indice multidimensionnel tridimensionnel
Loc_focus=[m,n,k];%indice de la valeur maximale
focus_coordinate=[x_vec(m),y_vec(n),z_vec(k)];
focal_length=norm(focus_coordinate-bowl_pos)*1000;%Calculer la distance du point focal au sommet du transducteur,mm

%% Dessinez une image contenant les sommets du transducteur et la cible
field_im1=amp(:,:,targrt_position_new(3));


% Un graphique contenant uniquement la distribution du champ sonore
figure;
imagesc(1e3 * x_vec, 1e3 * y_vec, field_im1'/1e6);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Pressure Field');
h=colorbar;
set(get(h,'title'),'string','MPa');
%Le résultat est tourné de 90 degrés dans le sens inverse des aiguilles d'une montre
view(-90,-90)

% Contient la distribution du champ sonore ainsi que les points cibles et incidents et les médias
figure;
imagesc(1e3 * x_vec, 1e3 * y_vec, field_im1'/1e6);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title(['Transcranial Pressure Field at real focal point is ',num2str(amp_max/1e6),'MPa'],'FontSize',16)
h=colorbar;
set(get(h,'title'),'string','MPa');
hold on
plot(bowl_pos(1)*1000,bowl_pos(2)*1000,'r*')
hold on
plot(focus_coordinate(1)*1000,focus_coordinate(2)*1000,'r*')
hold on
plot([bowl_pos(1)*1000,focus_coordinate(1)*1000],[bowl_pos(2)*1000,focus_coordinate(2)*1000],'r')
view(-90,-90)


%% Dessiner des images comparatives de CT, IRM et champ sonore
amp_focus_coordinate_oral=field_im1(targrt_position_new(1),targrt_position_new(2));%Trouver la pression acoustique du point cible d'origine（Pa）
focus_coordinate_oral=[x_vec(targrt_position_new(1)),y_vec(targrt_position_new(2)),z_vec(targrt_position_new(3))];
figure
subplot(1,3,1)
imagesc(squeeze(MRI_final(:,:,targrt_position_new(3))))
axis equal
axis tight
hold on
plot([v11(1),v22(1)],[v11(2),v22(2)],'r');
title(['MRI image'])
subplot(1,3,2)
imagesc(CT_final_Sagittla)
axis equal
axis tight
hold on
plot([v11(1),v22(1)],[v11(2),v22(2)],'r');
title(['CT image in Sagittla Plan'])
subplot(1,3,3)
imagesc(1e3 * x_vec, 1e3 * y_vec, field_im1'/1e6);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title(['Transcranial Pressure Field at targrt (ACC) is ',num2str(amp_focus_coordinate_oral/1e6),'MPa'],'FontSize',16)
h=colorbar;
set(get(h,'title'),'string','MPa');
hold on
plot(bowl_pos(1)*1000,bowl_pos(2)*1000,'r*')
hold on
plot(focus_coordinate_oral(1)*1000,focus_coordinate_oral(2)*1000,'r*')
hold on
plot([bowl_pos(1)*1000,focus_coordinate_oral(1)*1000],[bowl_pos(2)*1000,focus_coordinate_oral(2)*1000],'r')
view(-90,-90)

%% Trouvez la mise au point d'origine et calculez la profondeur de mise au point
% %Exclure les points de valeur maximale sur le support

focus_coordinate_oral=[x_vec(targrt_position_new(1)),y_vec(targrt_position_new(2)),z_vec(targrt_position_new(3))];

% Contient la distribution du champ sonore ainsi que les points cibles et incidents et les médias
amp_focus_coordinate_oral=field_im1(targrt_position_new(1),targrt_position_new(2));%Trouver la pression acoustique du point cible d'origine（Pa）
figure;
imagesc(1e3 * x_vec, 1e3 * y_vec, field_im1'/1e6);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title(['Transcranial Pressure Field at targrt (ACC) is ',num2str(amp_focus_coordinate_oral/1e6),'MPa'],'FontSize',16)
h=colorbar;
set(get(h,'title'),'string','MPa');
hold on
plot(bowl_pos(1)*1000,bowl_pos(2)*1000,'r*')
hold on
plot(focus_coordinate_oral(1)*1000,focus_coordinate_oral(2)*1000,'r*')
hold on
plot([bowl_pos(1)*1000,focus_coordinate_oral(1)*1000],[bowl_pos(2)*1000,focus_coordinate_oral(2)*1000],'r')
view(-90,-90)



% %% Trouvez la ligne de contour de la pression acoustique focale -3dB
% field_im1_db=20*log10(field_im1/amp_max);
% db_number=-3;%Le nombre de dB à évaluer, tel que -3dB, -6dB
% %flipud Matrice retournée verticalement
% %fliplr Retournez la matrice à gauche et à droite
% figure()
% [M,c]=contour(flipud(field_im1_db'),[db_number db_number],'ShowText','on');%
% axis image;
% 
% M_exculude=M;
% M_exculude_index=find((M_exculude(1,:)<start_n));%Exclure les données en dehors du crâne
% M_exculude(:,M_exculude_index)=[];
% 
% figure
% plot(M_exculude(1,:),M_exculude(2,:))
% 
% 
% %% Effectuez un ajustement de courbe sur l'ellipse pour obtenir les axes majeur et mineur
% MM=M';
% p0=[0.005 0.005 0.005 0.005 0.005 0.005];%Valeur initiale ajustée
% F=@(p,x)p(1)*MM(:,1).^2+p(2)*MM(:,1).*MM(:,2)+p(3)*MM(:,2).^2+p(4)*MM(:,1)+p(5)*MM(:,2)+p(6);
% % Coefficients d'ajustement, méthode des moindres carrés
% p=nlinfit(MM,zeros(size(MM,1),1),F,p0);
% A=p(1)/p(6);
% B=p(2)/p(6);
% C=p(3)/p(6);
% D=p(4)/p(6);
% E=p(5)/p(6);
% %%Centre de l'ellipse
% X_center = (B*E-2*C*D)/(4*A*C - B^2);
% Y_center = (B*D-2*A*E)/(4*A*C - B^2);
% fprintf(' X_center=%g, Y_center=%g\n',X_center,Y_center);
% %%axe long et court
% a= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C+sqrt(((A-C)^2+B^2))));
% b= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C-sqrt(((A-C)^2+B^2))));

%% résultat de sortie
clc
fprintf('The wavelength is %g mm\n', wave_length);
fprintf('The number of points per wave is %g \n', ppw);
fprintf('The k-wave step is %g mm\n', kwave_step);
fprintf('The surface pressure  is %6.2f MPa\n',source_amp/1e6);
fprintf('The  pressure at focal point is %6.2f MPa\n',amp_max/1e6);
fprintf('The   pressure at targrt is %6.2f MPa\n',amp_focus_coordinate_oral/1e6);
fprintf('The focal length  is %6.2f mm\n',focal_length);
% fprintf(' Largeur axe long à mi-hauteur=%g mm\n',b*dx*1000);%Considérez l'unité de longueur de la grille, multipliée par la taille du pas
% fprintf(' Largeur du petit axe à mi-hauteur=%g mm\n',a*dx*1000);