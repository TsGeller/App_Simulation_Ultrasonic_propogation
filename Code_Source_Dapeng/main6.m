%% 经颅超声声场仿真主程序，主要包括：靶点和入射点的选择，颅骨声学参数的获取，K_WAVE仿真
%%%版本：2022.9.21，修改人：李大鹏，使用程序：AB数据，猕猴E(已配准)

%% 初始化
clc;
clear;
close all;

%% 打开CT图像（CT已经配准）
[filename,filepath1]=uigetfile('E:\data\配准后数据\E\CT\.dcm'); 
%选择任意一个DIOCM文件即可，来取得CT图像的信息
filepath2=[filepath1 filename];
info1 = dicominfo(filepath2); %读取该DICOM数据的信息
[V1,spatial1,dim1] = dicomreadVolume(filepath1); 
%读取数据，V1为3维矩阵，代表CT原始数据（相对MRI配准）
%spatial:空间分辨率及坐标，dim:维度
V1 = squeeze(V1);%消除长度为1的维度
D1=single (V1);%转化为单精度
clear V1;
%%灰度转CT值
CT = D1.* info1.RescaleSlope + info1.RescaleIntercept;

%% 打开MRI图像
% filename2 就是你要用的DICOM文件的名字
[filename,filepath1]=uigetfile('E:\data\配准后数据\E\MRI\.dcm'); 
filepath2=[filepath1 filename];
info2 = dicominfo(filepath2);
[V2,spatial2,dim2] = dicomreadVolume(filepath1);
%V:原始数据（4D）,spatial:空间分辨率及坐标，dim:维度
V2 = squeeze(V2);%消除长度为1的维度
D2=single (V2);%转化为单精度
clear V2;
MRI=D2;


%% 获得CT和MRI图像的分辨率
MRI_space=[];%MRI的三维分辨率（mm）
MRI_space(1:2)=info2.PixelSpacing;
MRI_space(3)=info2.SliceThickness;

CT_space=MRI_space;%CT的三维分辨率（mm）,CT和MRI已配准，具有相同的空间分辨率

%% 设置颅骨CT数据的上下限及分层（水，软组织和猴子颅骨），不同数据阈值不同
CT_th1=0;%设置水的CT阈值为0 HU，低于0的全部置为水
CT_th2=1800;%设置颅骨的CT阈值下限为1800HU，0-1800HU的认为是软组织
CT_th3=7000;%设置颅骨的CT阈值上限为7000HU,1800-7000HU的认为是颅骨，大于7000的置为2400HU
CT(CT<CT_th1)=CT_th1;
CT(CT<CT_th2&CT>CT_th1)=1800;
CT(CT>CT_th3)=CT_th3;

%% 选择同一位置（水平面的中心）MRI和CT图像查看配准效果（可选的）
axial_plane_center1 = round(size(CT,3)/2);
axial_plane_center2 = round(size(MRI,3)/2);
CT_im=CT(:,:,axial_plane_center1);
MRI_im=MRI(:,:,axial_plane_center2);
CTandMRI=CT_im*0.3+MRI_im;%叠加图像
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


%% 计算超声的波长的k-wave的空间步长（3维仿真每个步长至少包含3个波长）
wave_speed=1500;%水中的波速 m/s
wave_frequency=0.5; %频率 MHz
wave_length=wave_speed/(wave_frequency*1e6)*1000;%波长 mm
ppw=5;% 每个网格点多少个波长（k-wave）
kwave_step=wave_length/ppw;%k-wave中每个网格点的步长



%% 对原始CT和MRI数据进行插值
[m,n,k]=size(CT);
[x,y,z] = meshgrid(1:n,1:m,1:k);
%%%注意：Vq = interp3(V,Xq,Yq,Zq) 假定一个默认的样本点网格。默认网格点覆盖区域 X=1:n、Y=1:m 和 Z=1:p，其中 [m,n,p] = size(V)
[xq,yq,zq] = meshgrid(1:kwave_step/CT_space(1):m,1:kwave_step/CT_space(2):n,1:kwave_step/CT_space(3):k);
CT_inter = interp3(x,y,z,CT,xq,yq,zq);%插值后的CT数据
MRI_inter = interp3(x,y,z,MRI,xq,yq,zq);%插值后的MRI数据

%% 设置颅骨CT数据的上下限及分层（水，软组织和猴子颅骨），插值后需要再次分层
CT_th1=0;%设置水的CT阈值为0 HU，低于0的全部置为水
CT_th2=1800;%设置颅骨的CT阈值下限为1800HU，0-1800HU的认为是软组织
CT_th3=7000;%设置颅骨的CT阈值上限为7000HU,1800-7000HU的认为是颅骨，大于7000的置为2400HU
CT_inter(CT_inter<CT_th1)=CT_th1;
CT_inter(CT_inter<CT_th2&CT_inter>CT_th1)=1800;
CT_inter(CT_inter>CT_th3)=CT_th3;

%% 对原始数据进行裁剪或者扩增，为换能器预留位置
%%需要引入换能器后才知道扩增大小
% %数据扩增
[m,n,k]=size(MRI_inter);%当前数据大小
add_number=30;%要增加的层数
%MRI数据扩增
MRI_add=zeros(m+add_number,n,k);
MRI_add(add_number+1:end,:,:)=MRI_inter;
%CT数据扩增
CT_add=zeros(m+add_number,n,k);
CT_add(add_number+1:end,:,:)=CT_inter;
%观察扩增后图像大小
% VolumeViewer3D(MRI_add)%工具箱
MRI_inter=[];
MRI_inter=MRI_add;
CT_inter=[];
CT_inter=CT_add;

%% 找到靶点所在平面，根据平面图像找到其余的坐标
VolumeViewer3D(MRI_inter);%需要工具箱，Coordinate 分别对应矩阵的行，列和页
%找到靶点所在的平面，确定是矢状面，横截面还是冠状面
target_plane_position=98;%找到ACC所在的平面

%% 画出靶点所在的平面
target_plane=squeeze(MRI_inter(:,:,target_plane_position));
figure
imagesc(target_plane)
axis equal
axis tight
colormap gray;

%% 确定靶点的其余坐标
ACC_target_xy=[120,120];%靶点
figure
imagesc(target_plane)
axis equal
axis tight
colormap gray;
hold on
plot((ACC_target_xy(1)),(ACC_target_xy(2)),'r+')


%% 计算颅骨的轮廓作为其切线方向
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
plot(n,m,'b')
hold on
plot(x_convhull,y_convhull,'r')
axis equal
axis tight
%% 确定换能器的位置，并且计算入射角度和路径
focol_distance=63.6;%聚焦换能器的真实焦距（从换能器顶点开始计算）
%%在矢状面中变换换能器，保证焦距
%%生成以靶点为圆心，焦距为半径的圆
theta = 1:0.5:360;
x0 = (ACC_target_xy(1));
y0 = (ACC_target_xy(2));
r = focol_distance/kwave_step;%半径占多少个网格
x_circle= x0 + r*cosd(theta);
% x_circle=round(x_circle);
y_circle = y0 + r*sind(theta);
% y_circle=round(y_circle);
%%保留在颅骨外部的点
in=[];%保存位置的特征，0表示在多边形外部，1表示在多边形边界和内部
for i=1:1:length(x_circle)
    in(i)=inpolygon(x_circle(i),y_circle(i),x_convhull,y_convhull);
end
x_circle1=x_circle(in==0);%取得内部的点
y_circle1=y_circle(in==0);%取得内部的点

%% 对待选择的点进行限制，排除不合理的点
x_circle1_cut=[];
y_circle1_cut=[];
i=find((x_circle1>50 & x_circle1<190));
x_circle1_cut= x_circle1(i);%取得内部的点
y_circle1_cut=y_circle1(i);%取得内部的点


%% 画图
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

%% 取得声束经颅的颅骨信息，x、y坐标以及HU值(修改)
pos_final_all={};%所有点穿过颅骨的坐标集合
skull_profile_all={};%所有点穿过颅骨的像素集合
angle_inside_all=[];%所有点穿过颅骨的角度集合
x2=x0;y2=y0;%靶点坐标
for i1=1:length(x_circle1_cut)
% x1=199;y1=121;
x1=x_circle1_cut(i1);%换能器顶点坐标
y1=y_circle1_cut(i1);

%%算出直线经过的图像像素坐标
pos1=[];
pos2=[];
pos3=[];
pos_final=[];
pos1 = [x1,y1;x2,y2];%取得直线的端点坐标
dx=pos1(2,1)-pos1(1,1);%第一列为x坐标
dz=pos1(2,2)-pos1(1,2);%第二列为y坐标
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
    %生成直线：(1)z=kx+b;(2)x=(z-b)/k;
%用x坐标计算
x11=min(pos1(2,1),pos1(1,1));x22=max(pos1(2,1),pos1(1,1));
    pos2(:,1)=x11:0.25:x22;
    pos2(:,2)=pos2(:,1).*k+b;
    pos2=round(pos2);
    pos2=unique(pos2,'rows','stable');
%用z坐标计算
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
%%根据直线经过的像素得到颅骨轮廓,找到颅骨进入点
skull_profile=[];
for i2=1:1:length(pos_final)
    skull_profile(i2)=target_plane_CT(pos_final(i2,2),pos_final(i2,1));
end
skull_profile_all{i1}=skull_profile;
%%计算一些颅骨信息（厚度，孔隙率，速度变化等等）
a=find(skull_profile~=0);%寻找颅骨中所有不为0的点
x_inside=pos_final(a(1),1);
y_inside=pos_final(a(1),2);%找到入射点的坐标

%%计算两个直线的夹角
%%找到入射点两侧对应的切线端点坐标
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

%%利用向量的点积公式计算角度 a.b=|a|*|b|*con(l);
v1 = [x1,y1] - [x2,y2];
v2 = [x4,y4] - [x3,y3];
angle_inside=acos(dot(v1,v2)/(norm(v1)*norm(v2)))*(180/pi);
%弧度转角度*（180/pi）
angle_inside_all(i1)=angle_inside;
end

%% 画出角度分布直方图
figure()
histogram(angle_inside_all)
xlabel('Angle range  ')
ylabel('Number ')
title('Potential ultrasound angle of incidence')
%% 找到与angel°最接近的入射角的顶点位置及其角度(x为所需要的角度，如90度)
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
% view(-90,90);%图像视野旋转90°

% %% 找到angel°入射角的经颅轮廓（斜入射不准确）
% skull_profile_1=skull_profile_all{k};
% skull_profile_x=0:kwave_step:(length(skull_profile_1)-1)*kwave_step;%mm
% skull_profile_2=skull_profile_1(skull_profile_1>1800);
% Skull_thickness=(length(skull_profile_2)-1)*kwave_step;
% figure()
% plot(skull_profile_x,skull_profile_1)
% xlabel('distance mm')
% ylabel('HU ')
% title(['Skull profile on ultrasound path,Skull thickness is ',num2str(Skull_thickness),' mm'])


%% 画一条直线连接靶点和入射点，确定旋转角度（使得声束垂直与图像，便于观察）
transducer_xy=[x_angel(1),y_angel(1)];%靶点
v1=ACC_target_xy;%靶点
v2=transducer_xy;%换能器顶点
v3=[];%垂直点
v3(1)=ACC_target_xy(1);v3(2)=transducer_xy(2);
distance_tran_tar=round(norm(v2-v1));%换能器顶点到靶点的距离（网格）
%计算超声声束方向和法线夹角
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

%% 根据包含靶点的平面图像确定靶点的坐标和超声入射点的坐标，此处为示例
transducer_position=[v2(2),v2(1),target_plane_position];
target_position=[v1(2),v1(1),target_plane_position];
transform_matrx=MRI_inter*0;%用来保存换能器的顶点和靶点
transform_matrx(target_position(1)-2:target_position(1)+2,target_position(2)-2:target_position(2)+2,target_position(3)-2:target_position(3)+2)=200;

% VolumeViewer3D(transform_matrx)

%% 根据计算的角度对数据进行旋转（此处绕着z轴旋转,volumeViewer显示靶区在xy切片）

MRI_rotate=imrotate3(MRI_inter,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);
CT_rotate=imrotate3(CT_inter,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);
transform_matrx_rotate=imrotate3(transform_matrx,-rotate_angle,[0 0 1],'linear','loose','FillValues',0);

%% 观察旋转后的数据
VolumeViewer3D(MRI_rotate);%需要工具箱，Coordinate 分别对应矩阵的行，列和页

%% 对数据进行裁剪(包括CT,MRI和位置矩阵)，减少运算量
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

%% 找到变换后的换能器顶点和靶点
Lax=find(transform_matrx_final==200);%
s=size(transform_matrx_final);%计算三维维数组的大小
[m,n,k]=ind2sub(s,Lax(round(length(Lax)/2)));%将最大值单下标转为三维多下标
targrt_position_new=[m,n,k];
transducer_position_new=[m-distance_tran_tar,n,k];


%% 在新的平面观察变换后的刺激靶点和换能器入射点的相对位置
v11=[targrt_position_new(2),targrt_position_new(1)];%靶点
v22=[transducer_position_new(2),transducer_position_new(1)];%换能器顶点
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

%% 计算出声束经过的颅骨厚度及CT值分布
Transcranial_HU=CT_final_Sagittla(v22(2):v11(2),v11(1));
position_skll=(0:abs(v11(2)-v22(2))).*kwave_step;

skull_index=find(Transcranial_HU>1800);%找到颅骨的截止点
skull_length=abs((skull_index(1)-skull_index(end)))*kwave_step;

figure()
plot(position_skll,Transcranial_HU)
xlabel('Ditance (mm)','FontSize',16);
ylabel('HU','FontSize',16);
title(['The thickness of the skull crossed by the ultrasound is ',num2str(skull_length),'mm'],'FontSize',16)


%% 设置颅骨CT数据的上下限及分层（水，软组织和颅骨),插值后需要再次分层
CT_th1=0;%设置水的CT阈值为0 HU，低于0的全部置为水
CT_th2=1800;%设置颅骨的CT阈值下限为150HU，0-150HU的认为是软组织
CT_th3=7000;%设置颅骨的CT阈值上限为2400HU,150-2400HU的认为是颅骨，大于2400的置为2400HU
CT_final(CT_final<CT_th1)=CT_th1;
CT_final(CT_final<CT_th2&CT_final>CT_th1)=1800;
CT_final(CT_final>CT_th3)=CT_th3;

%% 定义声学常数
HU_max=max(CT_final(:));
HU_water=0;
swater=wave_speed;%设置水中的声速度，对应于最小的HU
dwater=1000;%s水的密度
stissue=1560;%软组织的声速
dtissue=1030;%软组织的密度
sbone=3100;%假定为皮质骨中的声速度，对应于最大的HU
dbone=1900;%设置为皮质骨中的密度，对应于最大的HU
f=wave_frequency;%该衰减系数的测试超声频率（MHz）
b=1.1;%吸收系数a=a0*f^y;a0要化为[dB/(MHz^y cm)]
amin=0.2;%Attenuation  [dB/(MHz^y cm)]
amax=8;%%Attenuation  [dB/(MHz^y cm)]
atissue=0.6;%Attenuation  [dB/(MHz^y cm)]
%alpha_dB = 20*log10 （exp（alpha_Neper）） 
%alpha_dB= alpha_Neper * 20*log10 （exp） = 8.686 * alpha_Neper

%% 计算声场常数，包括声速、密度以及衰减
axial_plane_center = round(size(CT_final,3)/2);
% figure()
% imagesc(CT_final(:,:,axial_plane_center))
% axis equal
% axis tight
% colormap;
% h=colorbar;
% set(get(h,'title'),'string','HU');

%利用HU和密度及速度成线性关系求解
density=(CT_final-HU_water).*((dbone-dwater)/(HU_max-HU_water))+dwater;%计算密度
speed=(CT_final-HU_water).*((sbone-swater)/(HU_max-HU_water))+swater;%计算速度
k=(dbone-density)./(dbone-dwater);%计算孔隙率
atta=amin+(amax-amin).*(k.^0.5);%计算衰减系数
density(CT_final==1800)=dtissue;%软组织密度
speed(CT_final==1800)=stissue;%软组织速度
atta(CT_final==1800)=atissue;%软组织衰减系数
atta(CT_final==0)=amin;%设置水的衰减为amin

%% 显示声学参数的计算结果(可选的)
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


%% 设置k-wave仿真条件 %%

% 定义计算网格
[Nx,Ny,Nz]=size(CT_final);
dx = kwave_step/1e3;    	% x方向的步长
dy = dx ;            % y方向的步长
dz=dx ; % z方向的步长
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%定义介质的属性
medium.sound_speed = swater * ones(Nx, Ny, Nz);	% [m/s]
medium.density = dwater * ones(Nx, Ny, Nz);       % [kg/m^3]
medium.alpha_coeff = 0 * ones(Nx, Ny, Nz);  % [dB/(MHz^y cm)]
medium.alpha_power=b;

%设置衰减层
medium.sound_speed = speed;	% [m/s]
medium.density = density;       % [kg/m^3]
medium.alpha_coeff=atta; % [dB/(MHz^y cm)]

source_f0       = wave_frequency*1e6;      % source frequency [Hz]
%定义时间网格
% cfl             = 0.5;      % CFL 数目
cfl             = 0.2;      % CFL 数目
PPP = round(ppw / cfl);%每点多少个周期
dt = 1 / (PPP * source_f0);
%(1)自动定义时间
% kgrid.makeTime(medium.sound_speed);
%(2)手动定义时间
% cfl             = 0.5;      % CFL number
max_distance=sqrt(Nx^2+Ny^2+Nz^2)*dx;%超声传播的最大距离
t_end           = max_distance/1500;    % 总的计算时间(必须大于系统的稳定时间)
% PPP = round(ppw / cfl);%计算每个时间周期多少个点
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

%定义信号源
% 定义换能器的性质
source_f0       = wave_frequency*1e6;      % source frequency [Hz]
source_roc      = 66e-3;    % 换能器的曲率半径 [m]
source_diameter = 64e-3;    % 换能器的有效直径 [m]
source_amp      = 1e6;      % 换能器的表面声压 [Pa]

%定义换能器的类型，此处为聚焦的碗状结构
% 定义碗的空间位置和旋转
bowl_pos = [kgrid.x_vec(transducer_position_new(1)), kgrid.y_vec(transducer_position_new(2)), kgrid.z_vec(transducer_position_new(3))];
focus_pos = [kgrid.x_vec(targrt_position_new(1)), kgrid.y_vec(targrt_position_new(2)), kgrid.z_vec(targrt_position_new(3))];

% 定义一个空的 kWaveArray
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 10);
% 增加碗的形状的换能器
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);

%% 定义 binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);
% 定义信号时间序列
source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);
% 定义信号源，信号重新分布，这一步需要长时间和大的内存
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
% VolumeViewer3D(int8(source.p_mask));
% 
% VolumeViewer3D(medium.sound_speed);
% VolumeViewer3D(medium.density);
% VolumeViewer3D(medium.alpha_coeff);
% 定义传感器
sensor.mask = [1,1,1, Nx, Ny,Nz].';%要记录的网格为整个网格 
sensor.record = {'p_max'};%要记录的参数为声压的最大值
% sensor.record = {'p_min'};%要记录的参数为声压的最小值
% 只记录稳定时候的值
record_periods  = 3;        % 要记录的周期数目
sensor.record_start_index = kgrid.Nt - record_periods * PPP + 1;

%% 定义仿真输出条件
% display_mask = source.p_mask;%显示换能器网格
display_mask = CT_final*0;%显示颅骨
display_mask(CT_final>1800)=1;
% VolumeViewer3D(display_mask);
% input_args = {'DisplayMask', display_mask, 'PlotLayout', true, 'PMLInside', false, 'PlotPML', false,...
%     'DataCast', 'single', 'DataRecast', true,'PlotScale', ...
%     [-1, 1] * source_amp};%cpu运算
input_args = {'DisplayMask', display_mask, 'RecordMovie',true,'MovieName','1','PlotLayout', true,  'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10},...
    'PMLInside', false, 'PlotPML', false,...
    'DataCast', 'gpuArray-single', 'DataRecast', true,'PlotScale', ...
    [-1, 1] * source_amp};%gpu运算,并且记录视频
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});%开始3维k-wave仿真

%% 仿真结果评估 %%
%% 从传感器器中提取信号幅度
% amp = extractAmpPhase(sensor_data.p_max, 1/kgrid.dt, source_f0, ...
%     'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);%第二维表示时间
% % 重构矩阵
% amp = reshape(amp, Nx, Ny,Nz);
 amp=sensor_data.p_max;


%% 获取3轴坐标
x_vec = kgrid.x_vec;
y_vec = kgrid.y_vec;
z_vec = kgrid.z_vec;

%% 找到焦点，计算焦点深度
% %排除介质上的极大值点
Estimated_focal_length=10/1000;%目测的焦点距离
start_n=round(Estimated_focal_length/dx);
amp_exculude=amp;
amp_exculude(1:start_n,:,:)=0;

amp_max=max(amp_exculude(:));%计算三维维数组的最大值,Pa
s=size(amp);%计算三维维数组的大小
Lax=find(amp==amp_max);%计算最大值位置的单下标
[m,n,k]=ind2sub(s,Lax);%将最大值单下标转为三维多下标
Loc_focus=[m,n,k];%最大值位置下标
focus_coordinate=[x_vec(m),y_vec(n),z_vec(k)];
focal_length=norm(focus_coordinate-bowl_pos)*1000;%算出焦点距离换能器顶点距离,mm

%% 画出包含换能器顶点和靶点的图像
field_im1=amp(:,:,targrt_position_new(3));


% 只包含声场分布的图
figure;
imagesc(1e3 * x_vec, 1e3 * y_vec, field_im1'/1e6);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Pressure Field');
h=colorbar;
set(get(h,'title'),'string','MPa');
%结果逆时针旋转90度
view(-90,-90)

% 包含声场分布以及靶点和入射点以及介质
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


%% 画出CT、MRI与声场的比较图像
amp_focus_coordinate_oral=field_im1(targrt_position_new(1),targrt_position_new(2));%求出原来靶点的声压（Pa）
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

%% 找到原焦点，计算焦点深度
% %排除介质上的极大值点

focus_coordinate_oral=[x_vec(targrt_position_new(1)),y_vec(targrt_position_new(2)),z_vec(targrt_position_new(3))];

% 包含声场分布以及靶点和入射点以及介质
amp_focus_coordinate_oral=field_im1(targrt_position_new(1),targrt_position_new(2));%求出原来靶点的声压（Pa）
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



% %% 求焦点声压-3dB的等高线
% field_im1_db=20*log10(field_im1/amp_max);
% db_number=-3;%要评估的dB数，如-3dB，-6dB
% %flipud 矩阵垂直翻转
% %fliplr 矩阵左右翻转
% figure()
% [M,c]=contour(flipud(field_im1_db'),[db_number db_number],'ShowText','on');%
% axis image;
% 
% M_exculude=M;
% M_exculude_index=find((M_exculude(1,:)<start_n));%排除颅骨外的数据
% M_exculude(:,M_exculude_index)=[];
% 
% figure
% plot(M_exculude(1,:),M_exculude(2,:))
% 
% 
% %% 对椭圆做曲线拟合，得到长短轴
% MM=M';
% p0=[0.005 0.005 0.005 0.005 0.005 0.005];%拟合的初值
% F=@(p,x)p(1)*MM(:,1).^2+p(2)*MM(:,1).*MM(:,2)+p(3)*MM(:,2).^2+p(4)*MM(:,1)+p(5)*MM(:,2)+p(6);
% % 拟合系数，最小二乘方法
% p=nlinfit(MM,zeros(size(MM,1),1),F,p0);
% A=p(1)/p(6);
% B=p(2)/p(6);
% C=p(3)/p(6);
% D=p(4)/p(6);
% E=p(5)/p(6);
% %%椭圆中心
% X_center = (B*E-2*C*D)/(4*A*C - B^2);
% Y_center = (B*D-2*A*E)/(4*A*C - B^2);
% fprintf(' X_center=%g, Y_center=%g\n',X_center,Y_center);
% %%长短轴
% a= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C+sqrt(((A-C)^2+B^2))));
% b= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C-sqrt(((A-C)^2+B^2))));

%% 输出结果
clc
fprintf('The wavelength is %g mm\n', wave_length);
fprintf('The number of points per wave is %g \n', ppw);
fprintf('The k-wave step is %g mm\n', kwave_step);
fprintf('The surface pressure  is %6.2f MPa\n',source_amp/1e6);
fprintf('The  pressure at focal point is %6.2f MPa\n',amp_max/1e6);
fprintf('The   pressure at targrt is %6.2f MPa\n',amp_focus_coordinate_oral/1e6);
fprintf('The focal length  is %6.2f mm\n',focal_length);
% fprintf(' 长轴半高宽=%g mm\n',b*dx*1000);%考虑网格的单位长度，乘以步长
% fprintf(' 短轴半高宽=%g mm\n',a*dx*1000);