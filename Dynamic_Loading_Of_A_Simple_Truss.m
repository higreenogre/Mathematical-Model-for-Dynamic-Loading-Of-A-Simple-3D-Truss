%                       Submitted By George Henri T20CE011 

% 3D truss program Mini Project CE 555                                                           

clc
clear all % clear memory
tic       % starts a stopwatch timer
fig=figure;


scale=500;
% elements nodes
ele_nod=[1 5;2 5;3 5;4 5];

%number of elements
num_ele=size(ele_nod,1);


%number of nodes
num_nod=max(max(ele_nod));


% nodes coordinates
nod_coor=[-6 0 8;12 0 8;6 0 -8;-12 0 -8;0 24 0]*12;

%% elements degree of freedom (DOF) 
for i=1:num_ele
    for j= 1:2
        ele_dof(i,3*j-2)=ele_nod(i,j)*3-2;
        ele_dof(i,3*j-1)=ele_nod(i,j)*3-1;
        ele_dof(i,3*j)=ele_nod(i,j)*3;
    end
end
 

% A, E, L are cross sectional area, Young's modulus, length of elements,respectively.

A= ones(num_ele,1)*8.4;

E= ones(num_ele,1)*10000;

% initial zero matrix for all matrices
displacement=zeros(3*num_nod,1);

stiffness=zeros(3*num_nod);
%% 
force={};
cyclelength=100;
cycle=linspace(1,6*pi,cyclelength);
for e=1:3*num_nod
    force{e}=zeros(1,cyclelength);
end
static=ones(1,cyclelength);
%% applied loads at DOFs

force{13}= sin(cycle+1)*50;
force{14}= -10*sin(cycle-2);
force{15}= sin(cycle)*50;

p=plot3(0,0,0);

daspect([1,1,1])

hold on


for times=1:cyclelength

%Boundary conditions

displacement (1,1)=0.0;
displacement (2,1)=0.0;
displacement (3,1)=0.0;
displacement (4,1)=0.0;



%% computation of the system stiffness matrix
for e=1:num_ele
    member = strcat("Member ",num2str(e));
    
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2+...
      (nod_coor(ele_nod(e,2),3)-nod_coor(ele_nod(e,1),3))^2);
  
 Cx=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 Cy=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 Cz=(nod_coor(ele_nod(e,2),3)-nod_coor(ele_nod(e,1),3))/L(e);
 
 T=[Cx Cy Cz 0 0 0;0 0 0 Cx Cy Cz];
 t(:,:,e)=T;
 
    localk=A(e)*E(e)/L(e)*[1 -1;-1 1];
 
 k=(T'*localk*T);
   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);
    for i=1:6
        for j=1:6
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);
        end
    end
end

a=[13,14,15]; %input
targetstiffness = stiffness(a,a);

%% known force array
known_f_a=[13,14,15]';
for i=1:size(known_f_a,1)
   dis_new(i,1)=displacement(known_f_a(i,1),1);
   forcecycle=force{known_f_a(i,1)};
   force_new(i,1)=forcecycle(times);
end

for i=1:size(known_f_a,1)
    for j=1:size(known_f_a,1)
    stiff_new(i,j)=stiffness(known_f_a(i,1),known_f_a(j,1));
end
end

% solving the partitioned matrix 
dis_new=stiff_new\force_new;

for i=1:size(known_f_a,1)
  displacement(known_f_a(i,1),1)=dis_new(i,1);
end

% known dicplacement array

known_dis_a=[1:12]';

for i=1:size(known_dis_a,1)
    forcecycle=force{known_dis_a(i,1)};
    forcecycle(times)=stiffness(known_dis_a(i,1),:)*displacement;
    force{known_dis_a(i,1)}=forcecycle;
end
targetdisplacement =displacement(a);


% stress in elements
for e=1:num_ele
    member = strcat("Member ",num2str(e));
    T=t(:,:,e);
    localk=A(e)*E(e)/L(e)*[1 -1;-1 1];
    
 member_displacement=T*displacement(ele_dof(e,:));
 member_localforces=localk*member_displacement;
 stress(e)=member_localforces(2)/A(e); %Calculates Axial Stress of each members
 member_globalforces=T'*member_localforces; %transforms forces to global
 
end
k=0;
    for i=1:num_nod
        for j=1:3
            k=k+1;
            nod_coor_def(times,i,j)=nod_coor(i,j)+scale*displacement(k,1);%Calculates Deformed Displacements to scale
        end
    end


stresscycle(times,:)=stress;
end
%% 
for times=1:cyclelength % Solves for Each Frame
    clf
    maxstress=max(stresscycle(times,:));minstress=min(stresscycle(times,:));

for e=1:num_ele
    red(e)=(stresscycle(times,e)-minstress)/(maxstress-minstress); % Interpolates the red value of RGB for each member
    blue(e)=(stresscycle(times,e)-maxstress)/(minstress-maxstress);% Interpolates the red value of RGB for each member
end
%% 
for e=1:num_ele
    x=[nod_coor_def(times,ele_nod(e,1),1) nod_coor_def(times,ele_nod(e,2),1)];
    y=[nod_coor_def(times,ele_nod(e,1),2) nod_coor_def(times,ele_nod(e,2),2)];
    z=[nod_coor_def(times,ele_nod(e,1),3) nod_coor_def(times,ele_nod(e,2),3)];
    
    p=plot3(x,z,y,'LineWidth',3);
    legend( strcat ('stress 1 = ', num2str (round(stresscycle(times,1)*1000)/1000)), strcat ('stress 2 = ', num2str (round(stresscycle(times,2)*1000)/1000)), strcat ('stress 3 = ', num2str (round(stresscycle(times,3)*1000)/1000)), strcat ('stress 4 = ', num2str (round(stresscycle(times,4)*1000)/1000)))

    p.Color=[red(e) 0 blue(e)];
    xlabel('X-Axis')
    ylabel('Z-Axis')
    zlabel('Y-Axis')
    
    hold on
    
    daspect([1,1,1])
    
    
end
title("Stress Deformation")
%% 
for node=1:num_nod
    x=nod_coor_def(times,node,1);
    y=nod_coor_def(times,node,2);
    z=nod_coor_def(times,node,3);
    text(x,z,y,num2str(node));
end
%% 
movieVector(times) = getframe(fig,[5 5 530 400]); % Captures the frames from the figure
end 
% Writes the  Video File
mywriter = VideoWriter('Stress_Deformation', 'MPEG-4');
    mywriter.FrameRate = 20;
    
    open(mywriter);
    writeVideo(mywriter, movieVector);
    close(mywriter);
toc 