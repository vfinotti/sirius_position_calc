close all
clear all
clc

big_fonts = 1; % set fonts to big size;

%% Estimate real and measured beam positions

% Parameters

button_r = 3; %3
chamber_r = 12; %12

Kx = 8.5785; % mm;

x_array_length = 10; %length in mm
y_array_length = 10; %length in mm
array_size = 10;

% Create xy vector

xa = linspace(0, x_array_length, array_size); % axis coordinates
ya = zeros(1,array_size); % axis coordinates

xd = linspace(0, x_array_length/sqrt(2), array_size); % diagonal coordinates
yd = linspace(0, y_array_length/sqrt(2), array_size); % diagonal coordinates

xy = [xa xd; ya yd]';

%% Plot Matrix

% Create chamber plot

theta = linspace(0,2*pi); % Chamber draw

x_chamber = chamber_r*cos(theta);
y_chamber = chamber_r*sin(theta);

[x_button,y_button] = button_draw(chamber_r,button_r,4,pi/4);

figure
plot(xy(:,1),xy(:,2),'r*') % Plot data
hold on
plot(x_chamber,y_chamber,'k--') % Plot draws
for i=1:size(x_button,1)
    plot(x_button(i,:),y_button(i,:),'k.')
end
hold off
% axis([-chamber_r chamber_r -chamber_r chamber_r]*1.1)
axis equal
ll = legend('Chosen Points','location','Southeast');
tl = title('Chosen points for error analysis');
xl = xlabel('Real Beam Position (mm)');
yl = ylabel('Estimated Beam Position (mm)');

if big_fonts
    set(gca,'FontSize', 24);
    set(xl,'FontSize', 20);
    set(yl,'FontSize', 20);
    set(tl,'FontSize', 24);
end

grid on

print -depsc 1 % plotting figure

%% Convert to abcd coordinates
[abcd] = pos2abcd(xy,button_r,chamber_r);

%% Creating symbolic functions and Alfa

syms a b c d alfa;

% creating function
% f_ds    = '((a-c)-(b-d))/(a+b+c+d)';
f_ds    = '(a-c)/(a+c+alfa*(b+d))-(b-d)/((a+c)/alfa+b+d)';
f_pds   = '(a-c)/(a+c)+(d-b)/(d+b)';

da_ds=diff(f_ds,a);
db_ds=diff(f_ds,b);
dc_ds=diff(f_ds,c);
dd_ds=diff(f_ds,d);

da_pds=diff(f_pds,a);
db_pds=diff(f_pds,b);
dc_pds=diff(f_pds,c);
dd_pds=diff(f_pds,d);

% creating T gains and alfa
T = [100 100 50 50];
for i = 1:3
    alfav(i) = (T(2)+T(4))/(T(1)+T(3))*(10^(i-1));
end

%% Analisys

e = 1;

for j = 1:3
    for i= 1:size(abcd,1)

        da1 = da_ds;
        db1 = db_ds;
        dc1 = dc_ds;
        dd1 = dd_ds;
        
        da2 = da_pds;
        db2 = db_pds;
        dc2 = dc_pds;
        dd2 = dd_pds;
        
        % substitute values a
        da1=subs(da1,a,abcd(i));
        da1=subs(da1,b,abcd(i,2));
        da1=subs(da1,c,abcd(i,3));
        da1=subs(da1,d,abcd(i,4));
        da1=subs(da1,alfa,alfav(j));
        
        da2=subs(da2,a,abcd(i,1));
        da2=subs(da2,b,abcd(i,2));
        da2=subs(da2,c,abcd(i,3));
        da2=subs(da2,d,abcd(i,4));
        da2=subs(da2,alfa,alfav(j));
        
        % substitute values b
        db1=subs(db1,a,abcd(i,1));
        db1=subs(db1,b,abcd(i,2));
        db1=subs(db1,c,abcd(i,3));
        db1=subs(db1,d,abcd(i,4));
        db1=subs(db1,alfa,alfav(j));
        
        db2=subs(db2,a,abcd(i,1));
        db2=subs(db2,b,abcd(i,2));
        db2=subs(db2,c,abcd(i,3));
        db2=subs(db2,d,abcd(i,4));
        db2=subs(db2,alfa,alfav(j));
        
        % substitute values c
        dc1=subs(dc1,a,abcd(i,1));
        dc1=subs(dc1,b,abcd(i,2));
        dc1=subs(dc1,c,abcd(i,3));
        dc1=subs(dc1,d,abcd(i,4));
        dc1=subs(dc1,alfa,alfav(j));
        
        dc2=subs(dc2,a,abcd(i,1));
        dc2=subs(dc2,b,abcd(i,2));
        dc2=subs(dc2,c,abcd(i,3));
        dc2=subs(dc2,d,abcd(i,4));
        dc2=subs(dc2,alfa,alfav(j));
        % substitute values d
        dd1=subs(dd1,a,abcd(i,1));
        dd1=subs(dd1,b,abcd(i,2));
        dd1=subs(dd1,c,abcd(i,3));
        dd1=subs(dd1,d,abcd(i,4));
        dd1=subs(dd1,alfa,alfav(j));
        
        dd2=subs(dd2,a,abcd(i,1));
        dd2=subs(dd2,b,abcd(i,2));
        dd2=subs(dd2,c,abcd(i,3));
        dd2=subs(dd2,d,abcd(i,4));
        dd2=subs(dd2,alfa,alfav(j));
        
      
        % computing error
        e_ds(j,i)     = sqrt((da1)^2*(e)^2+(db1)^2*(e)^2+(dc1)^2*(e)^2+(dd1)^2*(e)^2);
        e_pds(j,i)    = sqrt((da2)^2*(e)^2+(db2)^2*(e)^2+(dc2)^2*(e)^2+(dd2)^2*(e)^2);
        
    end
end

% Separating Values

e_a_ds = e_ds(:,1:size(e_ds,2)/2);
e_d_ds = e_ds(:,size(e_ds,2)/2+1:end);

e_a_pds = e_pds(:,1:size(e_ds,2)/2);
e_d_pds = e_pds(:,size(e_ds,2)/2+1:end);

% axis plot
figure
plot(xa,e_a_ds(1,:),'r',xa,e_a_ds(2,:),'g',xa,e_a_ds(3,:),'b',xa,e_a_pds(1,:),'r--',xa,e_a_pds(2,:),'g--',xa,e_a_pds(3,:),'b--');
tl = title('Error Propatation - Axis');
yl = ylabel('Propagated Error');
xl = xlabel('Distance from axis [mm]');
legend(['DS - alfa = ' num2str(alfav(1))],['DS - alfa = ' num2str(alfav(2))],['DS - alfa = ' num2str(alfav(3))],['PDS - alfa = ' num2str(alfav(1))],['PDS - alfa = ' num2str(alfav(2))],['PDS - alfa = ' num2str(alfav(3))],'location','best')

if big_fonts
    set(gca,'FontSize', 24);
    set(xl,'FontSize', 20);
    set(yl,'FontSize', 20);
    set(tl,'FontSize', 24);
end

grid on

print -depsc 2 % plotting figure

% diagonal plot
figure
plot(xd,e_d_ds(1,:),'r',xd,e_d_ds(2,:),'g',xd,e_d_ds(3,:),'b',xd,e_d_pds(1,:),'r--',xd,e_d_pds(2,:),'g--',xd,e_d_pds(3,:),'b--');
tl = title('Error Propatation - Diagonal');
yl = ylabel('Propagated Error');
xl = xlabel('Distance from axis [mm]');
legend(['DS - alfa = ' num2str(alfav(1))],['DS - alfa = ' num2str(alfav(2))],['DS - alfa = ' num2str(alfav(3))],['PDS - alfa = ' num2str(alfav(1))],['PDS - alfa = ' num2str(alfav(2))],['PDS - alfa = ' num2str(alfav(3))],'location','best')

if big_fonts
    set(gca,'FontSize', 24);
    set(xl,'FontSize', 20);
    set(yl,'FontSize', 20);
    set(tl,'FontSize', 24);
end

grid on

print -depsc 3 % plotting figure