function [] = plotS809(S)

%input data
x(1) = str2double(get(S.ed(1),'String'));
x(2) = str2double(get(S.ed(3),'String'));
x(3) = str2double(get(S.ed(5),'String'));
x(4) = str2double(get(S.ed(7),'String'));
x(5) = str2double(get(S.ed(9),'String'));

hParam = get(S.bg(1),'SelectedObject');
if(hParam == S.rd(1))
    param = 1;
elseif(hParam == S.rd(2))
    param = 2;
else
    param = 3;
end
        
hSelectBc = get(S.bg(2),'SelectedObject');
if(hSelectBc == S.rd(4))
    selectbc = 1;
else
    selectbc = 2;
end

np = str2double(get(S.ed(11),'String'));
tol = str2double(get(S.ed(12),'String'));
iter = str2double(get(S.ed(13),'String'));

x(2) = x(2)*(pi()/180);

myData = guidata(S.fh);

P = myData.pts;
n = length(P);

k = myData.k;
clamped = myData.clamped;

% c2 cubic interpolation of airfoil points

u = findknot(param,n,P);
[c,d] = c2cubic(selectbc,u,n,P);


% find the bezier segment that contains x(1), and the segment before it
for i = (n+1)/2:n
    if((P(i,1) <= x(1)) && (x(1) <= P(i+1,1)))
        index = i;
    end
end
for i = 1:4
    for j = 1:2
        ni = (3*index)-(3-i);
        part(i,j) = c(ni,j);
        bpart(i,j) = c(ni-3,j);
    end
end

val = bisection(x(1),part,tol,iter); % find parameter value that matches x(1)
[c1,unused] = split(val,4,part); % split found segment, keeping the first part
[dx,dy] = calcnewd(c1,bpart); % calculate new control point

% find the right place to split the d vector and insert the new nurbs 
% control points
bcurr = 1;
for i = 1:3:length(c)
    if((c(i,1) <= x(1)) && (x(1) <= c(i+3,1)))
        index = bcurr+1;
    end
    bcurr = bcurr + 1;
end

% find vector of first segment control points
for i = 1:index-1
    for j = 1:2
        d1(i,j) = d(i,j);
    end
end
d1(index,1) = dx;
d1(index,2) = dy;
d1(index+1,1) = c1(3,1);
d1(index+1,2) = c1(3,2);
d1(index+2,1) = c1(4,1);
d1(index+2,2) = c1(4,2);

% find bezier segment that contains x(5), and the segment after it
for i = (n+1)/2:n
    if((P(i,1) <= x(5)) && (x(5) <= P(i+1,1)))
        index = i;
    end
end
for i = 1:4
    for j = 1:2
        ni = (3*index)-(3-i);
        part(i,j) = c(ni,j);
        apart(i,j) = c(ni+3,j);
    end
end

val = bisection(x(5),part,tol,iter); % find parameter value that matches x(1)
[unused,c2] = split(val,4,part); % split found segment, keeping the first part
[dx,dy] = calcnewd(c2,apart); % calculate new control point

% find the right place to split the d vector and insert the new nurbs 
% control points
bcurr = 1;
for i = 1:3:length(c)
    if((c(i,1) <= x(5)) && (x(5) <= c(i+3,1)))
        index = bcurr+1;
    end
    bcurr = bcurr + 1;
end

% find vector of second segment control points
for i = index+1:length(d)
    for j = 1:2
        d2(i-(index)+3,j) = d(i,j);
    end
end

d2(1,1) = c2(1,1);
d2(1,2) = c2(1,2);
d2(2,1) = c2(2,1);
d2(2,2) = c2(2,2);
d2(3,1) = dx;
d2(3,2) = dy;

% determine cubic curve segment
cc(1,1) = x(1);
cc(1,2) = c1(4,2);
cc(4,1) = x(5);
cc(4,2) = c2(1,2);
dx1 = cc(4,1)-cc(1,1);
dx2 = cc(4,2)-cc(1,2);
mag = sqrt((dx1*dx1)+(dx2*dx2));
len = sqrt(((cc(4,1)-cc(1,1))^2)+((cc(4,2)-cc(1,2))^2));
cc(2,1) = cc(1,1) - (x(3)*(len/3)*(cos(x(2))));
cc(2,2) = cc(1,2) - (x(3)*(len/3)*(sin(x(2))));
cc(3,1) = cc(4,1) - (x(4)*(len/3)*(dx1/mag));
cc(3,2) = cc(4,2) - (x(4)*(len/3)*(dx2/mag));

% calculate the knot vectors for the nurbs curves
u1 = findnurbsknot(param,clamped,length(d1),k,d1);
ucc = [0;0;0;0;1;1;1;1];
u2 = findnurbsknot(param,clamped,length(d2),k,d2);

% save shape for later writing to file
myData.u1 = u1;
myData.ucc = ucc;
myData.u2 = u2;

myData.d1 = d1;
myData.cc = cc;
myData.d2 = d2;

guidata(S.fh,myData)

% parameterize airfoil curves based on number of points
for i = 1:np+1
	t(i) = (i-1)/np;
end

% find wieght vectors
for i = 1:length(d1)
    w1(i) = 1.0;
end
wcc = [1;1;1;1];
for i = 1:length(d2)
    w2(i) = 1.0;
end

% find points to plot
segg1 = deboor(length(d1),k,np,t,u1,w1,d1);
ccc = deboor(4,k,np,t,ucc,wcc,cc);
segg2 = deboor(length(d2),k,np,t,u2,w2,d2);

% calc step height and width
width = x(5)-x(1)
height = ccc(1,2)
for i = 2:length(ccc)
    if(ccc(i,2) < height)
        height = ccc(i,2);
    end
end
height
s1len = length(segg1(:,1))
height =  segg1(s1len,2) - height

set(S.tx(33),'String',height);
set(S.tx(35),'String',width);

% print graph to screen
cla(S.ax);
hold on
plot(S.ax,segg1(:,1),segg1(:,2))
plot(S.ax,ccc(:,1),ccc(:,2))
plot(S.ax,segg2(:,1),segg2(:,2))

val = get(S.ch(5),'Value');
if(val == 1)
    plot(S.ax,cc(:,1),cc(:,2),'o')
    plot(S.ax,cc(:,1),cc(:,2),'green')
end

hold off
axis square;
end