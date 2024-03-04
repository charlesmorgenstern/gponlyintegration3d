using LinearAlgebra
using Printf
using Delaunay
using Plots
#################################################################
#################################################################
#Original code and recurrence relation
#################################################################
#################################################################
function getpaths(n,h,r)

f(x,y,z) = 3.0 + cos(x) + cos(y) + cos(z)
g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

points=points_on_sphere_regular(n,[0.0;0.0;0.0],r)
points=reduce(hcat,points)'
points=copy(points)

np=length(points[:,1])

gps=[]
remove=[]

 for i=1:np

  gp=[0.0 0.0 0.0; points[i,:]']
  check=true

   while check

    xn=gp[end,1]
    yn=gp[end,2]
    zn=gp[end,3]
    grad=-g(xn,yn,zn)
    grad=grad./norm(grad)
    kx1=grad[1]
    ky1=grad[2]
    kz1=grad[3]
    grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
    grad=grad./norm(grad)
    kx2=grad[1]
    ky2=grad[2]
    kz2=grad[3]
    grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
    grad=grad./norm(grad)
    kx3=grad[1]
    ky3=grad[2]
    kz3=grad[3]
    grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
    grad=grad./norm(grad)
    kx4=grad[1]
    ky4=grad[2]
    kz4=grad[3]


    xx=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
    yy=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
    zz=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
    gp=[gp;xx yy zz]

    temp=abs.([xx yy zz])
    hd=h/2
     if norm(temp-[0 0 pi])<hd || norm(temp-[0 pi 0])<hd || norm(temp-[pi 0 0])<hd || norm(temp-[0 pi pi])<hd || norm(temp-[pi 0 pi])<hd || norm(temp-[pi pi 0])<hd
      push!(remove,i)
      check=false
     elseif norm(temp-[pi pi pi])<1.5*h
      v=[xx yy zz]./temp
      v=[pi pi pi].*v
      gp=[gp;v]
      check=false
      push!(gps,gp)
     end
   end
 end

remove=reverse(remove,dims=1)
nr=length(remove)
 
 for i=1:nr
  points=[points[1:remove[i]-1,:];points[remove[i]+1:end,:]]
 end

mesh=delaunay(points)

areas=[spherical_triangle_area(mesh.points[v, :],[0.0;0.0;0.0]) for v in eachrow(mesh.convex_hull)]

surface_area=4.0*pi*r^2

np=length(points[:,1])
nf=length(mesh.convex_hull[:,1])

area=Vector{Float64}(undef,np)

 for i=1:np
  area[i]=0.0
   for j=1:nf
    if i==mesh.convex_hull[j,1] || i==mesh.convex_hull[j,2] || i==mesh.convex_hull[j,3]

     area[i]=area[i]+areas[j]/3

    end
   end
 end

volume=4/3*pi*r^3

mesh=Dict("points"=>mesh.points,
          "faces"=>mesh.convex_hull,
          "areas"=>areas,
          "area"=>area,
          "surface_area"=>surface_area)

sys=Dict("mesh"=>mesh,
         "f"=>f,
         "g"=>g,
         "gps"=>gps,
         "sphere_volume"=>volume)

return sys
end
#################################################################
#################################################################
function points_on_sphere_regular(N::Int, center::Vector{Float64}, radius::Float64)
    points = []
    Ncount = 0

    a = 4 * pi / N
    d = sqrt(a)
    M_theta = round(pi / d)
    d_theta = pi / M_theta

    for m in 0:(M_theta - 1)
        theta = pi * (m + 0.5) / M_theta
        M_phi = round(2 * pi * sin(theta) / d_theta)

        for n in 0:(M_phi - 1)
            phi = 2 * n * pi / M_phi
            x = radius * sin(theta) * cos(phi) + center[1]
            y = radius * sin(theta) * sin(phi) + center[2]
            z = radius * cos(theta) + center[3]
            push!(points, [x; y; z])
            Ncount += 1
        end
    end

    return points
end
#################################################################
#################################################################
function plotmesh(sys)

p=sys["mesh"]["points"]
f=sys["mesh"]["faces"]
n=length(f[:,1])

fig=Plots.plot()

 for i=1:n
  x=p[f[i,:],1]
  y=p[f[i,:],2]
  z=p[f[i,:],3]
  Plots.plot!([x; x[1]],[y; y[1]],[z; z[1]],legend=false)
 end

return fig
end
#################################################################
#################################################################
function spherical_triangle_area(vertices::Matrix{Float64}, o::Vector{Float64})
    a = zeros(3)

    for i = 1:3
        j = mod1(i + 1, 3)
        k = mod1(i + 2, 3)
        po = vertices[i, :] .- o

        v1 = vertices[j, :] .- vertices[i, :]
        vtmp = cross(po, v1)
        v1 = normalize(cross(vtmp, po))

        v2 = vertices[k, :] .- vertices[i, :]
        vtmp = cross(po, v2)
        v2 = normalize(cross(vtmp, po))

        a[i] = acos(clip(dot(v1, v2), -1.0, 1.0))  # Vector angle
    end

    r = norm(vertices[1, :] .- o)  # Distance
    e = sum(a) - π
    area = r^2 * e

    return area
end
#################################################################
#################################################################
function clip(v, l, h)
    if v < l
        return l
    elseif v > h
        return h
    else
        return v
    end
end
#################################################################
#################################################################
function plotgps(sys)

fig=Plots.plot()

ngp=length(sys["gps"])

 for i=1:ngp
  Plots.plot!(sys["gps"][i][:,1],sys["gps"][i][:,2],sys["gps"][i][:,3],legend=false)
 end


return fig
end
#################################################################
#################################################################
function curv(r)

x=r[1]
y=r[2]
z=r[3]

f(x,y,z) = 3.0 + cos(x) + cos(y) + cos(z)
g(x,y,z) = [-sin(x), -sin(y), -sin(z)]
h(x,y,z)=[-cos(x) 0 0; 0 -cos(y) 0;0 0 -cos(z)]

g1=g(x,y,z)
gn=norm(g1)
N=g1./gn

P=[1 0 0; 0 1 0; 0 0 1]-N*N'

H=h(x,y,z)

S=P*(H./gn)

E=eigen(S)

T=abs.(E.values)

d1=findmax(T)
d1=d1[2]
k1=E.values[d1]
v1=E.vectors[:,d1]
T[d1]=0

d2=findmax(T)
d2=d2[2]
k2=E.values[d2]
v2=E.vectors[:,d2]

vals=real(-[k1,k2])
vectors=[]
push!(vectors,real(v1))
push!(vectors,real(v2))

out=Dict("eigen_values"=>vals,
         "eigen_vectors"=>vectors)

return out
end
#################################################################
#################################################################
function getcurvs(sys)

ngp=length(sys["gps"])

k=[]

 for i=1:ngp
  kgp=[]
  np=length(sys["gps"][i][:,1])
  temp=curv(sys["gps"][i][2,:])
  push!(kgp,temp)
   for j=3:np

    w1=temp["eigen_vectors"][1]
    w2=temp["eigen_vectors"][2]

    temp=curv(sys["gps"][i][j,:])

    v1=temp["eigen_vectors"][1]
    v2=temp["eigen_vectors"][2]

    e1=temp["eigen_values"][1]
    e2=temp["eigen_values"][2]

     if norm(cross(w1,v1))>norm(cross(w1,v2))
      temp["eigen_vectors"][1]=v2
      temp["eigen_vectors"][2]=v1
      temp["eigen_values"][1]=e2
      temp["eigen_values"][2]=e1
     end

    v1=temp["eigen_vectors"][1]
    v2=temp["eigen_vectors"][2]


     if dot(v1,w1)<0 #norm(v1-w1)>norm(-v1-w1)
      temp["eigen_vectors"][1]=-v1
     end

     if dot(v2,w2)<0 #norm(v2-w2)>norm(-v2-w2)
      temp["eigen_vectors"][2]=-v2
     end

    push!(kgp,temp)
   end
  push!(k,kgp)
 end


sys["k"]=k
return sys
end
#################################################################
#################################################################
function getgpsegmentlengths(sys)

ngp=length(sys["gps"])
dgp=[]

 for i=1:ngp
  d=[]
  np=length(sys["gps"][i][:,1])
   for j=1:np-1
    temp=norm(sys["gps"][i][j,:]-sys["gps"][i][j+1,:])
    push!(d,temp)
   end
  push!(dgp,d)
 end

sys["dgp"]=dgp
return sys
end
#################################################################
#################################################################
function getarclengths(sys)

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]

 for i=1:ngp

  a1=[]
  a2=[]
  sa=[]
  vol=[]

  start_area=sys["mesh"]["area"][i]
  push!(a1,sqrt(start_area))
  push!(a2,sqrt(start_area))
  np=length(sys["gps"][i][:,1])

   for j=3:np

    k1=sys["k"][i][j-2]["eigen_values"][1]
    k2=sys["k"][i][j-2]["eigen_values"][2]

    ds=sys["dgp"][i][j-1]

    temp1=a1[end]*(1+k1*ds)
    temp2=a2[end]*(1+k2*ds)

    push!(sa,temp1*temp2)
    push!(vol,temp1*temp2*ds)

    push!(a1,temp1)
    push!(a2,temp2)

   end

  push!(arc1,a1)
  push!(arc2,a2)
  push!(surfarea,sa)
  push!(volumes,vol)

 end

sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes

return sys
end
#################################################################
#################################################################
function getvolume(sys)

ngp=length(sys["volumes"])
volume=0.0

 for i=1:ngp
  volume=volume+sum(sys["volumes"][i])
 end

volume=volume+sys["sphere_volume"]

sys["volume"]=volume
return sys
end
#################################################################
#################################################################
function runvolume(n,h,r)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengths(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolume(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolume(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
#First correction term code, fdm, normal point sampling
#################################################################
#################################################################
function getarclengthsc1v1(sys)

h=0.00001

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]


for i=1:ngp
 a1=[]
 a2=[]
 sa=[]
 vol=[]
 start_area=sys["mesh"]["area"][i]
 push!(a1,sqrt(start_area/pi))
 push!(a2,sqrt(start_area/pi))
 np=length(sys["gps"][i][:,1])
  for j=3:np

   k1=sys["k"][i][j-2]["eigen_values"][1]
   k2=sys["k"][i][j-2]["eigen_values"][2]
   ds=sys["dgp"][i][j-1]

   d1=sys["k"][i][j-2]["eigen_vectors"][1]
   d2=sys["k"][i][j-2]["eigen_vectors"][2]
   start=start=sys["gps"][i][j-1,:]
   p1=start+.5*a1[end]*d1
   p2=start-.5*a1[end]*d1
   p3=start+.5*a2[end]*d2
   p4=start-.5*a2[end]*d2

   temp1=curv(p1)
    check1=abs(temp1["eigen_values"][1]-k1)
    check2=abs(temp1["eigen_values"][2]-k1)
    if check1<check2
    k11=temp1["eigen_values"][1]
    else
    k11=temp1["eigen_values"][2]
    end

   temp2=curv(p2)
    check1=abs(temp2["eigen_values"][1]-k1)
    check2=abs(temp2["eigen_values"][2]-k1)
    if check1<check2
    k12=temp2["eigen_values"][1]
    else
    k12=temp2["eigen_values"][2]
    end

   temp3=curv(p3)
    check1=abs(temp3["eigen_values"][1]-k2)
    check2=abs(temp3["eigen_values"][2]-k2)
    if check1<check2
    k23=temp3["eigen_values"][1]
    else
    k23=temp3["eigen_values"][2]
    end

   temp4=curv(p4)
    check1=abs(temp4["eigen_values"][1]-k2)
    check2=abs(temp4["eigen_values"][2]-k2)
    if check1<check2
    k24=temp4["eigen_values"][1]
    else
    k24=temp4["eigen_values"][2]
    end

    dist1=.5*a1[end]
    dist2=.5*a1[end]
    dist3=.5*a2[end]
    dist4=.5*a2[end]


   dk1=(k11-k1)/dist1
   dk2=(k12-k1)/dist2
   dk3=(k23-k2)/dist3
   dk4=(k24-k2)/dist4


   temp1=a1[end]*(1+k1*ds+dk1/2*a1[end]*ds+dk2/2*a1[end]*ds)
   temp2=a2[end]*(1+k2*ds+dk3/2*a2[end]*ds+dk4/2*a2[end]*ds)

   push!(sa,pi*temp1*temp2)
   push!(vol,pi*temp1*temp2*ds)

   push!(a1,temp1)
   push!(a2,temp2)
  end
 push!(arc1,a1)
 push!(arc2,a2)
 push!(surfarea,sa)
 push!(volumes,vol)
end


sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes


return sys
end
#################################################################
#################################################################
function runvolumec1v1(n,h,r)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengthsc1v1(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumec1v1(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumec1v1(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
#First correction term code, fdm, curv based point sampling
#################################################################
#################################################################
function getarclengthsc1v2(sys)

g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

h=0.00001

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]


for i=1:ngp
 a1=[]
 a2=[]
 sa=[]
 vol=[]
 start_area=sys["mesh"]["area"][i]
 push!(a1,sqrt(start_area/pi))
 push!(a2,sqrt(start_area/pi))
 np=length(sys["gps"][i][:,1])
  for j=3:np

   k1=sys["k"][i][j-2]["eigen_values"][1]
   k2=sys["k"][i][j-2]["eigen_values"][2]
   ds=sys["dgp"][i][j-1]

   d1=sys["k"][i][j-2]["eigen_vectors"][1]
   d2=sys["k"][i][j-2]["eigen_vectors"][2]
   start=start=sys["gps"][i][j-1,:]

   r1=1/abs(k1)
   r2=1/abs(k2)

   grad=-g(start[1],start[2],start[3])
   grad=grad./norm(grad)

   if k1>0
   c1=start-r1.*grad
   v1=-grad
   else
   c1=start+r1.*grad
   v1=grad
   end

   if k2>0
   c2=start-r2.*grad
   v2=-grad
   else
   c2=start+r2.*grad
   v2=grad
   end

   f1(t)=c1+r1*cos(t)*d1+r1*sin(t)*v1
   f2(t)=c2+r2*cos(t)*d2+r2*sin(t)*v2

   t1=1.5*pi
   t2=1.5*pi

   dt1=.5*a1[end]/r1
   dt2=.5*a2[end]/r2

   p1=f1(t1+dt1)
   p2=f1(t1-dt1)
   p3=f2(t2+dt2)
   p4=f2(t2-dt2)

   temp1=curv(p1)
    check1=abs(temp1["eigen_values"][1]-k1)
    check2=abs(temp1["eigen_values"][2]-k1)
    if check1<check2
    k11=temp1["eigen_values"][1]
    else
    k11=temp1["eigen_values"][2]
    end

   temp2=curv(p2)
    check1=abs(temp2["eigen_values"][1]-k1)
    check2=abs(temp2["eigen_values"][2]-k1)
    if check1<check2
    k12=temp2["eigen_values"][1]
    else
    k12=temp2["eigen_values"][2]
    end

   temp3=curv(p3)
    check1=abs(temp3["eigen_values"][1]-k2)
    check2=abs(temp3["eigen_values"][2]-k2)
    if check1<check2
    k23=temp3["eigen_values"][1]
    else
    k23=temp3["eigen_values"][2]
    end

   temp4=curv(p4)
    check1=abs(temp4["eigen_values"][1]-k2)
    check2=abs(temp4["eigen_values"][2]-k2)
    if check1<check2
    k24=temp4["eigen_values"][1]
    else
    k24=temp4["eigen_values"][2]
    end


    dist1=.5*a1[end]
    dist2=.5*a1[end]
    dist3=.5*a2[end]
    dist4=.5*a2[end]

   dk1=(k11-k1)/dist1
   dk2=(k12-k1)/dist2
   dk3=(k23-k2)/dist3
   dk4=(k24-k2)/dist4


   temp1=a1[end]*(1+k1*ds+dk1/2*a1[end]*ds+dk2/2*a1[end]*ds)
   temp2=a2[end]*(1+k2*ds+dk3/2*a2[end]*ds+dk4/2*a2[end]*ds)

   push!(sa,pi*temp1*temp2)
   push!(vol,pi*temp1*temp2*ds)

   push!(a1,temp1)
   push!(a2,temp2)
  end
 push!(arc1,a1)
 push!(arc2,a2)
 push!(surfarea,sa)
 push!(volumes,vol)
end


sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes


return sys
end
#################################################################
#################################################################
function runvolumec1v2(n,h,r)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengthsc1v2(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumec1v2(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumec1v2(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
# Original reccurence relation, all paths included
#################################################################
#################################################################
function getallpaths(n,h,r)

f(x,y,z) = 3.0 + cos(x) + cos(y) + cos(z)
g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

points=points_on_sphere_regular(n,[0.0;0.0;0.0],r)
points=reduce(hcat,points)'
points=copy(points)

np=length(points[:,1])
gps=[]
remove=[]
for i=1:np
gp=[0.0 0.0 0.0; points[i,:]']
check=true
 while check

 xn=gp[end,1]
 yn=gp[end,2]
 zn=gp[end,3]
 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]


 xx=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 yy=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 zz=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
 gp=[gp;xx yy zz]
 temp=abs.([xx yy zz])
 v=[xx yy zz]
 hd=h*1.5


  if norm(temp-[0 0 pi])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[0 0 pi].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[0 pi 0])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[0 pi 0].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[pi 0 0])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[pi 0 0].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[0 pi pi])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[0 pi pi].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[pi 0 pi])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[pi 0 pi].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[pi pi 0])<hd
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[pi pi 0].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  elseif norm(temp-[pi pi pi])<1.5*h
   for k=1:3
    if temp[k]>0
     v[k]=v[k]/temp[k]
    end
   end
   v=[pi pi pi].*v
   gp=[gp;v]
   check=false
push!(gps,gp)
  end
 end
end

mesh=delaunay(points)
areas=[spherical_triangle_area(mesh.points[v, :],[0.0;0.0;0.0]) for v in eachrow(mesh.convex_hull)]
surface_area=4.0*pi*r^2

np=length(points[:,1])
nf=length(mesh.convex_hull[:,1])
area=Vector{Float64}(undef,np)
for i=1:np
area[i]=0.0
 for j=1:nf
  if i==mesh.convex_hull[j,1] || i==mesh.convex_hull[j,2] || i==mesh.convex_hull[j,3]
   area[i]=area[i]+areas[j]
  end
 end
end

total=sum(area)
area=area./total
area=area.*surface_area

area2=Vector{Float64}(undef,np)
for i=1:np
area2[i]=0.0
 for j=1:nf
  if i==mesh.convex_hull[j,1] || i==mesh.convex_hull[j,2] || i==mesh.convex_hull[j,3]
   area2[i]=area2[i]+areas[j]/3
  end
 end
end

volume=4/3*pi*r^3

mesh=Dict("points"=>mesh.points,
          "faces"=>mesh.convex_hull,
          "areas"=>areas,
          "area"=>area,
          "area2"=>area2,
          "surface_area"=>surface_area)

sys=Dict("mesh"=>mesh,
         "f"=>f,
         "g"=>g,
         "gps"=>gps,
         "sphere_volume"=>volume,
         "h"=>h)

return sys
end
#################################################################
#################################################################
function getcurvsallpaths(sys)

ngp=length(sys["gps"])

k=[]

 for i=1:ngp
  kgp=[]
  np=length(sys["gps"][i][:,1])
  temp=curv(sys["gps"][i][2,:])
  push!(kgp,temp)
   for j=3:np-1

    w1=temp["eigen_vectors"][1]
    w2=temp["eigen_vectors"][2]

    temp=curv(sys["gps"][i][j,:])

    v1=temp["eigen_vectors"][1]
    v2=temp["eigen_vectors"][2]

    e1=temp["eigen_values"][1]
    e2=temp["eigen_values"][2]

     if norm(cross(w1,v1))>norm(cross(w1,v2))
      temp["eigen_vectors"][1]=v2
      temp["eigen_vectors"][2]=v1
      temp["eigen_values"][1]=e2
      temp["eigen_values"][2]=e1
     end

    v1=temp["eigen_vectors"][1]
    v2=temp["eigen_vectors"][2]


     if dot(v1,w1)<0 
      temp["eigen_vectors"][1]=-v1
     end

     if dot(v2,w2)<0 
      temp["eigen_vectors"][2]=-v2
     end

    push!(kgp,temp)
   end
  push!(k,kgp)
 end


sys["k"]=k
return sys
end
#################################################################
#################################################################
function runvolumeallpaths(n,h,r)

sys=getallpaths(n,h,r)
sys=getcurvsallpaths(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengths(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumeallpaths(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumeallpaths(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
#Second correction term code, normal point sampling
#################################################################
#################################################################
function getarclengthsc2v1(sys)

h=0.00001

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]


for i=1:ngp
 a1=[]
 a2=[]
 sa=[]
 vol=[]
 start_area=sys["mesh"]["area"][i]
 push!(a1,sqrt(start_area/pi))
 push!(a2,sqrt(start_area/pi))
 np=length(sys["gps"][i][:,1])
  for j=3:np

   k1=sys["k"][i][j-2]["eigen_values"][1]
   k2=sys["k"][i][j-2]["eigen_values"][2]
   ds=sys["dgp"][i][j-1]

   d1=sys["k"][i][j-2]["eigen_vectors"][1]
   d2=sys["k"][i][j-2]["eigen_vectors"][2]
   start=start=sys["gps"][i][j-1,:]
   p1=start+.5*a1[end]*d1
   p2=start-.5*a1[end]*d1
   p3=start+.5*a2[end]*d2
   p4=start-.5*a2[end]*d2
  
   pk1=pathcurv(p1,h)
   pk2=pathcurv(p2,h)
   pk3=pathcurv(p3,h)
   pk4=pathcurv(p4,h)

   s1=getsign(start,p1,h,a1[end]*.5)
   s2=getsign(start,p2,h,a1[end]*.5)
   s3=getsign(start,p3,h,a2[end]*.5)
   s4=getsign(start,p4,h,a2[end]*.5)


   temp1=a1[end]*(1+k1*ds)+ds^2/2*(s1*pk1+s2*pk2)
   temp2=a2[end]*(1+k2*ds)+ds^2/2*(s3*pk3+s4*pk4)

   push!(sa,pi*temp1*temp2)
   push!(vol,pi*temp1*temp2*ds)

   push!(a1,temp1)
   push!(a2,temp2)
  end
 push!(arc1,a1)
 push!(arc2,a2)
 push!(surfarea,sa)
 push!(volumes,vol)
end


sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes


return sys
end
#################################################################
#################################################################
function getsign(r1,r2,h,d)

g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

 xn=r1[1]
 yn=r1[2]
 zn=r1[3]
 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]



 x1=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 y1=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 z1=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
 p1=[x1 y1 z1]

 h1=sqrt((xn-x1)^2+(yn-y1)^2+(zn-z1)^2)


 xn=r2[1]
 yn=r2[2]
 zn=r2[3]
 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]


 x2=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 y2=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 z2=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
 p2=[x2 y2 z2]

 h2=sqrt((xn-x2)^2+(yn-y2)^2+(zn-z2)^2)

 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 x3=xn+h2*grad[1]
 y3=yn+h2*grad[2]
 z3=zn+h2*grad[3]
 p3=[x3 y3 z3]

 h3=sqrt((xn-x3)^2+(yn-y3)^2+(zn-z3)^2)

length1=norm(p1-p3)
length2=norm(p1-p2)


if length2>length1
s=1
elseif length2<length1
s=-1
else
s=0
end

temp=abs.(r2')
 hd=d
  if norm(temp-[0 0 pi])<hd || norm(temp-[0 pi 0])<hd || norm(temp-[pi 0 0])<hd || norm(temp-[0 pi pi])<hd || norm(temp-[pi 0 pi])<hd || norm(temp-[pi pi 0])<hd
s=0
end

return s
end
#################################################################
#################################################################
function runvolumec2v1(n,h,r)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengthsc2v1(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumec2v1(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumec2v1(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
function pathcurv(r,h)

x1=r[1]
y1=r[2]
z1=r[3]

g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

 grad=-g(x1,y1,z1)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(x1+h*kx1/2,y1+h*ky1/2,z1+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(x1+h*kx2/2,y1+h*ky2/2,z1+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(x1+h*kx3,y1+h*ky3,z1+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]


 x2=x1+h/6*(kx1+2*kx2+2*kx3+kx4)
 y2=y1+h/6*(ky1+2*ky2+2*ky3+ky4)
 z2=z1+h/6*(kx1+2*kz2+2*kz3+kz4)

 grad=g(x1,y1,z1)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=g(x1+h*kx1/2,y1+h*ky1/2,z1+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=g(x1+h*kx2/2,y1+h*ky2/2,z1+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=g(x1+h*kx3,y1+h*ky3,z1+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]


 x0=x1+h/6*(kx1+2*kx2+2*kx3+kx4)
 y0=y1+h/6*(ky1+2*ky2+2*ky3+ky4)
 z0=z1+h/6*(kx1+2*kz2+2*kz3+kz4)


dx=(x2-x0)/2
dy=(y2-y0)/2
dz=(z2-z0)/2

ddx=x2-2*x1+x0
ddy=y2-2*y1+y0
ddz=z2-2*z1+z0

t0=(dx^2+dy^2+dz^2)^(3/2)
t1=(ddz*dy-ddy*dz)^2
t2=(ddx*dz-ddz*dx)^2
t3=(ddy*dx-ddx*dy)^2

k=sqrt(t1+t2+t3)/t0

return k
end
#################################################################
#################################################################
#Second correction term code, curv based point sampling
#################################################################
#################################################################
function getarclengthsc2v2(sys)

g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

h=0.00001

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]


for i=1:ngp
 a1=[]
 a2=[]
 sa=[]
 vol=[]
 start_area=sys["mesh"]["area"][i]
 push!(a1,sqrt(start_area/pi))
 push!(a2,sqrt(start_area/pi))
 np=length(sys["gps"][i][:,1])
  for j=3:np

   k1=sys["k"][i][j-2]["eigen_values"][1]
   k2=sys["k"][i][j-2]["eigen_values"][2]
   ds=sys["dgp"][i][j-1]

   d1=sys["k"][i][j-2]["eigen_vectors"][1]
   d2=sys["k"][i][j-2]["eigen_vectors"][2]
   start=start=sys["gps"][i][j-1,:]

   r1=1/abs(k1)
   r2=1/abs(k2)

   grad=-g(start[1],start[2],start[3])
   grad=grad./norm(grad)

   if k1>0
   c1=start-r1.*grad
   v1=-grad
   else
   c1=start+r1.*grad
   v1=grad
   end

   if k2>0
   c2=start-r2.*grad
   v2=-grad
   else
   c2=start+r2.*grad
   v2=grad
   end

   f1(t)=c1+r1*cos(t)*d1+r1*sin(t)*v1
   f2(t)=c2+r2*cos(t)*d2+r2*sin(t)*v2

   t1=1.5*pi
   t2=1.5*pi

   dt1=.5*a1[end]/r1
   dt2=.5*a2[end]/r2

   p1=f1(t1+dt1)
   p2=f1(t1-dt1)
   p3=f2(t2+dt2)
   p4=f2(t2-dt2)

   pk1=pathcurv(p1,h)
   pk2=pathcurv(p2,h)
   pk3=pathcurv(p3,h)
   pk4=pathcurv(p4,h)

   s1=getsign(start,p1,h,a1[end]*.5)
   s2=getsign(start,p2,h,a1[end]*.5)
   s3=getsign(start,p3,h,a2[end]*.5)
   s4=getsign(start,p4,h,a2[end]*.5)


   temp1=a1[end]*(1+k1*ds)+ds^2/2*(s1*pk1+s2*pk2)
   temp2=a2[end]*(1+k2*ds)+ds^2/2*(s3*pk3+s4*pk4)

   push!(sa,pi*temp1*temp2)
   push!(vol,pi*temp1*temp2*ds)

   push!(a1,temp1)
   push!(a2,temp2)
  end
 push!(arc1,a1)
 push!(arc2,a2)
 push!(surfarea,sa)
 push!(volumes,vol)
end


sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes


return sys
end
#################################################################
#################################################################
#################################################################
#################################################################
function runvolumec2v2(n,h,r)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengthsc2v2(sys)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumec2v2(nruns)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumec2v2(n,h,r)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
#Second correction term code, normal point sampling, path curv weighted distance
#################################################################
#################################################################
function getarclengthsc2v1weight(sys,aweight)

h=0.00001

ngp=length(sys["gps"])

arc1=[]
arc2=[]
surfarea=[]
volumes=[]


for i=1:ngp
 a1=[]
 a2=[]
 sa=[]
 vol=[]
 start_area=sys["mesh"]["area"][i]
 push!(a1,sqrt(start_area/pi))
 push!(a2,sqrt(start_area/pi))
 np=length(sys["gps"][i][:,1])
  for j=3:np

   k1=sys["k"][i][j-2]["eigen_values"][1]
   k2=sys["k"][i][j-2]["eigen_values"][2]
   ds=sys["dgp"][i][j-1]


   d1=sys["k"][i][j-2]["eigen_vectors"][1]
   d2=sys["k"][i][j-2]["eigen_vectors"][2]
   start=start=sys["gps"][i][j-1,:]

   kweight=aweight^-abs(pathcurv(start,h))

   p1=start+.5*a1[end]*d1*kweight
   p2=start-.5*a1[end]*d1*kweight
   p3=start+.5*a2[end]*d2*kweight
   p4=start-.5*a2[end]*d2*kweight

   pk1=pathcurv(p1,h)
   pk2=pathcurv(p2,h)
   pk3=pathcurv(p3,h)
   pk4=pathcurv(p4,h)

   s1=getsignweight(start,p1,h,a1[end]*.5*kweight)
   s2=getsignweight(start,p2,h,a1[end]*.5*kweight)
   s3=getsignweight(start,p3,h,a2[end]*.5*kweight)
   s4=getsignweight(start,p4,h,a2[end]*.5*kweight)

   temp1=a1[end]*(1+k1*ds)+ds^2/2*(s1*pk1+s2*pk2)
   temp2=a2[end]*(1+k2*ds)+ds^2/2*(s3*pk3+s4*pk4)

   push!(sa,pi*temp1*temp2)
   push!(vol,pi*temp1*temp2*ds)

   push!(a1,temp1)
   push!(a2,temp2)
  end
 push!(arc1,a1)
 push!(arc2,a2)
 push!(surfarea,sa)
 push!(volumes,vol)
end


sys["arc1"]=arc1
sys["arc2"]=arc2
sys["surfarea"]=surfarea
sys["volumes"]=volumes


return sys
end
#################################################################
#################################################################
function getsignweight(r1,r2,h,d)

g(x,y,z) = [-sin(x), -sin(y), -sin(z)]

 xn=r1[1]
 yn=r1[2]
 zn=r1[3]
 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]



 x1=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 y1=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 z1=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
 p1=[x1 y1 z1]

 h1=sqrt((xn-x1)^2+(yn-y1)^2+(zn-z1)^2)


 xn=r2[1]
 yn=r2[2]
 zn=r2[3]
 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 kx1=grad[1]
 ky1=grad[2]
 kz1=grad[3]
 grad=-g(xn+h*kx1/2,yn+h*ky1/2,zn+h*kz1/2)
 grad=grad./norm(grad)
 kx2=grad[1]
 ky2=grad[2]
 kz2=grad[3]
 grad=-g(xn+h*kx2/2,yn+h*ky2/2,zn+h*kz2/2)
 grad=grad./norm(grad)
 kx3=grad[1]
 ky3=grad[2]
 kz3=grad[3]
 grad=-g(xn+h*kx3,yn+h*ky3,zn+h*kz3)
 grad=grad./norm(grad)
 kx4=grad[1]
 ky4=grad[2]
 kz4=grad[3]


 x2=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 y2=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 z2=zn+h/6*(kx1+2*kz2+2*kz3+kz4)
 p2=[x2 y2 z2]

 h2=sqrt((xn-x2)^2+(yn-y2)^2+(zn-z2)^2)

 grad=-g(xn,yn,zn)
 grad=grad./norm(grad)
 x3=xn+h2*grad[1]
 y3=yn+h2*grad[2]
 z3=zn+h2*grad[3]
 p3=[x3 y3 z3]

 h3=sqrt((xn-x3)^2+(yn-y3)^2+(zn-z3)^2)

length1=norm(p1-p3)
length2=norm(p1-p2)


if length2>length1
s=1
elseif length2<length1
s=-1
else
s=0
end

temp=abs.(r2')	
 hd=d
  if norm(temp-[0 0 pi])<hd || norm(temp-[0 pi 0])<hd || norm(temp-[pi 0 0])<hd || norm(temp-[0 pi pi])<hd || norm(temp-[pi 0 pi])<hd || norm(temp-[pi pi 0])<hd
display("error: sampled point too close to saddle")
end

return s
end
#################################################################
#################################################################
function runvolumec2v1weight(n,h,r,aweight)

sys=getpaths(n,h,r)
sys=getcurvs(sys)
sys=getgpsegmentlengths(sys)
sys=getarclengthsc2v1weight(sys,aweight)
sys=getvolume(sys)
volume=sys["volume"]

return sys
end
#################################################################
#################################################################
function errortablevolumec2v1weight(nruns,aweight)

h=.005
r=0.1
ex=(2*pi)^3

@printf "\n Relative error of volume of cube"
@printf "\n Exact volume is (2*pi)^3=248.05"
@printf "\n ngps is number of gradient paths"
@printf "\n Radius of starting sphere is r=0.1"
@printf "\n ngps        vol           rel err"

 for i=1:nruns

  n=100*i
  sys=runvolumec2v1weight(n,h,r,aweight)
  vol=sys["volume"]
  err=abs(vol-ex)/ex
  ngp=length(sys["gps"])
  @printf "\n %g       %g         %g" ngp  vol  err

 end

end
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
