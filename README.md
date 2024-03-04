# gponlyintegration3d
Gradient path only volume integration methods in 3D.

f(x,y,z) = 3.0 + cos(x) + cos(y) + cos(z)

Functions:
-----------------------------------------

getpaths(n,h,r) generates the gradient paths. n is number of paths. h is step size for RK4. r is radius of starting circle. Gradient paths that intersect saddle points are excluded.

plotmesh(sys) plots the mesh where the gradient paths are seeded.

plotgps(sys) plots the gradient paths.

curv(r) calculates the principal curvatures and directions of the isosurface at point r.

getcurvs(sys) calculates the principal curvatures and directions along all of the gradient paths.

getgpsegmentlengths(sys) calculates the arc lengths between points along all the the gradient paths.

getarclengths(sys) uses the original recurrence relation to calculate all of the arc lengths, gradient bundle cross section areas, and segment volumes along each gradient path.

getvolume(sys) sums the segment volumes down each gradient path and add the starting sphere volume to get the approximate total volume of the cube.

runvolume(n,h,r) runs all steps to approximate the volume of the cube using the original recurrence relation without corrections. n is number of paths. h is step size for RK4. r is radius of starting circle.

errortablevolume(nruns) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses the original recurrence relation without corrections. Cases start with 100 gradient paths and increase by 100 paths each time.

getarclengthsc1v1(sys) uses the recurrence relation including the first correction term (the change in curvature of the isosurfaces in the principal directions) to calculate all of the arc lengths, gradient bundle cross section areas, and segment volumes along each gradient path. The derivative of curvature is approximated with FDM. The sampled off path points are calculated in the normal direction to the gradient path.

runvolumec1v1(n,h,r) runs all steps to approximate the volume of the cube using getarclengthsc1v1(sys). n is number of paths. h is step size for RK4. r is radius of starting circle.

errortablevolumec1v1(nruns) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses getarclengthsc1v1(sys). Cases start with 100 gradient paths and increase by 100 paths each time.

getarclengthsc1v2(sys) uses the recurrence relation including the first correction term (the change in curvature of the isosurfaces in the principal directions) to calculate all of the arc lengths, gradient bundle cross section areas, and segment volumes along each gradient path. The derivative of curvature is approximated with FDM. The sampled off path points are calculated using the curvature of the isosurface in the principal directions.


runvolumec1v2(n,h,r) runs all steps to approximate the volume of the cube using getarclengthsc1v2(sys). n is number of paths. h is step size for RK4. r is radius of starting circle.

errortablevolumec1v2(nruns) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses getarclengthsc1v2(sys). Cases start with 100 gradient paths and increase by 100 paths each time.

getallpaths(n,h,r) generates the gradient paths. n is number of paths. h is step size for RK4. r is radius of starting circle. No paths are excluded. Gradient paths that intersect saddle points are included.

getcurvsallpaths(sys) is the same as getcurvs(sys), but it is modified to avoid returning NaN for gradient paths intersect saddles.

runvolumeallpaths(n,h,r) is the same as runvolume(n,h,r), but uses getallpaths(n,h,r) to include all paths.

errortablevolumeallpaths(nruns) is the same as errortablevolume(nruns), but it uses getallpaths(n,h,r) to include all paths.

getarclengthsc2v1(sys) uses the recurrence relation including the second correction term (the curvature of the gradient paths that bound the gradient bundle) to calculate all of the arc lengths, gradient bundle cross section areas, and segment volumes along each gradient path. The curvature of the gradient paths is approximated with FDM. The sampled off path points are calculated in the normal direction to the gradient path.

getsign(r1,r2,h,d) determines if the curvature of the bounding gradient path through sampled point r2 increases or decreases the volume of the segment compared to assuming the path is straight. r1 is where the original gradient path intersects the isosurface. h is a step size for the test comparing Euler's method to RK4. d is the minimum distance allowed for the proximity of the sampled point to a sampled point. It returns 1 if the volume is increased due to curvature, -1 if the volume is decreased due to curvature, and 0 if the sampled point is too close to a saddle.

runvolumec2v1(n,h,r) runs all steps to approximate the volume of the cube using getarclengthsc2v1(sys). n is number of paths. h is step size for RK4. r is radius of starting circle.

errortablevolumec2v1(nruns) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses getarclengthsc2v1(sys). Cases start with 100 gradient paths and increase by 100 paths each time.

pathcurv(r,h) calculates the curvature of the gradient path at point r using step size h for the FDM.

getarclengthsc2v2(sys) uses the recurrence relation including the second correction term (the curvature of the gradient paths that bound the gradient bundle) to calculate all of the arc lengths, gradient bundle cross section areas, and segment volumes along each gradient path. The curvature of the gradient paths is approximated with FDM. The sampled off path points are calculated using the curvature of the isosurface in the principal directions.

runvolumec2v1(n,h,r) runs all steps to approximate the volume of the cube using getarclengthsc2v2(sys). n is number of paths. h is step size for RK4. r is radius of starting circle.

errortablevolumec2v2(nruns) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses getarclengthsc2v2(sys). Cases start with 100 gradient paths and increase by 100 paths each time.

getarclengthsc2v1weight(sys,aweight) is the same as getarclengthsc2v1(sys), but it uses the weight function f(k)=aweight^-abs(k) to adjust the distance for sampling off path points where k is the curvature of the gradient path.

getsignweight(r1,r2,h,d) is the same as getsign(r1,r2,h,d), but instead of returning 0 if the sampled point is too close to a saddle it returns an error message. If the error message is given, aweight must be increased to insure sampled points don't get too close to a saddle.

runvolumec2v1weight(n,h,r,aweight) runs all steps to approximate the volume of the cube using getarclengthsc2v1weight(sys). n is number of paths. h is step size for RK4. r is radius of starting circle. Uses weight function f(k)=aweight^-abs(k) to adjust the distance for sampling off path points.

errortablevolumec2v1weight(nruns,aweight) runs several cases of different numbers of gradient paths, calculates the approximate volume and relative error, produces a table of the results. Uses getarclengthsc2v1weight(sys). Cases start with 100 gradient paths and increase by 100 paths each time.
