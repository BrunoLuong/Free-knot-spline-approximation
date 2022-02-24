# Free-knot-spline-approximation

 The purpose of this function is to provide a flexible and robust fit to  one-dimensional data using free-knot splines. The knots are free and  able to cope with rapid change in the underlying model. Knot removal  strategy is used to fit with only a small number of knots. 
Optional L2-regularization on the derivative of the spline function can be used to enforce the smoothness. 
Shape preserving approximation can be enforced by specifying the  lower and upper bounds of the derivative(s) of the spline function on  sub-intervals. Furthermore specific values of the spline function and  its derivative can be specified on a set of discrete data points. 
I did not test QUADPROG engine, but I have implemented it. Any feedback is welcome. 