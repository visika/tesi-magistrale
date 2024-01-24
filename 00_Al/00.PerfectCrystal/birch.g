# define the Birch-Murnaghan equation of state birch(x):

birch(x) = emin + 1.5*v0*k0/1602 * ( 0.75*(1+2*xsi)*(v0/x)**(4.0/3.0) \
-xsi/2*(v0/x)**2 - 3.0/2.0*(1+xsi)*(v0/x)**(2.0/3.0) + 0.5*(xsi+1.5) )

dbirch(x) = 1.5*v0*k0 * ( -(v0/x)**(4.0/3.0)/x + (v0/x)**(2.0/3.0)/x )

# set initial guess parameters for the fit:

v0=16
emin=-4
k0=100
xsi=0.0

# now perform the fit with the command:

# fit birch(x) 'filename' via emin,v0,k0

# where filename is the name of the file containig the energies

