module Ctes
    implicit none
    integer, parameter:: N=10                           !Number of PP
    real, parameter:: Df=2.1                            !Fractal dimension
    real, parameter:: kf=.9                         !Fractal prefactor
    real, parameter:: rp_g=5.0                          !Geometric mean PP
    real, parameter:: rp_gstd=1.00                       !Geometric PP standard deviation
    integer, parameter:: Quantity_aggregates = 100       !Quantity of aggregates to be generated
    integer, parameter:: Ext_case = 0                    !Activate extreme cases
    real, parameter:: Nsubcl_perc = .1                  !controls the subclusters size (keep always in 10%, only increase when PC is not able to work)
    real, parameter:: tol_ov=10.**(-1.)                  !Tolerance to overlapping
    real, parameter :: pi=4.0*atan(1.0)
    real R(N)                                            !Primary particles radii and mass
    real X(N), Y(N), Z(N)                                !Coordinates of primary particles
    integer:: iter, i, j, k
end module Ctes
