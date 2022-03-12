##### Create a stable VAR(2) #### 

## Experiment setting

mc     <- 1000 ## Monte Carlo replications
n_lhs  <- 500 ## Number of parameters' configurations for initialization


angles <- seq(from = 0.01, to = (pi/2-0.01),length.out = 29) ## Rotation angles in Givens Rotation matrices
comb.angles <- permutations(n = 29,r = 3,v = angles,repeats.allowed = T)

## An arbitrary 2x2 Mixing Matrix
Btrue2Dz <- matrix(c(1.14, -0.38,
                     0,    1.26), nrow = 2, ncol = 2)

## An arbitrary 3x3 Mixing Matrix
## This Btrue3D has an inverse (the unmixing matrix containing contemporaneous relationships) 
## with an entry equal to 0
Btr3D  <- matrix(c(0.9,   0.150, 0.65,
                   -0.75,   1.13, 0.22,
                   0.21,  -0.53, 1.5)  , nrow = 3, ncol = 3)

Wtrue3D   <- t(solve(Btr3D))
Wtrue3D[3,1]<-0
Btrue3D   <- solve(Wtrue3D)


## Empty objects for finding best initial conditions

winit         <- vector('list', n_lhs)
dc_init       <- vector('list', n_lhs)
ica_init      <- vector('list', n_lhs)
pml_init      <- vector('list', n_lhs)
ilm.index     <- double(length = n_lhs)
B.init        <- vector("list",n_lhs)
FR_ica_nt     <- double(length = n_lhs)
FR_dc_nt      <- double(length = n_lhs)
Wica          <- vector("list", n_lhs)
Wscaled_ica  <- vector("list", n_lhs)
Wscaled_dc   <- vector("list", n_lhs)
Wscal_dc     <- vector("list", mc)
Wscal_ica    <- vector("list", mc)
FR_pml_nt     <- double(length = n_lhs)

B.init          <- vector("list", n_lhs)
DC_MD          <- vector("list", n_lhs)
icares_MD      <- vector("list", n_lhs)

A.ica        <- vector("list", mc)
A.dc         <- vector("list", mc)
A.cvm        <- vector("list", mc)
A.pml        <- vector("list", mc)
A.ID.ica     <- vector('list',mc)
A.ID.dc      <- vector('list',mc)
A.ID.cvm     <- vector("list", mc)
A.ID.pml     <- vector("list", mc)
perfm_ica    <- double(length = mc)
perfm_dc     <- double(length = mc)
perfm_cvm    <- double(length = mc)
pvalJB1      <- double(length = mc)
pvalJB2      <- double(length = mc)
pvalJB3      <- double(length = mc)
boot.ica     <- vector("list", mc)
boot.dc      <- vector("list", mc)
boot.cvm     <- vector("list", mc)
Frb_m_ica            <- double(length = length(SEQ))
Frb_m_dc             <- double(length = length(SEQ))
Frb_m_cvm            <- double(length = length(SEQ))
Frb_sd_dc            <- double(length = length(SEQ))
Frb_sd_ica           <- double(length = length(SEQ))
Frb_sd_cvm           <- double(length = length(SEQ))
ress                 <- vector('list',length = length(SEQ))
pippo                <- vector('list',length = length(SEQ)) 


## PML initialization
p.start    <- vector("list",n_lhs)
ilm.index        <- double(length = n_lhs)
## PML estimation ####
C          <- vector("list",mc)
u_chol     <- vector("list",mc)
Orth.true  <- vector("list",mc)
Orth.max   <- vector("list",mc)
var.matrix <- vector("list",mc)
t.stat     <- double(length = mc)

mean.perf.ml <- double(length = length(SEQ))
sd.perf      <- double(length = length(SEQ))
