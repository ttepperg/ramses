module metal_yields
   use amr_commons, only: dp
   implicit none

   !---------------------------------------
   ! Facilitates yield tables and interpolation function. Includes
   ! yield tables for C,N,O,Al,Si,Ca and Fe from AGB winds and SNII.
   ! Yield tables taken from NuGrid II, Ritter et al. (2018).
   !---------------------------------------
   integer, parameter :: nmet=5
   integer, parameter :: nmass_AGB=8
   integer, parameter :: nmass_SNII=4
  
   ! New yields taken from Seitenzahl et al. (2013). Sum of all isotopes using model N100, see table 2. 
   ! Compare to old yields (Raiteri et al. 1996) 0.63 Msun Fe + 0.13Msun O (e.g. Theilemann 1986)
   ! Mass (Msun) of enriched material from SNIa 
   real(dp) :: SNIaFe=7.40d-01  ! Iron; old value 0.7d0*scale_m  !see Hopkins, 0.63d0 Raiteri
   real(dp) :: SNIaO=1.01d-01   ! Oxygen, old value 0.13d0*scale_m
   real(dp) :: SNIaC=3.04d-03   ! Carbon
   real(dp) :: SNIaN=3.21d-06   ! Nitrogen
   real(dp) :: SNIaMg=1.54d-02  ! Magnesium
   real(dp) :: SNIaAl=6.74d-04  ! Aluminium
   real(dp) :: SNIaSi=2.87d-01  ! Silicon
   real(dp) :: SNIaCa=1.48d-02  ! Calcium
   real(dp) :: SNIaEu=0.d0      ! Europium (not known)
  
   ! Mass (Msun) of enriched material from neutron star mergers.
   real(dp) :: MEuNSNS=1.0d-5   ! Europium. See Cote et al. (2018).

   ! Tables from Nugrid.
   real(dp), dimension(1:nmet) :: table_met = (/ 0.0001, 0.001, 0.006, 0.01, 0.02/)
   real(dp), dimension(1:nmass_AGB) :: table_mass_AGB = (/ 1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 /)
   real(dp), dimension(1:nmass_SNII) :: table_mass_SNII = (/ 12.0, 15.0, 20.0, 25.0/)

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Fe_AGB = transpose(reshape( &
      (/ &
         6.453d-07, 1.451d-06, 1.881d-06, 3.073d-06, 4.429d-06, 5.754d-06, 7.012d-06, 8.238d-06, &
         5.977d-06, 1.473d-05, 1.931d-05, 3.116d-05, 4.448d-05, 5.764d-05, 7.016d-05, 8.259d-05, &
         3.733d-05, 8.917d-05, 1.186d-04, 1.963d-04, 2.701d-04, 3.493d-04, 4.266d-04, 5.037d-04, &
         3.169d-04, 7.400d-04, 9.961d-04, 1.676d-03, 2.260d-03, 2.933d-03, 3.590d-03, 4.227d-03, &
         6.325d-04, 1.437d-03, 1.982d-03, 3.379d-03, 4.543d-03, 5.909d-03, 7.231d-03, 8.526d-03 &
      /), [nmass_AGB, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_O_AGB = transpose(reshape(&
      (/ &
         1.835d-03, 8.626d-03, 9.952d-03, 2.781d-03, 4.080d-03, 1.833d-04, 2.239d-04, 3.249d-03, &
         5.767d-04, 5.937d-03, 1.184d-02, 4.024d-03, 6.280d-03, 1.217d-03, 1.252d-03, 3.756d-03, &
         1.921d-03, 7.652d-03, 1.330d-02, 1.577d-02, 1.502d-02, 1.426d-02, 1.410d-02, 1.550d-02, &
         2.127d-03, 5.908d-03, 1.389d-02, 2.232d-02, 1.982d-02, 1.997d-02, 1.943d-02, 1.835d-02, &
         4.243d-03, 1.001d-02, 1.911d-02, 3.604d-02, 3.844d-02, 3.923d-02, 3.799d-02, 4.214d-02 &
      /), [nmass_AGB, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_N_AGB = transpose(reshape(&
      (/ &
         7.763d-05, 5.710d-05, 3.870d-05, 3.408d-05, 1.019d-02, 4.692d-03, 2.360d-03, 8.018d-03, &
         1.985d-05, 1.036d-04, 1.711d-04, 2.213d-04, 7.812d-03, 6.323d-03, 5.150d-03, 7.922d-03, &
         8.937d-05, 4.067d-04, 7.848d-04, 1.694d-03, 2.854d-03, 7.439d-03, 1.453d-02, 1.699d-02, &
         3.519d-04, 1.248d-03, 1.931d-03, 3.817d-03, 5.838d-03, 1.559d-02, 2.974d-02, 3.151d-02, &
         6.713d-04, 2.328d-03, 3.561d-03, 7.173d-03, 1.041d-02, 1.678d-02, 4.112d-02, 4.746d-02 &
      /), [nmass_AGB, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Si_AGB = transpose(reshape(&
      (/ &
         7.645d-07, 1.791d-06, 2.559d-06, 3.970d-06, 6.714d-06, 7.136d-06, 1.495d-05, 9.418d-06, &
         6.568d-06, 1.642d-05, 2.203d-05, 3.545d-05, 5.208d-05, 6.518d-05, 8.007d-05, 9.254d-05, &
         4.094d-05, 9.836d-05, 1.316d-04, 2.184d-04, 3.014d-04, 3.845d-04, 4.729d-04, 5.591d-04, &
         1.787d-04, 4.176d-04, 5.643d-04, 9.529d-04, 1.293d-03, 1.671d-03, 2.055d-03, 2.421d-03, &
         3.568d-04, 8.105d-04, 1.121d-03, 1.919d-03, 2.603d-03, 3.362d-03, 4.098d-03, 4.834d-03 &
      /), [nmass_AGB, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Ca_AGB = transpose(reshape(&
       (/ &
         5.590d-08, 1.240d-07, 1.607d-07, 2.630d-07, 3.790d-07, 4.920d-07, 5.993d-07, 7.038d-07, &
         5.108d-07, 1.257d-06, 1.647d-06, 2.666d-06, 3.806d-06, 4.929d-06, 5.997d-06, 7.057d-06, &
         3.192d-06, 7.613d-06, 1.011d-05, 1.678d-05, 2.310d-05, 2.986d-05, 3.647d-05, 4.305d-05, &
         1.647d-05, 3.844d-05, 5.161d-05, 8.681d-05, 1.174d-04, 1.524d-04, 1.865d-04, 2.196d-04, &
         3.288d-05, 7.464d-05, 1.027d-04, 1.749d-04, 2.359d-04, 3.070d-04, 3.758d-04, 4.430d-04 &
       /), [nmass_AGB, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Al_AGB = transpose(reshape(&
       (/ &
         3.506d-08, 3.905d-07, 8.674d-07, 6.366d-07, 2.057d-06, 2.275d-06, 1.841d-06, 1.037d-06, &
         2.997d-07, 7.956d-07, 1.619d-06, 2.118d-06, 3.509d-06, 4.639d-06, 7.882d-06, 5.905d-06, &
         1.682d-06, 4.163d-06, 5.781d-06, 1.022d-05, 1.527d-05, 1.637d-05, 2.226d-05, 2.862d-05, &
         1.429d-05, 3.373d-05, 4.624d-05, 7.894d-05, 1.151d-04, 1.432d-04, 1.770d-04, 2.175d-04, &
         2.851d-05, 6.518d-05, 9.123d-05, 1.576d-04, 2.284d-04, 2.850d-04, 3.422d-04, 4.029d-04 &
       /), [nmass_AGB, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Mg_AGB = transpose(reshape(&
       (/ &
         1.005d-06, 5.403d-05, 9.078d-05, 3.960d-05, 1.193d-04, 2.429d-05, 1.937d-05, 4.907d-05, &
         7.271d-06, 3.665d-05, 1.105d-04, 6.605d-05, 1.570d-04, 9.500d-05, 1.083d-04, 1.736d-04, &
         4.351d-05, 1.191d-04, 1.905d-04, 3.405d-04, 4.288d-04, 4.160d-04, 5.799d-04, 1.017d-03, &
         1.657d-04, 3.927d-04, 5.835d-04, 1.099d-03, 1.549d-03, 1.784d-03, 2.296d-03, 3.465d-03, &
         3.308d-04, 7.573d-04, 1.126d-03, 2.166d-03, 3.451d-03, 3.713d-03, 4.115d-03, 4.980d-03 &
       /), [nmass_AGB, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Eu_AGB = transpose(reshape(&
       (/ &
         3.857E-12, 6.504E-13, 9.461E-13, 1.386E-12, 2.736E-12, 1.848E-12, 1.967E-12, 2.290E-12, &
         1.927E-12, 4.903E-12, 7.150E-12, 9.971E-12, 2.064E-11, 1.750E-11, 1.980E-11, 2.305E-11, &
         1.113E-11, 2.769E-11, 3.974E-11, 6.642E-11, 9.837E-11, 1.028E-10, 1.263E-10, 1.434E-10, &
         9.444E-11, 2.202E-10, 3.170E-10, 5.874E-10, 7.074E-10, 8.725E-10, 1.037E-09, 1.202E-09, &
         1.886E-10, 4.264E-10, 5.879E-10, 1.051E-09, 1.432E-09, 1.779E-09, 2.114E-09, 2.441E-09 &
       /), [nmass_AGB, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_C_AGB = transpose(reshape(&
       (/ &
         8.708E-03, 2.147E-02, 2.356E-02, 8.884E-03, 1.988E-03, 7.856E-04, 2.937E-04, 1.797E-03, &
         7.063E-04, 1.337E-02, 2.608E-02, 8.723E-03, 4.042E-03, 1.234E-03, 4.238E-04, 1.284E-03, &
         2.736E-04, 7.715E-03, 1.786E-02, 1.856E-02, 9.100E-03, 8.749E-04, 1.498E-03, 1.240E-03, &
         7.294E-04, 4.655E-03, 1.852E-02, 2.990E-02, 1.991E-02, 9.310E-03, 3.304E-03, 2.154E-03, &
         1.353E-03, 4.750E-03, 1.894E-02, 3.956E-02, 3.097E-02, 2.335E-02, 1.977E-03, 2.462E-03 &
       /), [nmass_AGB, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Fe_SNII = transpose(reshape(&
       (/ &
         2.443d-01, 5.073d-02, 7.646d-05, 2.518d-05, &
         1.692d-01, 5.407d-02, 4.089d-02, 2.326d-04, &
         1.863d-01, 1.739d-01, 3.275d-01, 1.036d-03, &
         2.032d-01, 7.515d-02, 2.745d-01, 8.288d-03, &
         1.531d-01, 5.561d-02, 1.493d-02, 7.470d-03 &
       /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_O_SNII = transpose(reshape(&
       (/ &
         2.828d-01, 1.176d+00, 2.043d+00, 7.062d-01, &
         3.459d-01, 1.148d+00, 1.871d+00, 7.875d-01, &
         1.879d-01, 8.612d-01, 1.442d+00, 8.904d-01, &
         1.194d-01, 8.829d-01, 1.314d+00, 8.494d-01, &
         1.205d-01, 8.970d-01, 1.504d+00, 8.446d-01 &
       /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_N_SNII = transpose(reshape(&
       (/ &
         2.160d-04, 3.047d-04, 2.320d-02, 2.813d-02, &
         2.047d-03, 2.664d-03, 3.984d-03, 5.828d-03, &
         1.172d-02, 1.417d-02, 2.151d-02, 2.986d-02, &
         1.993d-02, 2.453d-02, 3.150d-02, 4.422d-02, &
         4.055d-02, 4.235d-02, 5.627d-02, 3.436d-02 &
       /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Si_SNII = transpose(reshape(&
       (/ &
          9.289d-02, 1.078d-01, 1.138d-01, 2.063d-02, &
          9.255d-02, 9.346d-02, 2.482d-01, 3.416d-02, &
          8.842d-02, 1.468d-01, 3.725d-01, 2.865d-02, &
          8.212d-02, 2.248d-01, 8.855d-02, 3.778d-02, &
          5.401d-02, 1.902d-01, 1.487d-01, 5.818d-02 &
       /), [nmass_SNII, nmet] &
   ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Ca_SNII = transpose(reshape(&
      (/ &
         8.686d-03, 5.152d-03, 1.363d-03, 2.154d-06, &
         6.655d-03, 5.118d-03, 1.390d-02, 3.186d-05, &
         8.256d-03, 7.418d-03, 2.550d-02, 9.611d-05, &
         8.661d-03, 1.717d-03, 5.985d-04, 7.348d-04, &
         5.745d-03, 2.537d-03, 1.102d-03, 4.162d-04 &
      /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Al_SNII = transpose(reshape(&
      (/ &
         6.219d-05, 1.157d-03, 1.505d-03, 1.315d-04, &
         2.466d-04, 9.427d-04, 1.540d-03, 3.695d-04, &
         2.530d-04, 1.299d-03, 1.517d-03, 1.826d-03, &
         4.564d-04, 4.246d-03, 3.031d-03, 2.227d-03, &
         7.659d-04, 2.780d-03, 5.394d-03, 3.650d-03 &
      /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Mg_SNII = transpose(reshape(&
      (/ &
         1.638E-02, 7.562E-02, 1.160E-01, 2.178E-02, &
         2.463E-02, 8.823E-02, 1.019E-01, 2.791E-02, &
         1.122E-02, 6.267E-02, 7.786E-02, 4.062E-02, &
         1.109E-02, 6.908E-02, 7.642E-02, 4.494E-02, &
         1.314E-02, 3.892E-02, 1.366E-01, 9.823E-02 &
      /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Eu_SNII = transpose(reshape(&
      (/ &
         3.856E-12, 4.834E-12, 5.922E-12, 7.524E-12, &
         3.948E-11, 4.893E-11, 6.355E-11, 8.307E-11, &
         2.071E-10, 2.645E-10, 3.369E-10, 4.567E-10, &
         1.804E-09, 2.176E-09, 2.504E-09, 2.892E-09, &
         3.718E-09, 3.717E-09, 4.520E-09, 2.807E-09 &
      /), [nmass_SNII, nmet] &
   ))
   
   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_C_SNII = transpose(reshape(&
       (/ &
         1.226E-01, 1.646E-01, 1.956E-01, 2.037E-01, &
         1.200E-01, 1.537E-01, 2.222E-01, 2.115E-01, &
         1.015E-01, 1.655E-01, 1.473E-01, 2.206E-01, &
         8.914E-02, 1.575E-01, 2.294E-01, 2.238E-01, &
         7.646E-02, 1.472E-01, 1.923E-01, 2.035E-01 &
       /), [nmass_SNII, nmet] &
   ))
contains

!---------------------------------------
function AGB_Fe_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Fe_yield
   AGB_Fe_yield = interp_yield(table_mass_AGB, table_met, table_yield_Fe_AGB, mass, met)

END FUNCTION AGB_Fe_yield

!---------------------------------------
FUNCTION AGB_O_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_O_yield
   AGB_O_yield = interp_yield(table_mass_AGB, table_met, table_yield_O_AGB, mass, met)
END function AGB_O_yield

!---------------------------------------
FUNCTION AGB_N_yield(mass, met)
   !---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_N_yield
   AGB_N_yield = interp_yield(table_mass_AGB, table_met, table_yield_N_AGB, mass, met)
END FUNCTION AGB_N_yield

!---------------------------------------
FUNCTION AGB_Si_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Si_yield
   AGB_Si_yield = interp_yield(table_mass_AGB, table_met, table_yield_Si_AGB, mass, met)
END FUNCTION AGB_Si_yield

!---------------------------------------
FUNCTION AGB_Ca_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Ca_yield
   AGB_Ca_yield = interp_yield(table_mass_AGB, table_met, table_yield_Ca_AGB, mass, met)
END FUNCTION AGB_Ca_yield

!---------------------------------------
FUNCTION AGB_Al_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Al_yield
   AGB_Al_yield = interp_yield(table_mass_AGB, table_met, table_yield_Al_AGB, mass, met)
END FUNCTION AGB_Al_yield

!---------------------------------------
FUNCTION AGB_Mg_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Mg_yield
   AGB_Mg_yield = interp_yield(table_mass_AGB, table_met, table_yield_Mg_AGB, mass, met)
END FUNCTION AGB_Mg_yield

!---------------------------------------
FUNCTION AGB_Eu_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_Eu_yield
   AGB_Eu_yield = interp_yield(table_mass_AGB, table_met, table_yield_Eu_AGB, mass, met)
END FUNCTION AGB_Eu_yield

!---------------------------------------
FUNCTION AGB_C_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: AGB_C_yield
   AGB_C_yield = interp_yield(table_mass_AGB, table_met, table_yield_C_AGB, mass, met)
END FUNCTION AGB_C_yield

!---------------------------------------
FUNCTION SNII_Fe_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Fe_yield
   SNII_Fe_yield = interp_yield(table_mass_SNII, table_met, table_yield_Fe_SNII, mass, met)
END FUNCTION SNII_Fe_yield

!---------------------------------------
FUNCTION SNII_O_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_O_yield

   SNII_O_yield = interp_yield(table_mass_SNII, table_met, table_yield_O_SNII, mass, met)
END FUNCTION SNII_O_yield

!---------------------------------------
FUNCTION SNII_N_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_N_yield

   SNII_N_yield = interp_yield(table_mass_SNII, table_met, table_yield_N_SNII, mass, met)
END FUNCTION SNII_N_yield

!---------------------------------------
FUNCTION SNII_Si_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Si_yield

   SNII_Si_yield = interp_yield(table_mass_SNII, table_met, table_yield_Si_SNII, mass, met)
END FUNCTION SNII_Si_yield

!---------------------------------------
FUNCTION SNII_Ca_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Ca_yield

   SNII_Ca_yield = interp_yield(table_mass_SNII, table_met, table_yield_Ca_SNII, mass, met)
END FUNCTION SNII_Ca_yield

!---------------------------------------
FUNCTION SNII_Al_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Al_yield

   SNII_Al_yield = interp_yield(table_mass_SNII, table_met, table_yield_Al_SNII, mass, met)
END FUNCTION SNII_Al_yield

!---------------------------------------
FUNCTION SNII_Mg_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Mg_yield

   SNII_Mg_yield = interp_yield(table_mass_SNII, table_met, table_yield_Mg_SNII, mass, met)
END FUNCTION SNII_Mg_yield

!---------------------------------------
FUNCTION SNII_Eu_yield(mass, met)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_Eu_yield

   SNII_Eu_yield = interp_yield(table_mass_SNII, table_met, table_yield_Eu_SNII, mass, met)
END FUNCTION SNII_Eu_yield

!---------------------------------------
FUNCTION SNII_C_yield(mass, met)
!---------------------------------------
   real(dp), intent(in) :: mass, met
   real(dp) :: SNII_C_yield
   SNII_C_yield = interp_yield(table_mass_SNII, table_met, table_yield_C_SNII, mass, met)
END FUNCTION SNII_C_yield

function interp_yield(table_mass, table_met, table_yield, mass, met)
   real(dp), dimension(:), intent(in) :: table_met
   real(dp), dimension(:), intent(in) :: table_mass
   real(dp), dimension(:, :), intent(in) :: table_yield
   real(dp), intent(in) :: met, mass
   real(dp) :: interp_yield

   integer :: nmet, nmass
   integer :: imet, imass

   real(dp) :: x1, x2, y1, y2, f11, f21, f22, f12, a1, a2, a3, a4
   nmet = size(table_met)
   nmass = size(table_mass)

   ! Find index right of metallicity.
   imet = 2
   do while ((table_met(imet) < met) .AND. (imet <= nmet-1))
      imet = imet+1
   end do
    
   ! Find index right of mass.
   imass = 2
   do while ((table_mass(imass) < mass) .AND. (imass <= nmass-1))
      imass = imass+1
   end do

   ! ----------------
   ! Interpolate
   ! Store values around point.
   x1 = table_met(imet-1)
   x2 = table_met(imet)
   y1 = table_mass(imass-1)
   y2 = table_mass(imass)
 
   f11 = LOG10(table_yield(imet-1, imass-1))
   f21 = LOG10(table_yield(imet, imass-1))
   f22 = LOG10(table_yield(imet, imass))
   f12 = LOG10(table_yield(imet-1, imass))

   ! Coefficients of linear system.
   a1 = f11*x2*y2 / ((x1-x2)*(y1-y2)) + f12*x2*y1 / ((x1-x2)*(y2-y1)) + &
     &  f21*x1*y2 / ((x1-x2)*(y2-y1)) + f22*x1*y1 / ((x1-x2)*(y1-y2))
   a2 = f11*y2 / ((x1-x2)*(y2-y1)) + f12*y1 / ((x1-x2)*(y1-y2)) + &
     &  f21*y2 / ((x1-x2)*(y1-y2)) + f22*y1 / ((x1-x2)*(y2-y1))
   a3 = f11*x2 / ((x1-x2)*(y2-y1)) + f12*x2 / ((x1-x2)*(y1-y2)) + &
     &  f21*x1 / ((x1-x2)*(y1-y2)) + f22*x1 / ((x1-x2)*(y2-y1))
   a4 = f11 / ((x1-x2)*(y1-y2)) + f12 / ((x1-x2)*(y2-y1)) + &
     &  f21 / ((x1-x2)*(y2-y1)) + f22 / ((x1-x2)*(y1-y2))
 
   ! Compute bilinear interpolation.
   interp_yield = a1 + a2*met + a3*mass + a4*met*mass
   interp_yield = 10**interp_yield

end function interp_yield

!---------------------------------------
FUNCTION OBwind_Fe_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.2262d0
   real(dp) :: a1=-1.1031d0
   real(dp) :: OBwind_Fe_yield
 
   OBwind_Fe_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_Fe_yield

!---------------------------------------
FUNCTION OBwind_O_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=1.8280d0
   real(dp) :: a1=-0.8183d0
   real(dp) :: OBwind_O_yield
 
   OBwind_O_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_O_yield

!---------------------------------------
FUNCTION OBwind_N_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.0899d0
   real(dp) :: a1=-0.8936d0
   real(dp) :: OBwind_N_yield
 
   OBwind_N_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_N_yield

!---------------------------------------
FUNCTION OBwind_Mg_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.0709d0
   real(dp) :: a1=-1.5738d0
   real(dp) :: OBwind_Mg_yield

   OBwind_Mg_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_Mg_yield

!---------------------------------------
FUNCTION OBwind_Al_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.2067d0
   real(dp) :: a1=-2.4871d0
   real(dp) :: OBwind_Al_yield
 
   OBwind_Al_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_Al_yield

!---------------------------------------
FUNCTION OBwind_Si_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.1027d0
   real(dp) :: a1=-1.4978d0
   real(dp) :: OBwind_Si_yield
 
   OBwind_Si_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_Si_yield

!---------------------------------------
FUNCTION OBwind_Eu_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=2.2263d0
   real(dp) :: a1=-7.6284d0
   real(dp) :: OBwind_Eu_yield
 
   OBwind_Eu_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_Eu_yield

!---------------------------------------
FUNCTION OBwind_C_yield(met)
!---------------------------------------
   real(dp),intent(in) :: met
   real(dp) :: a0=1.9689d0
   real(dp) :: a1=-1.1375d0
   real(dp) :: OBwind_C_yield
 
   OBwind_C_yield = OB_wind_fit(met,a0,a1)
END FUNCTION OBwind_C_yield

function OB_wind_fit(met, a0, a1)
   real(dp),intent(in) :: met,a0,a1
   real(dp) :: OB_wind_fit
   ! ---------
   ! Fit returns IMF weighted mass fraction of the yield.
   ! The mass limits used for OB stars are 8-60 Msun and
   ! yields are taken as imf weighted averages with yield
   ! data from NuGrid (Ritter et al., 2018). 
   real(dp) :: logZ,logy

   ! Limit metallicity to range.
   if (met > 0.02)then
      logZ = -1.69897d0
   else if (met < 1.0d-4)then
      logZ = -4.d0
   else
      logZ = log10(met)
   endif

   logy = a0*logZ + a1
   OB_wind_fit = 10**logy

end function OB_wind_fit

end module metal_yields
