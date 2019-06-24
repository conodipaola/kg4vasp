# *__KG4VASP__*

__version__: v1.0 stable (stand alone)

## Acknowledgements:

I would like to thank:

* kgec developers [^1] for their inspiring work and effort to deliver a nice tool for *QE*
* for the mathematical derivation of various quantities from the generalised Kubo-Greenwood formula [^2]
* *Quantum espresso* [^3]
* *VASP* [^4]


## The code
This code is meant to be as a VASP post-processing tool for people that want to calculate transport properties.

The quantities calculated at the moment are:

1. Electrical conductivity/resistivity

2. Seebeck coefficient

3. Electron thermal conductivity

4. Possibly, a flavour of the relaxation time $\tau$ can be obtained by fitting $\sigma (\omega) = \sigma_0 / 1+\omega^2 \tau^2$ at $\omega \rightarrow 0$ (Drude Model)

    ####   The INPUT file
    *kg4vasp* input file is called INPUT and needs that  the following FLAGS have to be addressed:

    * SYSTEM (character): usually name of the system in study
    * NBANDS (integer): number of bands
    * FILE_NABLA (character): name of file where the DFT code stores nabla matrix elements
    * NKPOINTS (integer): number of K-POINTS
    * NP_FREQ (integer): number of frequency between min and max
    * DELTAE (real): energy threshold for degeneracy between two contiguous energies
    * MAX_FREQ (real): max value for frequency in eV
    * MIN_FREQ (real): min value for frequency in eV
    * GAUS_SIGMA (real): spreading in eV for delta function that is approximated by 1 gaussian function  
    * VOLUME (real): cell volume in $Ang ^3$
    * DEGAUSS (real): electronic temperature $K_b \cdot  T$ (degauss in quantum espresso and sigma in VASP if the smearing has been set up on ) in eV
    * INTRA_BAND (boolean): **TRUE**=in the calculation intraband and degerate energies are included; **FALSE**=other way round
    * FERMI_ENERGY (real): fermi energy from SCF (STATIC DFT CALCULATION) in eV
    * FROM_VASP (boolean): **TRUE**=post-processing of VASP calculations; **FALSE**=post-processing of QE calculations
    * OCC_FROM_CODE (boolean):
    * N_ELECTRON (integer): total number of electrons in the system
    * NABLA_BIN (boolean, default=**.FALSE.**): **TRUE**: the nabla file is in binary format; **FALSE** the nabla file is in ASCII formatted file.









## The repository
The structure is the following:

1. src: source code in fortran90
2. examples: test drive for kg4vasp
    a) Simple VASP calculations
    b) from kgec code simple QE/kgec calculations
    c) comparison between them and validation
3) PATCHES: necessary kgec v2.1 and VASP 5.4.4 to run kg4vasp and  it with QE output as well.


[^1]: https://doi.org/10.1016/j.cpc.2017.08.008 L.Calderín, V.V.Karasiev2, S.B.Trickey, Computer Physics Communications Volume 221, December 2017, Pages 118-142
[^2]: https://doi.org/10.1016/j.hedp.2009.06.004
S.Mazevet,M.Torrent,V.Recoules,F.Jollet, High Energy Density Physics Volume 6, Issue 1, January 2010, Pages 84-88

[^3]: http://dx.doi.org/10.1088/0953-8984/21/39/395502 P. Giannozzi, S. Baroni, N. Bonini, M. Calandra, R. Car, C. Cavazzoni, D. Ceresoli, G. L. Chiarotti, M. Cococcioni, I. Dabo, A. Dal Corso, S. Fabris, G. Fratesi, S. de Gironcoli, R. Gebauer, U. Gerstmann, C. Gougoussis, A. Kokalj, M. Lazzeri, L. Martin-Samos, N. Marzari, F. Mauri, R. Mazzarello, S. Paolini, A. Pasquarello, L. Paulatto, C. Sbraccia, S. Scandolo, G. Sclauzero, A. P. Seitsonen, A. Smogunov, P. Umari, R. M. Wentzcovitch, J.Phys.:Condens.Matter 21, 395502 (2009)

[^4]: G. Kresse and J. Hafner, Phys. Rev. B 47 , 558 (1993); ibid. 49 , 14 251 (1994).
