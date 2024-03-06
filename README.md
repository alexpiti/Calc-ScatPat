# calc-ScatPat
MATLAB implementation of the **Huygens-Fresnel principle (HFP)** for the calculation of the **scattering pattern** produced from the illumination of a rectangular metasurface.

This repositories includes the following three m-files:
1. The core-function ```calc_ScatPat.m``` 
2. A test-script (```test_script.m```) to showcase the execution
3. The ```plot_3D_Pattern.m``` function (can also be found in my [other repository](https://github.com/alexpiti/Plot-Pattern)), with modified parameters.

Here are the two figures produced by the execution of ```test_script.m```

![MS_config](https://user-images.githubusercontent.com/97299585/171837671-ff581032-ff32-41e3-9b42-edc133d4bc7f.png)

Fig.1 : Phase (left) and amplitude (right) profile of the reflection coefficient of each cell across the metasurface. Note that the reflection coefficient profile is assumed homogenized, i.e., any coupling between adjacent cells has been accounted for. 

![Pattern](https://user-images.githubusercontent.com/97299585/171837691-3dcbc20d-bfd5-4e19-845d-ae0f80162bdb.png)

Fig.2 : Scattering pattern produced from olbique (30 deg) illumination in the far-field. Left is in spherical coordinates, right is in cylindrical (also known as "U-V") coordinates. Note that the present HFP implementation is scalar, i.e., it works with only one polarization or, equivalently, it produces only the co-polarized scattering pattern.

The resulting scattering pattern can be post-processed to extract various metrics, e.g., directivity, HPBW, SLL, etc.

Extensions to this HFP implementation for the response in the transmissive (refractive) hemisphere or from non-uniformely arranged scatterers can be easily implemented. The most strong limitations for HFP validity are: 
* (i) The incident wavefront is spatially coherent, i.e., it has large curvature-over-lambda. Equivalently, it means that all (point) sources are far from the scatterer-array.
* (ii) The arrangement of the scatterers that make up the metasurface forms a convex, flat, or only mildly-concave aperture (so that there's no wave bouncing between scatterers).
* (iii) The obliquity is generally low, i.e., incidence/scattering is mostly from/towards the "broadside" of the aperture (angles up to 45-60 deg).

If you use these functions for research, **please cite my 2022 paper** ["Multi-functional metasurface architecture for amplitude, polarization and wavefront control"](https://doi.org/10.1103/PhysRevApplied.17.064060). These functions have been used also in other papers such as ["Multiwideband Terahertz Communications Via Tunable Graphene-Based Metasurfaces in 6G Networks"](https://arxiv.org/abs/2203.10298) and ["Scalability Analysis of Programmable Metasurfaces for Beam Steering"](https://arxiv.org/abs/2004.06917).

A similar physical-optics tool, **Fresnel-Kirchhoff diffraction (FKD)** from a finite-aperture metasurface, that can be used for imaging or for wireless path-loss calculation between two antennas with scattering off a metasurface, can be found in Section III.A-B of my **2023 paper** ["On the Mobility Effect in UAV-Mounted Absorbing Metasurfaces: A Theoretical and Experimental Study"](https://doi.org/10.1109/ACCESS.2023.3299379). A broader discussion/perspective on the physical/wave-optics usage in wireless-communications with metasurfaces can be found in Section III.C therein. I plan to add FKD to the present MATLAB toolpack, as it complements HFP.

Thanks for reading/visiting =)
