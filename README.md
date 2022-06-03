# calc-ScatPat
MATLAB implementation of the **Huygens-Fresnel principle** for the calculation of the **scattering pattern** produced from the illumination of a rectangular metasurface.

This repositories includes the following three m-files:
1. The core-function ```calc_ScatPat.m``` 
2. A test-script (```test_script.m```) to showcase the execution
3. The ```plot_3D_Pattern.m``` function (can also be found in my [other repository](https://github.com/alexpiti/Plot-Pattern)), with modified parameters.

Here are the two figures produced by the execution of ```test_script.m```

![MS_config](https://user-images.githubusercontent.com/97299585/171837671-ff581032-ff32-41e3-9b42-edc133d4bc7f.png)

Fig.1 : Phase (left) and amplitude (right) profile of the reflection coefficient of each cell across the metasurface.

![Pattern](https://user-images.githubusercontent.com/97299585/171837691-3dcbc20d-bfd5-4e19-845d-ae0f80162bdb.png)

Fig.2 : Scattering pattern produced from olbique (30 deg) illumination. Left is in spherical coordinates, right is in cylindrical coordinates.

The resulting scattering pattern can be post-processed to extract various metrics, e.g., directivity, HPBW, SLL, etc.

If you use these functions for research, **please cite my most recent paper** ["Multi-functional metasurface architecture for amplitude, polarization and wavefront control"](https://arxiv.org/abs/2204.03962). These functions have been used also in other papers such as ["Multiwideband Terahertz Communications Via Tunable Graphene-Based Metasurfaces in 6G Networks"](https://arxiv.org/abs/2203.10298) and ["Scalability Analysis of Programmable Metasurfaces for Beam Steering"](https://arxiv.org/abs/2004.06917).
