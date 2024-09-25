# modis_cloud_radiative_kernels

[![DOI](https://zenodo.org/badge/835898737.svg)](https://zenodo.org/doi/10.5281/zenodo.13839355)

Code for decomposing the shortwave effective radiative forcing from aerosol-cloud interactions (SW ERFaci) from liquid clouds into components associated with the Twomey effect and LWP and CF adjustments is provided along with the associated SW cloud radiative kernel. The provided environment.yml file should enable the creation of a conda environment that allows the notebook to be executed.

## References 
- Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: [Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part I: Cloud Radiative Kernels](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00248.1). J. Climate, 25, 3715-3735. 
    doi:10.1175/JCLI-D-11-00248.1.
- Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: [Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part II: Attribution to Changes in Cloud Amount, Altitude, and Optical Depth](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00249.1). 
    J. Climate, 25, 3736-3754. doi:10.1175/JCLI-D-11-00249.1.
- Zelinka, M.D., S.A. Klein, K.E. Taylor, T. Andrews, M.J. Webb, J.M. Gregory, and P.M. Forster, 2013: 
    [Contributions of Different Cloud Types to Feedbacks and Rapid Adjustments in CMIP5](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00555.1). 
    J. Climate, 26, 5007-5027. doi:10.1175/JCLI-D-12-00555.1.
- Scott, R. C., Myers, T. A., Norris, J. R., Zelinka, M. D., Klein, S. A., Sun, M., and Doelling, D. R., 2020: [Observed Sensitivity of Low-Cloud Radia-
tive Effects to Meteorological Perturbations over the Global Oceans](https://journals.ametsoc.org/view/journals/clim/33/18/jcliD191028.xml). J. Climate, 33, 7717 – 7734.
- Zelinka et al. (2022): [Evaluating climate models’ cloud feedbacks against expert judgement](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JD035198), _J. Geophys. Res._, 127, e2021JD035198, doi:10.1029/2021JD035198.
- Wall, C. J., Storelvmo, T., A. Possner, 2023: [Global observations of aerosol indirect effects from marine liquid clouds](https://acp.copernicus.org/articles/23/13125/2023/). Atmospheric Chemistry and Physics, 23, 13 125–13 141, doi:10.5194/acp-23-13125-2023.
- Zelinka, M., Chao, L.-W., Myers, T., Qin, Y., and Klein, S., preprint: [Technical Note: Recommendations for Diagnosing Cloud Feedbacks and Rapid
Cloud Adjustments Using Cloud Radiative Kernels](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-2782/). EGUsphere, 1–27.

## Input

The associated code requires the following inputs, including new effective radius x liquid water path cloud fraction joint histograms from the MODIS satellite simulator:

| Frequency | Name | Description | Unit | File Format |
|-----------|------|-------------|------|-------------|
| monthly mean | CLMODIS_LWPR | MODIS simulator liquid-cloud fraction histograms | % | nc |
| monthly mean | CLTMODIS | MODIS simulator total cloud fraction | % | nc |
| monthly mean | FSDSC | Clearsky downwelling solar flux at surface | W/m^2 | nc |
| monthly mean | FSNSC | Clearsky net solar flux at surface | W/m^2 | nc |
| monthly mean | TS     | surface temperature | K     | nc            |
| monthly mean | SWkernel | SW cloud radiative kernel | W/m^2/% | nc |

The SW cloud radiative kernel is available to download at https://github.com/brandonduran/modis_cloud_radiative_kernels/tree/main/data

- ensmean_SW_kernel.nc: SW cloud radiative kernel developed using zonal mean temperature and humidity profiles averaged across control runs of five CMIP6-era climate models as input to the RRMTG radiation code. These are best for diagnosing feedbacks / forcing relative to a modeled pre-industrial climate state. Please refer to Wall et al. (2023) and Duran et al. (in prep) for details.
