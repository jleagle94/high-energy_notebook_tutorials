# high-energy_notebook_tutorials
Notebook tutorials for open source data accessibility and usability. 

Click any project below to view:

- [Project 1: Chandra and Fermi-LAT Galactic Source Population View using 4FGL](heasarc_notebooks/all_sky_galactic_chandra_fermi.md)
  - Project 1 pulls CSC catalog data from heasarc Xamin using pyvo's TAP service. Fermi 4FGL data is extracted directly from FSSC website for latest catalog data.
  - [Project 1.5: Same thing but using same database throughout (restricted to 3FGL for now)](heasarc_notebooks/all_sky_galactic_chandra_fermi_v2.md)
  - Project 1.5 pulls both CSC and Fermi catalog data from heasarc Xamin using pyvo's TAP service, which only has up to 3FGL Fermi catalog data.
  - Benefit to using ADQL for both CSC and Fermi datasets is a more straight forward cross-match using Fermi uncertainty sizes directly in the query search.
- [Project 2: Expands Project 1 to include IRSA and MAST archives](fornax_notebooks/csc-fermi-gaia-2mass_cross-match/cxc-fermi-gaia-2mass_cross-match.md)
  - Pulls both CSC and Fermi (4FGL-DR4) catalog data from Xamin and cross-matches with IRSA-2MASS and MAST-GAIA datasets

Other material:
- [HEASARC Cloud Tutorial](heasarc_notebooks/heasarc_cloud_tutorial.md): Working version of the [HEASARC cloud tutorial](https://heasarc.gsfc.nasa.gov/docs/archive/cloud/heasarc_cloud_tutorial.html)
  - Quick fix: pyvo image search for heasarc chandra data must specify the format type
  - Temporary work around to: access_url in image search results is broken (there is a // just before obsid directory in path)
- Any notebooks not listed here are working (but not really guaranteed to work) draft versions of future notebooks. 
