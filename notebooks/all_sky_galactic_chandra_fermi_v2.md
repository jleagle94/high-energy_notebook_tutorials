---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.17.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# This notebook will do the following:

* Queries metadata for all Chandra observations using HEASARC VO services.

* Extracts their sky coordinates.

* Plots them on an all-sky map using matplotlib Aitoff projection.

* Queries all sources along the Galactic plane from the Chandra Source Catalog and the Fermi Source Catalog

* Performs a simple cross-match of the sources, saved to a new table

* Assesses the properties of the new sample of sources that have a CXC and Fermi counterpart


# Dependencies

* astropy
* pyvo
* matplotlib
* healpy
* os
* astroquery
* numpy


## Step 1: Imports

```python
import pyvo as vo
import healpy
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.heasarc import Heasarc
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import numpy as np
from astropy.table import unique, vstack
from astropy.wcs import WCS
%matplotlib inline
plt.style.use(astropy_mpl_style)
```

## Step 2: Use HEASARC VO Service to access Chandra observation metdata

```python
heasarc = vo.dal.TAPService("https://heasarc.gsfc.nasa.gov/xamin/vo/tap")
```

## Step 3: Run a ADQL query to get all observation positions

```python
query = """
SELECT TOP 99999999 obsid, ra, dec, exposure, detector, time, name
FROM chanmaster
WHERE ra IS NOT NULL AND dec IS NOT NULL
"""


results = heasarc.search(query)
table = results.to_table()
table[:5]
```

## Step 4: Define the coordinates using SkyCoord

```python
coords = SkyCoord(ra=table['ra'],
                  dec=table['dec'])
```

## Step 5: Plot an all-sky map in Aitoff projection

```python
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection="aitoff")

mask = table['exposure'] != 0
exp_log = np.log10(table['exposure'][mask]/1000)

# Convert to radians and flip for Aitoff
ra_rad = coords.ra.wrap_at(180*u.deg).radian[mask]
dec_rad = coords.dec.radian[mask]

sc = ax.scatter(ra_rad, dec_rad, c=exp_log,s=1, alpha=0.5, cmap='plasma')
cb = plt.colorbar(sc, orientation='horizontal', pad=0.05)
cb.set_label("log10([exposure time]) (ks)",fontsize=12)
ax.grid(True)
plt.title("Chandra Observations (All-Sky)",fontsize=12,pad=20)
plt.show()
```

## Step 6: Find Chandra source locations and fluxes for all-sky using the CSC

### Find those that lie along the Galactic plane and plot the sky map.

```python
default_cols = ['name', 'ra', 'dec', 'b_flux_ap']
columns_str = ', '.join(default_cols)
# construct a query for all entries; use TOP with a large number greater than the server's 1e5 LIMIT
query = f"SELECT TOP 9999999 {columns_str} FROM csc"
sample = heasarc.search(query).to_table()
```

```python
print(sample)
```

There are a total of 407806 sources in the CSC. Now let's grab the relevant properties we need.

```python
query = """
SELECT TOP 9999999 name, ra, dec, b_flux_ap
FROM csc
WHERE b_flux_ap > 0
"""

sample = heasarc.search(query)
sample_table = sample.to_table()
sample_table[:5]
```

```python
print(len(sample_table))
print(f"{sample_table['b_flux_ap'].min():.2e}")
print(sample_table['b_flux_ap'].max())
print(np.sum(sample_table['b_flux_ap'])/len(sample_table))
```

We can check the sample properties in the cell above to make sure we are getting the same output each time. See also [GitHub issue #3387 in Astroquery](https://github.com/astropy/astroquery/issues/3387).

In principle, one could try:

``sample = heasarc.query_region(
SkyCoord(0, 0, unit='deg',frame='galactic'),
spatial='all-sky',
catalog='csc',
maxrec=None
)``

Once the following bugs are resolved:
* fields parameter errors out every time.
* maxrec messes with the query, sending a random sample that matches the criteria for a maximum of 1e5 sources. 

Next we mask sources that meet the Galactic plane criteria: |b|< 10&deg;.

```python
# Create SkyCoord from the sample table
coords = SkyCoord(ra=sample_table['ra'], dec=sample_table['dec'], unit='deg')
gal_coords = coords.galactic

# Apply Galactic latitude filter: |b| <= 10°
mask = np.abs(gal_coords.b.deg) <= 10

# Filter coordinates and data
filtered_coords = gal_coords[mask]
filtered_flux = sample_table['b_flux_ap'][mask]
filtered_names = sample_table['name'][mask]

print(f"Total sources in catalog: {len(coords)}")
print(f"Number of filtered sources: {len(filtered_coords)}")
```

Next, we will wrap the longitude coordinates for Aitoff projection. We did this above for Chandra observations, too.

```python
print("Longitude range (unwrapped):", filtered_coords.l.min().deg, "to", filtered_coords.l.max().deg)
print("Longitude range (wrapped):", filtered_coords.l.wrap_at(180*u.deg).min().deg, "to", filtered_coords.l.wrap_at(180*u.deg).max().deg)
```

```python
flux_log = np.log10(filtered_flux)

l_rad = filtered_coords.l.wrap_at(180*u.deg).radian
b_rad = filtered_coords.b.radian

plt.figure(figsize=(10, 6))
ax = plt.subplot(111, projection='aitoff')

sc = ax.scatter(l_rad, b_rad, c=flux_log, cmap='plasma', s=5, alpha=0.7)
cb = plt.colorbar(sc, orientation='horizontal', pad=0.05)
cb.set_label("log10([b_flux_ap]) (erg/cm2/s)",fontsize=12)

ax.grid(True)
plt.title(r"Chandra Source Catalog Sources within $\pm$ 10$\circ$ of Galactic Plane",fontsize=12,pad=20)
plt.show()
```

We can plot a histogram of source count by longitude, too. Naturally the Galactic Center sees the largest number of sources. 

```python
l_deg = filtered_coords.l.deg

plt.hist(l_deg, bins=36, range=(0, 360))
plt.xlabel('Galactic Longitude (deg)',fontsize=10)
plt.ylabel('Number of sources',fontsize=10)
plt.title('Distribution of sources near Galactic plane by longitude',fontsize=12)
plt.show()
```

The code below plots the sources along the Galactic plane again, but using the Astropy WCS package so we can display only the Galactic plane.

```python
# 1. Build a WCS header manually for Galactic coordinates
wcs = WCS(naxis=2)
wcs.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']  # Galactic longitude/latitude
wcs.wcs.crval = [0, 0]  # Center of the plot in l and b
wcs.wcs.crpix = [900, 50]  # Reference pixel
wcs.wcs.cdelt = [-0.2, 0.2]  # Degrees per pixel
wcs.wcs.cunit = ['deg', 'deg']

# 2. Prepare your data (Galactic longitudes/latitudes in degrees)
l_deg = filtered_coords.l.wrap_at(180 * u.deg).deg
b_deg = filtered_coords.b.deg
flux_log = np.log10(filtered_flux)

# 3. Convert world coords (deg) to pixel coords
x_pix, y_pix = wcs.wcs_world2pix(l_deg, b_deg, 0)

# 4. Plot
fig = plt.figure(figsize=(9, 4))
ax = fig.add_subplot(111, projection=wcs)

sc = ax.scatter(x_pix, y_pix, c=flux_log, cmap='cool', s=5, alpha=0.7)

cb = plt.colorbar(sc, orientation='vertical', pad=0.05)
cb.set_label('log10(Broadband Flux)')

# 5. Set axis limits in pixel coords to match lon/lat range
ax.set_xlim(0, 1850)  # 360 deg / 0.2 deg per pixel = 1800 pixels
ax.set_ylim(0, 100)   # 20 deg / 0.2 deg per pixel = 100 pixels

# 6. Axis labels, ticks, grid (in world coordinates)
ax.coords.grid(True, color='grey', ls='-')
ax.coords[0].set_axislabel('Galactic Longitude (deg)')
ax.coords[1].set_axislabel('Galactic Latitude (deg)')
ax.coords[0].set_ticks(spacing=20 * u.deg)
ax.coords[1].set_ticks(spacing=5 * u.deg)

plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15)
plt.show()
```

## Step 7: Perform a cross-match of the Fermi and CSC Catalogs using ADQL

What the query will filter:
* $ |b| < 10\,$&deg; 
* CSC flux > 0
* A CSC source must fall within the error circle of the Fermi source defined by ``circle('ICRS', fsc.ra, fsc.dec, fsc.semi_minor_axis_95)``

```python
query = """
    select fsc.name,fsc.lii,fsc.bii, fsc.semi_minor_axis_95, fsc.energy_flux, fsc.source_type, csc.name, csc.lii, csc.bii, csc.b_flux_ap 
    from fermi3fgl as fsc, csc 
    where ( fsc.bii between -10 and 10 and csc.b_flux_ap > 0)
    and contains(
        point('ICRS', csc.ra, csc.dec),
        circle('ICRS', fsc.ra, fsc.dec, fsc.semi_minor_axis_95)
      ) = 1"""
result = heasarc.search(query)
```

```python
table = result.to_table()
new_names = [
    "fermi_name", "fermi_l", "fermi_b", "fermi_r95", "fermi_flux", "fermi_class",
    "csc_name", "csc_l", "csc_b", "csc_flux"]
table.rename_columns(table.colnames,new_names)
print(table)
```

There are 4,155 sources that have a CSC counterpart within the Fermi error circle. Note: This includes several of the same Fermi source overlapping with multiple CSC sources. We will filter this further below. 

```python
csc_coords = SkyCoord(l=table['csc_l'], b=table['csc_b'], frame='galactic')
fermi_coords = SkyCoord(l=table['fermi_l'], b=table['fermi_b'], frame='galactic')

# Compute angular separation
sep2d = csc_coords.separation(fermi_coords)

# Add the information as a new column to the table
table['sep_deg'] = sep2d.deg
```

### Plot a sky map of the matched table, plotting each source only once even if it appears more than once in the table just above. 

```python
# Unique CSC and Fermi source rows
unique_csc = unique(table, keys='csc_name')
unique_fermi = unique(table, keys='fermi_name')

# CSC SkyCoords
csc_coords = SkyCoord(l=unique_csc['csc_l'], b=unique_csc['csc_b'], unit='deg', frame='galactic')
csc_lons = csc_coords.l.wrap_at(180 * u.deg).radian
csc_lats = csc_coords.b.radian

# Fermi SkyCoords
fermi_coords = SkyCoord(l=unique_fermi['fermi_l'], b=unique_fermi['fermi_b'], unit='deg', frame='galactic')
fermi_lons = fermi_coords.l.wrap_at(180 * u.deg).radian
fermi_lats = fermi_coords.b.radian

# Fermi R95 sizes
marker_scale = 100
fermi_r95_deg = unique_fermi['fermi_r95']
fermi_sizes = (fermi_r95_deg * marker_scale)**2

# Plot
plt.figure(figsize=(12, 6))
ax = plt.subplot(111, projection='aitoff')
ax.grid(True)

ax.scatter(csc_lons, csc_lats, s=5, color='blue', alpha=0.6, label='CSC source')
ax.scatter(fermi_lons, fermi_lats, s=fermi_sizes, edgecolor='red', facecolor='none', alpha=0.7)
ax.scatter(0, 0, s=5, edgecolor='red', facecolor='none', alpha=0.7, label='Fermi R95')

ax.set_title("CSC sources with 4FGL counterparts (Galactic Coordinates)", pad=20, fontsize=12)
ax.legend(loc='best')
plt.subplots_adjust(bottom=0.1, top=0.95)
plt.show()
```

The code below plots the matched table along the Galactic plane again, but using the Astropy WCS package so we can display only the Galactic plane.

```python
# Log flux ratio
flux_ratio = np.log10(table['csc_flux'] / table['fermi_flux'])

# Build a WCS header manually for Galactic coordinates
wcs = WCS(naxis=2)
wcs.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']  # Galactic longitude/latitude
wcs.wcs.crval = [0, 0]  # Center of the plot in l and b
wcs.wcs.crpix = [900, 50]  # Reference pixel
wcs.wcs.cdelt = [-0.2, 0.2]  # Degrees per pixel
wcs.wcs.cunit = ['deg', 'deg']

# CSC SkyCoords
csc_coords = SkyCoord(l=table['csc_l'], b=table['csc_b'], unit='deg', frame='galactic')
csc_lons = csc_coords.l.wrap_at(180 * u.deg).deg
csc_lats = csc_coords.b.deg

# Fermi SkyCoords
fermi_coords = SkyCoord(l=unique_fermi['fermi_l'], b=unique_fermi['fermi_b'], unit='deg', frame='galactic')
fermi_lons = fermi_coords.l.wrap_at(180 * u.deg).deg
fermi_lats = fermi_coords.b.deg

# Fermi R95 sizes
marker_scale = 100
fermi_r95_deg = unique_fermi['fermi_r95']
fermi_sizes = (fermi_r95_deg * marker_scale)**2

# Convert to pixel coords
csc_x_pix, csc_y_pix = wcs.wcs_world2pix(csc_coords.l.wrap_at(180 * u.deg).deg, csc_coords.b.deg, 0)
fermi_x_pix, fermi_y_pix = wcs.wcs_world2pix(fermi_coords.l.wrap_at(180 * u.deg).deg, fermi_coords.b.deg, 0)

# Plot
fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111, projection=wcs)

sc = ax.scatter(csc_x_pix, csc_y_pix, c=flux_ratio, cmap='coolwarm', s=8, alpha=0.8, label='CSC')

ax.scatter(fermi_x_pix, fermi_y_pix, s=fermi_sizes, edgecolor='red', facecolor='none', alpha=0.7)
ax.scatter(0, 0, s=5, edgecolor='red', facecolor='none', alpha=0.7, label='Fermi R95')

cb = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.05)
cb.set_label(r'Log$_{10}$(X-ray to Gamma-ray Flux Ratio)',fontsize=10)

ax.set_xlim(0, 1800)
ax.set_ylim(0, 100)

ax.coords.grid(True, color='gray', ls='--')
ax.coords[0].set_axislabel('Galactic Longitude (deg)',fontsize=10)
ax.coords[1].set_axislabel('Galactic Latitude (deg)',fontsize=10)
ax.coords[0].set_ticks(spacing=30 * u.deg)
ax.coords[1].set_ticks(spacing=5 * u.deg)

ax.set_title("Associated Sources by Flux Ratio", pad=10,fontsize=12)
ax.legend(loc='upper right')

plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)
plt.show()
```

## Step 8: Build a final table consolidating repeating entries

Get rid of repeating Fermi entries. However, we want to retain the information that several CSC sources can overlap with the same Fermi source, so we build a nested table.

```python
# Step 1: Unique Fermi sources
unique_pairs = unique(table, keys=['csc_name', 'fermi_name'])
unique_fermi_names = list(set(unique_pairs['fermi_name']))

# Step 2: Prepare row containers
rows = []

for fname in unique_fermi_names:
    # Filter for this Fermi source
    mask = unique_pairs['fermi_name'] == fname
    group = unique_pairs[mask]

    # Use first row for Fermi info
    f_row = group[0]

    # Extract CSC info as Python lists
    csc_names = list(group['csc_name'])
    csc_fluxes = list(group['csc_flux'])
    csc_lons = list(group['csc_l'])
    csc_b_lats = list(group['csc_b'])

    n_matches = len(csc_names)

    # Append row
    rows.append((
        fname,
        f_row['fermi_l'],
        f_row['fermi_b'],
        f_row['fermi_flux'],
        f_row['fermi_r95'],
        n_matches,
        csc_names,
        csc_fluxes,
        csc_lons,
        csc_b_lats
    ))

# Step 3: Build final table
grouped_display_table = Table(rows=rows, names=[
    'fermi_name', 'fermi_l', 'fermi_b', 'fermi_flux', 'fermi_r95',
    'num_csc_matches', 'csc_names', 'csc_fluxes', 'csc_glons', 'csc_glats'
])

print(len(rows))
```

There are 184 Fermi sources with CSC counterparts. Save the final table as a CSV. 

```python
export_table = grouped_display_table.copy()
export_table.sort(['fermi_l', 'fermi_b'])

# Convert all list-like columns to strings
for col in export_table.colnames:
    if isinstance(export_table[col][0], (list, tuple)):
        export_table[col] = [', '.join(str(x) for x in row) for row in export_table[col]]

# Write to CSV
export_table.write("grouped_csc_fermi_pairs.csv", format="csv", overwrite=True)
```

## Step 9: Let's plot some of the properties of the final table

How many CSC matches per Fermi source as a function of longitude. Again we see the Galactic Center have the largest number of matches. 

```python
# Count CSC matches per Fermi source using a dict
fermi_counts = {}
for fname in table['fermi_name']:
    fermi_counts[fname] = fermi_counts.get(fname, 0) + 1

fermi_lons = []
num_csc_matches = []

for fname, count in fermi_counts.items():
    # Find fermi longitude in unique_fermi
    idx = (unique_fermi['fermi_name'] == fname).nonzero()[0]
    if len(idx) == 0:
        continue
    lon = unique_fermi['fermi_l'][idx[0]]
    fermi_lons.append(lon)
    num_csc_matches.append(count)

plt.figure(figsize=(8,5))
plt.bar(fermi_lons, num_csc_matches, 100*unique_fermi['fermi_r95'], alpha=0.7,lw=2)
plt.xlabel('Fermi Galactic Longitude (deg)',fontsize=10)
plt.ylabel('Number of CSC Matches',fontsize=10)
plt.title('Number of CSC Matches vs Fermi Longitude', fontsize=12)
plt.ylim(0,100)
plt.grid(True)
plt.show()
```

What is more interesting are properties like the X- to gamma-ray flux ratios, which can help identify a specific source class or subclass. 

```python
log_flux_ratio = np.log10(table['csc_flux'] / table['fermi_flux'])

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)  # 1 row, 2 columns

# First plot
axes[0].scatter(table['fermi_l'], log_flux_ratio, alpha=0.6, color='C1')
axes[0].set_xlabel('Fermi Galactic Longitude (deg)', fontsize=10)
axes[0].set_ylabel('Log10(CSC Flux / Fermi Flux)', fontsize=10)
axes[0].grid(True)

# Second plot
axes[1].scatter(table['fermi_b'], log_flux_ratio, alpha=0.6, color='C2')
axes[1].set_xlabel('Fermi Galactic Latitude (deg)', fontsize=10)
axes[1].grid(True)

plt.tight_layout()
plt.show()
```

Finally, we can plot the flux ratio as a function of source class (from the Fermi catalog). 

```python
fermi_classes_clean = [' '.join(c.split()) for c in table['fermi_class']]

class_counts = {cls: fermi_classes_clean.count(cls) for cls in set(fermi_classes_clean)}

# Define desired order explicitly - galactic to extragalactic
desired_order = ['PSR', 'psr', 'PWN', 'pwn', 'SNR', 'snr', 'SPP', 'spp','HMB', 'hmb', '',
                'glc', 'FSRQ', 'fsrq', 'sey', 'BCU', 'bcu', 'BLL', 'nlsy1', 'gal', 'sbg', 
                 'BIN', 'RDG', 'rdg', 'bll']

nonzero_classes = [cls for cls in desired_order if class_counts.get(cls, 0) > 0]
# Extract unique classes and map to numbers for plotting
class_to_num = {cls: i for i, cls in enumerate(nonzero_classes)}

x_vals = np.array([class_to_num[c] for c in fermi_classes_clean if c in class_to_num])
jitter = np.random.normal(0, 0.1, size=len(x_vals))
x_jittered = x_vals + jitter

plt.figure(figsize=(12,6))
plt.scatter(x_jittered, log_flux_ratio, alpha=0.6)
plt.xticks(range(len(nonzero_classes)), nonzero_classes, rotation=45)
plt.xlabel('Fermi Class', fontsize=10)
plt.ylabel('Log10(CSC Flux / Fermi Flux)', fontsize=10)
plt.grid(True)
plt.tight_layout()
plt.show()
```

With this cross-matched set you can do lots of things:

* Do brighter gamma-ray sources have more X-ray counterparts? Make a scatter plot of CSC matches as a function of Fermi flux.

* Assess association reliability. Here we only do a superficial cross-match using Fermi's 95% semiminor axis of each source. Plot the difference between Fermi and CSC coordinates as a histogram. 

* What regions of the sky show a clustering of sources? Plot a sky map where the point size is proportional to the number of CSC matches. 

* How many sources that lie along the Galactic plane are classified as extragalactic source classes? From this, what is the probability any given source is Galactic in nature? 

* Create a new table of the subsample (maybe prioritizing one or more similar source classes) to explore further trends, etc. 
