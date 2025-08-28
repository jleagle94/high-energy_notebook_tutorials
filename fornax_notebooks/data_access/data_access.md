# This notebook will do the following

* Access data in the AWS S3 bucket for all archives (HEASARC, IRSA, and MAST)
* Derived from: 
    * <a href="https://github.com/HEASARC/sciserver_cookbooks/blob/main/data-find-download.md">data-find-download-sciserver</a>
    * <a href="https://github.com/HEASARC/sciserver_cookbooks/blob/main/data-access.md">data-access-sciserver</a>
    
        * introduce cloud-specific options like:
        ```
        from astropy.io import fits
        hdul = fits.open(s3_uri, use_fsspec=True, fsspec_kwargs={"anon": True})
        to stream the file directly from S3. 
        ```
        and
        ```
        import heasoftpy as hsp
        hsp.fdump(
        infile="https://nasa-heasarc.s3.amazonaws.com/chandra/data/byobsid/5/4475/primary/acisf04475N004_full_img2.fits.gz",
        outfile='STDOUT',columns=' ',rows='1',prhead=False,more=False)
        ```
        which has the command line version
        ```
        > fdump infile='https://nasa-heasarc.s3.amazonaws.com/chandra/data/byobsid/5/4475/primary/acisf04475N004_full_img2.fits.gz' outfile=STDOUT columns=' ' rows=1 prhead=no more=no
        ```
    * <a href="https://heasarc.gsfc.nasa.gov/docs/archive/cloud.html">heasarc-cloud</a>
    * <a href="https://github.com/nasa-fornax/fornax-s3-subsets/blob/main/notebooks/astropy-s3-subsetting-demo.ipynb">astropy-s3-subsetting</a>


# Let's explore some of the tools and methods.

## Data in the cloud:

### AWS S3 cloud service
   <a href="https://heasarc.gsfc.nasa.gov/docs/archive/cloud.html">HEASARC Data in the Cloud</a>
   
   <a href="https://irsa.ipac.caltech.edu/cloud_access/">IRSA Data in the cloud</a>
   
   <a href="https://outerspace.stsci.edu/display/MASTDOCS/Public+AWS+Data"> MAST Data in the Cloud</a>
   
##  Python Tools: 
### 1. s3fs, fsspec, and boto3

### 2. pyvo

## Methods: 
### 1. ADQL query search with pyvo

### 2. TAP, SIA, SSA, SCS, and SLAP with pyvo

For more detailed information on the various tools and methods, see the ``data_access_advanced`` notebook. 


## Step 1: Imports

```python
import sys
import os
import s3fs
import pyvo as vo
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

import boto3
from botocore import UNSIGNED
from botocore.client import Config
s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))
```


```python
#use only s3fs, pyvo, and boto3!

#access mast via s3://stpubdata. list all directories here
#plot some file somewhere? 

#access irsa s3://nasa-irsa. list all directories here. 
#plot some file somewhere? 

#access heasarc s3:nasa-heasarc. list all directories here. 
#plot some file somewhere? 
```

```python
s3_client.list_objects_v2(Bucket="nasa-heasarc", Prefix="fermi/data/lat/", Delimiter="/")
```




    {'ResponseMetadata': {'RequestId': 'T9HS34B12MTX03V4',
      'HostId': 'Q3cfmmkQoPpTu912b7sBMSi4cXMiwlDrkZts/1+6TGPmSKe22ouhHziEdesdYfysUedJSt9yHFI=',
      'HTTPStatusCode': 200,
      'HTTPHeaders': {'x-amz-id-2': 'Q3cfmmkQoPpTu912b7sBMSi4cXMiwlDrkZts/1+6TGPmSKe22ouhHziEdesdYfysUedJSt9yHFI=',
       'x-amz-request-id': 'T9HS34B12MTX03V4',
       'date': 'Tue, 19 Aug 2025 19:47:39 GMT',
       'x-amz-bucket-region': 'us-east-1',
       'content-type': 'application/xml',
       'transfer-encoding': 'chunked',
       'server': 'AmazonS3'},
      'RetryAttempts': 0},
     'IsTruncated': False,
     'Name': 'nasa-heasarc',
     'Prefix': 'fermi/data/lat/',
     'Delimiter': '/',
     'MaxKeys': 1000,
     'CommonPrefixes': [{'Prefix': 'fermi/data/lat/mission/'},
      {'Prefix': 'fermi/data/lat/weekly/'}],
     'EncodingType': 'url',
     'KeyCount': 2}




```python
bucket="nasa-heasarc"
prefix="fermi/data/lat/"
```


```python
def list_s3_tree(bucket, prefix, indent=0):
    """Recursively list subdirectories under a given S3 prefix."""
    paginator = s3_client.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/")
    
    for page in pages:
        for p in page.get("CommonPrefixes", []):
            print(" " * indent + "- " + p["Prefix"])
            list_s3_tree(bucket, p["Prefix"], indent + 2)

print(f"Full directory tree under s3://{bucket}/fermi/data/lat/\n")
list_s3_tree(bucket, "fermi/data/lat/")
```

    Full directory tree under s3://nasa-heasarc/fermi/data/lat/
    
    - fermi/data/lat/mission/
      - fermi/data/lat/mission/spacecraft/
    - fermi/data/lat/weekly/
      - fermi/data/lat/weekly/1s_spacecraft/
      - fermi/data/lat/weekly/diffuse/
      - fermi/data/lat/weekly/extended/
      - fermi/data/lat/weekly/photon/
      - fermi/data/lat/weekly/spacecraft/



```python
paginator = s3_client.get_paginator("list_objects_v2")
pages = paginator.paginate(Bucket=bucket, Prefix=prefix)

print(f"Files under s3://{bucket}/{prefix}\n")
for page in pages:
    for obj in page.get("Contents", []):
        print("-", obj["Key"])
```

```python

```
