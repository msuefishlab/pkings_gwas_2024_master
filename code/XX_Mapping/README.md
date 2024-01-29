# Load Data
Need to add description of where data originally came from, etc.

# Watershed and Catchment Analysis

We're going to use the ArcHydro Tools Pro toolchain to perform these operations (see https://www.youtube.com/watch?v=06q7_90tk54)

## Step 1: Convert Small Rivers to Raster 
1. Set ndende_small_rivers id field to 1
2. ndende_small_rivers > Polyline to Raster
    - value field "id"
    - output coordinate ystem WGS_1984_UTM_Zone_32S
    - environments > extent "RawDEM"
    - cellsize = 25 (2* DEM cellsize)
3. Reclassify Raster Small Rivers
    - set "NODATA"=0

## Step 2: Create "Reconditioned" DEM

1. river_raster + RawDEM > DEM Reconditioning

## Step 3: Fill Sinks in DEM

AgreeDEM.tif > Fill Sinks

## Step 4: Flow Direction

Fil.tif > Flow Direction

## Step 5: Flow Accumulation

Fdr.tif > Flow Accumulation

## Step 6: Stream Definition

Fac.tif > Stream Definition
    - Area SqKm to define stream: 1 (6400 cells)

## Step 7: Stream Segmentation

Str.tif Fdr.tif > Stream Segmentation

## Step 8: Catchment Defintiion

Fdr.tif StrLnk.tif > Catchment Grid Deliniation

## Step 8: Catchment Polygon Processing

Cat.tif > Catchment Polygon Processing

## Step 9: Drainage Line Processing

StrLink.tif Fld.tif > Drainage Line Processing

## Step 10: Stream Order

Drainage Line > Assign River Order
    - Input River Order Field: order_
    - RIver Order Type: Strahler 

## Step 11: Adjust Symbology for Drainage Lines
    - Graduated Symbols
    -Field order_
    - Symbol: water
    - CLasses: 6

## Step 12: Adjoint Catchment Processing

Drainage LIne Catchment > Adjoint Catchment Processing

## Step 13: Clean Up Symbology for Catchments

## Step 14: Create Outlets Layer and Mark Outlets for Streams of Interest

Catalog View > New Shape File
    - Feature Class Name: Outlets
    - Coordinate System: WGS_1984_UTM_Zone_32S

# Process ALOS PALSAR Data
Following this data recipie:
https://asf.alaska.edu/how-to/data-recipes/how-to-map-regional-inundation-with-spaceborne-l-band-sar-using-arcgis/

## Organize Data
1. Download from Vertex (Alaska Satellite Facility) based on ROI.  Filter to RTC datasets.
2. Extract ZIP images, and obtain the HH image from each granuale, place in folder
3. Load these images into ArcGIS

## Convert to dB
Raster Calculator 
    - db=10.*Log10*("RTC_HH_granule")

Repeat for Each

## Filter Images to Reduce Speckle
Speckle 
    - Enhanced Lee
    - Keep other paramters default

## Copy Raster to Output Folder
Raster Dataset > 
    - Output RRaster Dataset Field:
    ./input-data/XX_Mapping/SourceData/db_converted_filtered_despeckled/

## Open A New Project (ALOS Processing)

## Create Mosaics

Notes on Filename:

AA BBBBB CCC DDDDD RTE

AA- Mission
BB - Orbit No
CC - Beam
DD - Frame
E - Processing Level (1=High Rez, 2=LowRez)

Mosaic Files by Date:
13519: August 7, 2008
10835: Feb 5, 2008
09070: Oct 7, 2007
08399: Aug 22, 2007
05715: Feb 19, 2007
05467: Feb 2, 2007
021857 June 22, 2006

Mosaic To New Raster > 
(always input the frames in order, 120,130,140)
- Output Location: ALOS_Mosacics_final
- Spatial Reference for Raster: WGS_1984_UTM_Zone_32S
- Pixel Type 32 bit float
- Cellsize 12.5
- Number of Bands 1
- Mosaic Operator Last
- Mosaic Colormap Mode First

## Import to Integrated Map

Convert Rasters to "Classified" with Breakpoints at -11, -4
Set "center" color to "nocolor" and low color to blue and high color to yellow

Batch Apply Symbology from Layer