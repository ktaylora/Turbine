#!/bin/bash

while read x; do
  echo $x | awk '{ print "gdal_rasterize -burn 1 -sql \"SELECT * FROM Transmission_Lines_TX_Albers WHERE (STATUS = 1) AND (VOLTAGE ="$1")\" -ot Float32 -tr 0.0002777778 0.0002777778 Transmission_Lines_TX_Albers.shp raster_"$1".tif" }' | /bin/bash
done < voltage_classes.txt
