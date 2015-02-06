#Name: intersect_areas.py
# Purpose: Intersect of 2 polygons (feature class in geodatabase) to get 
# 3 values: area of polygon A, area of polygon B, area of intersection
# by Viacheslav Shalisko 2015

# Import system modules
import sys
import arcpy
from arcpy import env
 
try:
    # Set the workspace
    env.workspace = "C:/Users/Viacheslav/Documents/Projects_GIS/Grammitidaceae_MaxEnt/Pruebas_201410/GIS/Moranopteris_areas_Lambert.gdb"    
    
    # The input parameters for intersect
    # inFeature1 = "Moranopteris_achilleifolia_u03_b1km"
    # inFeature2 = "Moranopteris_gradata_u015_b1km"
    # intersectOutput = "_test_output1"

    inFeature1 = sys.argv[1]
    inFeature2 = sys.argv[2]
    intersectOutput = sys.argv[3]

    # Get area of first polygon
    rows_A = arcpy.SearchCursor(inFeature1)
    area_A = 0
    for row in rows_A:
        area_A = row.getValue("Shape_Area") + area_A
    # Print the intersect area
    print area_A

    # Get area of second polygon
    rows_B = arcpy.SearchCursor(inFeature2)
    area_B = 0
    for row in rows_B:
        area_B = row.getValue("Shape_Area") + area_B
    # Print the intersect area
    print area_B

    # Intersect and save result in the same workspace
    # [todo: check existence of layer before make intersect]
    inFeatures = [ inFeature1, inFeature2 ]
    arcpy.Intersect_analysis(inFeatures, intersectOutput, "ONLY_FID")

    # Get area of resulting intersect polygon
    rows_C = arcpy.SearchCursor(intersectOutput)
    area_C = 0
    for row in rows_C:
        area_C = row.getValue("Shape_Area") + area_C 
    # Print the intersect area
    print area_C
 
 
except Exception, e:
    # If an error occurred, print line number and error message
    import traceback, sys
    tb = sys.exc_info()[2]
    print "Exception in intersect_areas.py in line %i" % tb.tb_lineno
    print e.message