#!/bin/sh

################################################################################
#
# NAME:         potential analysis for roof greening in Bonn (NRW)
#
# AUTHOR:    	  Meike Reimann (2019)
#
#
# PURPOSE:      analysis of LIDAR and other required data for rooftop greening
#               analysis using Open Data NRW
#
# REQUIREMENTS: PDAL (http://www.pdal.io), standard system tools
#               (basename, cat, ...), GRASS GIS (https://grass.osgeo.org/) ,
#               QGIS (https://www.qgis.org/de/site/)
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#
# PROCEDURE: 1. extracting areas for rooftop greening analysis by eliminating
#               wrong HU and overlapping trees as well as buffering buildings to
#               supplement small house foot prints
# PROCEDURE: 2. slope analysis and classification of roof greening potential
# PROCEDURE: 3. create mosaic for classification of roof greening potential
#
# comment: the script is based on script "data_download_import.sh".
#          The results of the analyses are mainly based on the year and
#          quality of aerial
#          survey of the digital surface model
################################################################################


#./Analyse_Bonn_forschleife
for i in `seq -s ' ' 1 51`; do
    echo "Processing district: $i"

    # extract area for analyses. In this example a district in Bonn
    g.mapset -c location=Dachbegruenung_BN mapset=ortsteil_${i} --o
    v.extract input=ortsteile_reproj@ortsteil out=ortsteil_${i} type=area \
        where="cat=${i}"
    g.region vector=ortsteil_${i}@ortsteil_${i} res=0.1 \
        align=dop3_mosaik@DOPBonn
    r.mask vector=ortsteil_${i}@ortsteil_${i}
    #connect ALKIS and OSM Buildings to complete  possible missing buildings
    v.overlay ainput=ALKIS@Bonn_buildings atype=area \
        binput=ortsteil_${i}@ortsteil_${i} btype=area operator=and \
        output=ALKIS_${i}
    v.overlay ainput=OSMBuildings@Bonn_buildings atype=area \
        binput=ortsteil_${i}@ortsteil_${i} btype=area \
        operator=and output=OSM_${i}
    v.overlay ainput=ALKIS_${i}@ortsteil_${i} atype=area \
        binput=OSM_${i}@ortsteil_${i} btype=area operator=or \
        output=ALKIS_OSM
    #remove monument protection areas
    v.overlay ainput=ALKIS_OSM@ortsteil_${i} atype=area \
        binput=denkmal_reproj@denkmal btype=area operator=not \
        out=ALKIS_OSM_Denkmal
    #clean small areas
    v.clean tool=rmarea th=10 input=ALKIS_OSM_Denkmal@ortsteil_${i} \
        out=ALKIS_OSM_clean
    #set buffer around the buildings to complete too small building outlines
    v.to.rast input=ALKIS_OSM_clean@ortsteil_${i} \
        type=area out=ALKIS_OSM_clean use=cat
    r.buffer input=ALKIS_OSM_clean@ortsteil_${i} output=buffer2m \
        distances=2 units=meters
    r.mask -r
    r.mask raster=buffer2m

    r.grow.distance input=ALKIS_OSM_clean@ortsteil_${i} \
        distance=distance value=value
    r.to.vect in=value@ortsteil_${i} out=value type=area
    v.generalize in=value@ortsteil_${i} type=area out=valuesmooth \
        method=snakes thres=0.1
    v.clean tool=rmarea th=10 input=valuesmooth@ortsteil_${i} \
        out=buildings_${i}
    v.to.rast input=buildings_${i}@ortsteil_${i} type=area \
        out=buildings_${i} use=cat

    r.mask -r
    r.mask raster=buildings_${i}@ortsteil_${i}

    i.vi red=dop3_mosaik@DOPBonn nir=dop4_mosaik@DOPBonn viname=ndvi \
        output=ndvi_${i}

    #transform floating point to integer to minimize size of the NDVI Raster
    r.mapcalc expression="ndvi_${i}int100 = int( ndvi_${i}*100)"

    #extract areas where ndvi <=0.18 -> removes trees
    r.mask -r
    r.mask raster=ndvi_${i}int100 where="value<=18"

    r.to.vect in=dgm_mosaik@DGMBonn \
        out=dgm_mosaik_${i}@ortsteil_${i} type=area
    r.to.vect in=dom_mosaik@DOMBonn \
        out=dom_mosaik_${i}@ortsteil_${i} type=area

    #analysis of high differences to detect wrong building outlines:

    v.db.addcolumn map=dom_mosaik_${i}@ortsteil_${i} columns="nr int"
    v.db.addcolumn map=dom_mosaik_${i}@ortsteil_${i} \
        columns="dgm double precision"
    v.db.addcolumn map=dom_mosaik_${i}@ortsteil_${i} \
        columns="diff double precision"

    v.what.vect map=dom_mosaik_${i}@ortsteil_${i} column=nr \
        query_map=buildings_${i}@ortsteil_${i} query_column=cat
    v.what.vect map=dom_mosaik_${i}@ortsteil_${i} \
        column=dgm query_map=dgm_mosaik_${i}@ortsteil_${i} \
        query_column=value
    v.db.update map=dom_mosaik_${i}@ortsteil_${i} col=diff \
        qcol="value-dgm"
    #extract building <2 m and <10sqm
    v.extract input=dom_mosaik_${i} out=buildings2m type=area \
        where="diff>=2"
    v.to.rast in=buildings2m out=buildings2m type=area use=attr \
        attribute_column=nr
    r.to.vect in=buildings2m output=buildings2m_nr type=area column=value

    v.db.addcolumn map=buildings2m_nr columns="area double precision"
    v.to.db map=buildings2m_nr option=area columns=area units=meters
    v.extract input=buildings2m_nr out=buildings10qm type=area where="area>=10"
    v.to.rast in=buildings10qm out=buildings10qm type=area use=attr \
        attribute_column=value

    r.mask -r
    r.mask raster=buildings10qm@ortsteil_${i}

    #set resolution to DOM resolution for slope analysis
    g.region res=1

    r.slope.aspect elevation=dom_mosaik@DOMBonn slope=slope_e aspect=aspect_e \
        pcurvature=pcurv_e tcurvature=tcurv_e -e

    #set resolution to 0.1 for interpolation
    g.region res=0.1
    r.resamp.interp input=slope_e@ortsteil_${i} out=slope_e_resamp \
        method=bilinear

    r.mapcalc expression="slope_e_resamp_int = int(slope_e_resamp*10)"

    r.reclass input=slope_e_resamp_int@ortsteil_${i} \
        out=slope_classification rules=/scratch/Meike/rules.txt
    r.to.vect in=slope_classification@ortsteil_${i} \
        out=roof_classification type=area -s
    v.extract input=roof_classification@ortsteil_${i} out=slope45 \
        type=area where="value<5"

    # connect buildings with building numbers -> a_cat is the generated
    # building number of buildings_${i}

    v.overlay ainput=buildings_${i}@ortsteil_${i} atype=area \
        binput=slope45@ortsteil_${i} btype=area olayer=1,0,0 \
        snap=0.00000001 out=slope45_id operator=and

    # preparation for classification of slope values
    v.db.renamecolumn map=slope45_id@ortsteil_${i} column=b_value,Klasse

    # class1 (0 - 2 degrees):
    v.extract input=slope45_id@ortsteil_${i} \
        out=Klasse1_${i}@ortsteil_${i} type=area where="Klasse=1"
    v.reclass input=Klasse1_${i}@ortsteil_${i} \
        out=RG1${i}@ortsteil_${i} column=a_cat
    v.db.addtable map=RG1${i}@ortsteil_${i} column="areacl DOUBLE"
    v.to.db map=RG1${i}@ortsteil_${i} option=area unit=meters \
        column=areacl
    #2D area multiplication with smallest roof pitch-factor for roof surface
    v.db.addcolumn map=RG1${i}@ortsteil_${i} column="areasl DOUBLE"
    v.db.update map=RG1${i}@ortsteil_${i} column=areasl \
        query_column="areacl * 1"

    v.db.addcolumn map=RG1${i}@ortsteil_${i} column="Klasse INT"
    v.db.update map=RG1${i}@ortsteil_${i} column=Klasse value=1

    #class2 (>2 - 15 degrees):
    v.extract input=slope45_id@ortsteil_${i} \
        out=Klasse2_${i}@ortsteil_${i} type=area where="Klasse=2"
    v.reclass input=Klasse2_${i}@ortsteil_${i} \
        out=RG2${i}@ortsteil_${i} column=a_cat
    v.db.addtable map=RG2${i}@ortsteil_${i} column="areacl DOUBLE"
    v.to.db map=RG2${i}@ortsteil_${i} \
        option=area unit=meters column=areacl
    #2D area multiplication with smallest roof pitch-factor for roof surface
    v.db.addcolumn map=RG2${i}@ortsteil_${i} column="areasl DOUBLE"
    v.db.update map=RG2${i}@ortsteil_${i} column=areasl \
        query_column="areacl * 1.0006"

    v.db.addcolumn map=RG2${i}@ortsteil_${i} column="Klasse INT"
    v.db.update map=RG2${i}@ortsteil_${i} column=Klasse value=2

    #class3 (>15 - 30 degrees):
    v.extract input=slope45_id@ortsteil_${i} \
        out=Klasse3_${i}@ortsteil_${i} type=area where="Klasse=3"
    v.reclass input=Klasse3_${i}@ortsteil_${i} \
        output=RG3${i}@ortsteil_${i} column=a_cat
    v.db.addtable map=RG3${i}@ortsteil_${i} column="areacl DOUBLE"
    v.to.db map=RG3${i}@ortsteil_${i} option=area unit=meters \
        column=areacl
    #2D area multiplication with smallest roof pitch-factor for roof surface
    v.db.addcolumn map=RG3${i}@ortsteil_${i} column="areasl DOUBLE"
    v.db.update map=RG3${i}@ortsteil_${i} column=areasl \
        query_column="areacl * 1.0353"

    v.db.addcolumn map=RG3${i}@ortsteil_${i} column="Klasse INT"
    v.db.update map=RG3${i}@ortsteil_${i} column=Klasse value=3

    #class4 (>30 - 45 degrees):
    v.extract input=slope45_id@ortsteil_${i} \
        out=Klasse4_${i}@ortsteil_${i} type=area where="Klasse=4"
    v.reclass input=Klasse4_${i}@ortsteil_${i} \
        output=RG4${i}@ortsteil_${i} column=a_cat
    v.db.addtable map=RG4${i}@ortsteil_${i} \
        column="areacl DOUBLE"
    v.to.db map=RG4${i}@ortsteil_${i} option=area unit=meters \
        column=areacl
    #2D area multiplication with smallest roof pitch-factor for roof surface
    v.db.addcolumn map=RG4${i}@ortsteil_${i} column="areasl DOUBLE"
    v.db.update map=RG4${i}@ortsteil_${i} column=areasl \
        query_column="areacl * 1.1547"

    v.db.addcolumn map=RG4${i}@ortsteil_${i} column="Klasse INT"
    v.db.update map=RG4${i}@ortsteil_${i} column=Klasse value=4

    g.mapset mapset=roofgreening
    g.copy vector=RG1${i}@ortsteil_${i},RG1${i}
    g.copy vector=RG2${i}@ortsteil_${i},RG2${i}
    g.copy vector=RG3${i}@ortsteil_${i},RG3${i}
    g.copy vector=RG4${i}@ortsteil_${i},RG4${i}

done



#patches all vector files of Klasse1
MAPSRG1=`g.list type=vector separator=comma pat="RG1*"`
g.region raster=dop4_mosaik@DOPBonn
v.patch in=${MAPSRG1} out=RG1_mosaik -e


#patches all vector files of Klasse2
MAPSRG2=`g.list type=vector separator=comma pat="RG2*"`
g.region raster=dop4_mosaik@DOPBonn
v.patch in=${MAPSRG2} out=RG2_mosaik -e

#patches all vector files of Klasse3
MAPSRG3=`g.list type=vector separator=comma pat="RG3*"`
g.region raster=dop4_mosaik@DOPBonn
v.patch in=${MAPSRG3} out=RG3_mosaik -e

#patches all vector files of Klasse4
MAPSRG4=`g.list type=vector separator=comma pat="RG4*"`
g.region raster=dop4_mosaik@DOPBonn
v.patch in=${MAPSRG4} out=RG4_mosaik -e

#copy patches to new mapset for better overview
g.mapset -c mapset=patches_classes
g.copy vector=RG1_mosaik@roofgreening,RG1_mosaik
g.copy vector=RG2_mosaik@roofgreening,RG2_mosaik
g.copy vector=RG3_mosaik@roofgreening,RG3_mosaik
g.copy vector=RG4_mosaik@roofgreening,RG4_mosaik

#remove duplicate area centroids
v.clean in=RG1_mosaik out=RG1_mosaik_clean tool=rmdac
v.clean in=RG2_mosaik out=RG2_mosaik_clean tool=rmdac
v.clean in=RG3_mosaik out=RG3_mosaik_clean tool=rmdac
v.clean in=RG4_mosaik out=RG4_mosaik_clean tool=rmdac

#mosaic for Bonn with all classes
v.patch in=RG1_mosaik_clean,RG2_mosaik_clean,RG3_mosaik_clean,RG4_mosaik_clean \
    out=green_roof_potential -e
v.clean in=green_roof_potential out=green_roof_potential_clean type=area \
    tool=rmdac


#connect green_roof_potential with "ortsteile" vector
v.overlay ainput=green_roof_potential_clean binput=ortsteile_reproj@ortsteil \
    atype=area btype=area operator=and out=green_roof_potential_ort

#rename attribute column, if it is to long for the later export as ESRI shape
# file

v.db.renamecolumn map=green_roof_potential_ort column=a_areacl,areacl
v.db.renamecolumn map=green_roof_potential_ort column=a_areasl,areasl
v.db.renamecolumn map=green_roof_potential_ort column=a_Klasse,Klasse
v.db.renamecolumn map=green_roof_potential_ort column=b_ortsteil,Ort
v.db.renamecolumn map=green_roof_potential_ort column=b_bezirk,Bezirk
v.db.renamecolumn map=green_roof_potential_ort column=b_ortsteilnm,Ortnm
v.db.renamecolumn map=green_roof_potential_ort column=b_bezirknm,Bezirknm
