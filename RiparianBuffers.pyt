import math
import os
import multiprocessing
import time

globTransectId = 0

def fieldExists(featureclass, fieldname):
    fieldList = arcpy.ListFields(featureclass, fieldname)
    fieldCount = len(fieldList)
    if (fieldCount == 1):
        return True
    else:
        return False

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Riparian Buffers"
        self.alias = "RiparianBuffers"

        # List of tool classes associated with this toolbox
        #self.tools = [LinearTransect, TransectWindows, CalculateForestedBuffer, TrimTransectWindows]
        self.tools = [ForestedBuffer]

def LinearTransect_execute(in_feature,length,dist,units,out_feature,out_riv):
    try:
        LineDissolve = "LineDissolve.shp"
        LineSplit = "LineSplit.shp"
        Azline1="Azline1.shp"
        Azline2="Azline2.shp"
        Azline="Azline.shp"        
        #Azline_Dissolve="Azline_Dissolve.shp"

        #1.Take the input and dissolve lines 
        arcpy.env.outputMFlag = "Disabled"
        arcpy.env.outputZFlag = "Disabled"
        arcpy.Dissolve_management(in_feature, LineDissolve,"","","SINGLE_PART")
        
        #2.Run splitlines to have a version with a linebreak at each selected interval
        #arcpy.AddMessage("Splitting input")
        splitline(LineDissolve, LineSplit, dist)
             
        #3.Create parallel lines
        arcpy.AddField_management(LineSplit,"TransectID","DOUBLE")
        if (fieldExists(LineSplit,"OBJECTID")):
            idField = "!OBJECTID!"
        elif (fieldExists(LineSplit, "FID")):
            idField = "!FID!"   
        else:
            arcpy.AddError("Cannot find id field in {0}".format(out_file));        
        
        global globTransectId
        
        CodeBlock_TransectID="""
rec={0}
def autoIncrement(): 
 global rec 
 pStart = 1  
 pInterval = 1 
 if (rec == 0):  
  rec = pStart  
 else:  
  rec += pInterval  
 return rec""".format(globTransectId)
 
        

        arcpy.CalculateField_management (LineSplit, "TransectID", "autoIncrement()", "PYTHON_9.3", CodeBlock_TransectID)
        if (out_riv):
            arcpy.CopyFeatures_management(LineSplit,out_riv)

        FieldsNames=["Direction", "Azimuth", "X_mid", "Y_mid", "AziLine_1", "AziLine_2", "Distance"]
        for fn in FieldsNames:
            arcpy.AddField_management (LineSplit, fn, "DOUBLE")
            
        CodeBlock_Direction="""def GetAzimuthPolyline(shape):
         num = (shape.lastpoint.x - shape.firstpoint.x)
         denom = (shape.lastpoint.y - shape.firstpoint.y)
         if (denom == 0):
          return 0
         radian = math.atan(num/denom)
         degrees = radian * 180 / math.pi
         return degrees"""
         
        CodeBlock_Azimuth="""def Azimuth(direction):
         if direction < 0:
          azimuth = direction + 360
          return azimuth
         else:
          return direction"""
          
        CodeBlock_NULLS="""def findNulls(fieldValue):
            if fieldValue is None:
                return 0
            elif fieldValue is not None:
                return fieldValue"""                                       

        arcpy.CalculateField_management (LineSplit, "Direction", "GetAzimuthPolyline(!Shape!)", "PYTHON_9.3", CodeBlock_Direction)
        #arcpy.AddMessage('Calc findNulls(direction)')
        arcpy.CalculateField_management (LineSplit, "Direction", "findNulls(!Direction!)", "PYTHON_9.3", CodeBlock_NULLS)
        arcpy.CalculateField_management (LineSplit, "Azimuth", "Azimuth(!Direction!)", "PYTHON_9.3", CodeBlock_Azimuth)
        arcpy.CalculateField_management (LineSplit, "X_mid", "!Shape!.positionAlongLine(0.5,True).firstPoint.X", "PYTHON_9.3")
        arcpy.CalculateField_management (LineSplit, "Y_mid", "!Shape!.positionAlongLine(0.5,True).firstPoint.Y", "PYTHON_9.3")
        CodeBlock_AziLine1="""def Azline1(azimuth):
         az1 = azimuth + 90
         if az1 > 360:
          az1-=360
          return az1
         else:
          return az1"""
        CodeBlock_AziLine2="""def Azline2(azimuth):
         az2 = azimuth - 90
         if az2 < 0:
          az2+=360
          return az2
         else:
          return az2"""
        arcpy.CalculateField_management (LineSplit, "AziLine_1", "Azline1(!Azimuth!)", "PYTHON_9.3", CodeBlock_AziLine1)
        
        arcpy.CalculateField_management (LineSplit, "AziLine_2", "Azline2(!Azimuth!)", "PYTHON_9.3", CodeBlock_AziLine2) 

        arcpy.CalculateField_management (LineSplit, "Distance", length, "PYTHON_9.3")
        
        #Generate Azline1 and Azline2
        spatial_reference=arcpy.Describe(in_feature).spatialReference
        arcpy.BearingDistanceToLine_management (LineSplit, Azline1, "X_mid", "Y_mid", "Distance", units, "AziLine_1", "DEGREES", "GEODESIC", "TransectID", spatial_reference)
        
        arcpy.BearingDistanceToLine_management (LineSplit, Azline2, "X_mid", "Y_mid", "Distance", units, "AziLine_2", "DEGREES", "GEODESIC", "TransectID", spatial_reference)
        
        #Create Azline and append Azline1 and Azline2
        arcpy.AddField_management(Azline1, "Side", "SHORT")
        arcpy.AddField_management(Azline2, "Side", "SHORT")
        arcpy.CalculateField_management(Azline1, "Side", 1, "PYTHON_9.3")
        arcpy.CalculateField_management(Azline2, "Side", 2, "PYTHON_9.3")

        arcpy.CreateFeatureclass_management(arcpy.env.workspace, Azline, "POLYLINE", "", "", "", spatial_reference)

        arcpy.AddField_management (Azline, "TransectID", "DOUBLE") 
        arcpy.AddField_management(Azline, "Side", "SHORT")    
        arcpy.Append_management ([Azline1, Azline2], Azline, "NO_TEST")
        arcpy.CopyFeatures_management(Azline,out_feature)
        
        '''#Dissolve Azline
        #arcpy.Dissolve_management (Azline, Azline_Dissolve,"TransectID", "", "SINGLE_PART")
        
        Azline_Dissolve = Azline #TEST! Don't change now because I want them to remain split at the middle
        
        #Add Fields to Azline_Dissolve
        FieldsNames2=["TransectID","x_start", "y_start", "x_end", "y_end"]
        for fn2 in FieldsNames2:
            arcpy.AddField_management (Azline_Dissolve, fn2, "DOUBLE")
            
        #Calculate Azline_Dissolve fields
        arcpy.CalculateField_management (Azline_Dissolve, "x_start", "!Shape!.positionAlongLine(0,True).firstPoint.X", "PYTHON_9.3") 
        arcpy.CalculateField_management (Azline_Dissolve, "y_start", "!Shape!.positionAlongLine(0,True).firstPoint.Y", "PYTHON_9.3")
        arcpy.CalculateField_management (Azline_Dissolve, "x_end", "!Shape!.positionAlongLine(1,True).firstPoint.X", "PYTHON_9.3")
        arcpy.CalculateField_management (Azline_Dissolve, "y_end", "!Shape!.positionAlongLine(1,True).firstPoint.Y", "PYTHON_9.3")
        
        #Generate output file
        arcpy.XYToLine_management (Azline_Dissolve, out_feature,"x_start", "y_start", "x_end","y_end", "", idField[], spatial_reference)
        '''

        #Clean up files
        arcpy.Delete_management(LineDissolve)
        arcpy.Delete_management(LineSplit)
        arcpy.Delete_management(Azline1)
        arcpy.Delete_management(Azline2)
        arcpy.Delete_management(Azline)
        #arcpy.Delete_management(Azline_Dissolve)

    #clean up in case of failure
    except Exception as e:
        arcpy.AddError(e.message)
        if arcpy.Exists(LineDissolve):
            arcpy.Delete_management(LineDissolve)
        if arcpy.Exists(LineSplit):
            arcpy.Delete_management(LineSplit)
        if arcpy.Exists(Azline1):
            arcpy.Delete_management(Azline1)
        if arcpy.Exists(Azline2):
            arcpy.Delete_management(Azline2)
        if arcpy.Exists(Azline):
            arcpy.Delete_management(Azline)
        
    return

#Source: http://nodedangles.wordpress.com/2011/05/01/quick-dirty-arcpy-batch-splitting-polylines-to-a-specific-length/
def splitline (inFC,FCName,alongDist):
    arcpy.env.outputMFlag = "Disabled"
    arcpy.env.outputZFlag = "Disabled"
    OutDir = arcpy.env.workspace
    outFCName = FCName
    outFC = OutDir+"/"+outFCName
    
    def distPoint(p1, p2):
        calc1 = p1.X - p2.X
        calc2 = p1.Y - p2.Y

        return math.sqrt((calc1**2)+(calc2**2))

    def midpoint(prevpoint,nextpoint,targetDist,totalDist):
        newX = prevpoint.X + ((nextpoint.X - prevpoint.X) * (targetDist/totalDist))
        newY = prevpoint.Y + ((nextpoint.Y - prevpoint.Y) * (targetDist/totalDist))
        return arcpy.Point(newX, newY)

    def splitShape(feat,splitDist):
        # Count the number of points in the current multipart feature
        #
        partcount = feat.partCount
        partnum = 0
        # Enter while loop for each part in the feature (if a singlepart feature
        # this will occur only once)
        #
        lineArray = arcpy.Array()

        while partnum < partcount:
              # Print the part number
              #
              #arcpy.AddMessage("Part {0}: ".format(str(partnum)));
              part = feat.getPart(partnum)
              #arcpy.AddMessage(part.count)

              totalDist = 0

              pnt = part.next()
              pntcount = 0

              prevpoint = None
              shapelist = []

              # Enter while loop for each vertex
              #
              while pnt:

                    if not (prevpoint is None):
                        thisDist = distPoint(prevpoint,pnt)
                        maxAdditionalDist = splitDist - totalDist

                        #arcpy.AddMessage("{0}, {1}, {2}".format(thisDist, totalDist, maxAdditionalDist))

                        if (totalDist+thisDist)> splitDist:
                              while(totalDist+thisDist) > splitDist:
                                    maxAdditionalDist = splitDist - totalDist
                                    #arcpy.AddMessage("{0}, {1}, {2}".format(thisDist, totalDist, maxAdditionalDist))
                                    newpoint = midpoint(prevpoint,pnt,maxAdditionalDist,thisDist)
                                    lineArray.add(newpoint)
                                    shapelist.append(lineArray)

                                    lineArray = arcpy.Array()
                                    lineArray.add(newpoint)
                                    prevpoint = newpoint
                                    thisDist = distPoint(prevpoint,pnt)
                                    totalDist = 0

                              lineArray.add(pnt)
                              totalDist+=thisDist
                        else:
                              totalDist+=thisDist
                              lineArray.add(pnt)
                              #shapelist.append(lineArray)
                    else:
                        lineArray.add(pnt)
                        totalDist = 0
                    prevpoint = pnt                
                    pntcount += 1

                    pnt = part.next()

                    # If pnt is null, either the part is finished or there is an
                    #   interior ring
                    #
                    if not pnt:
                        pnt = part.next()
                        #if pnt:
                              #print "Interior Ring:"
              partnum += 1

        if (lineArray.count > 1):
              shapelist.append(lineArray)
        return shapelist

    #arcpy.AddMessage("in splitline")
    if arcpy.Exists(outFC):
        arcpy.Delete_management(outFC)

    arcpy.Copy_management(inFC,outFC)
    
    #origDesc = arcpy.Describe(inFC)
    #sR = origDesc.spatialReference

    #revDesc = arcpy.Describe(outFC)
    #revDesc.ShapeFieldName

    deleterows = arcpy.UpdateCursor(outFC)
    for iDRow in deleterows:       
         deleterows.deleteRow(iDRow)

    try:
        del iDRow
        del deleterows
    except:
        #print "Couldn't delete updateCursor "
        pass

    inputRows = arcpy.SearchCursor(inFC)
    outputRows = arcpy.InsertCursor(outFC)
    fields = arcpy.ListFields(inFC)

    numRecords = int(arcpy.GetCount_management(inFC).getOutput(0))
    OnePercentThreshold = numRecords // 100

    #printit(numRecords)
    #arcpy.AddMessage("numRecords: {0}".format(numRecords))

    iCounter = 0
    iCounter2 = 0

    for iInRow in inputRows:
        inGeom = iInRow.shape
        iCounter+=1
        iCounter2+=1    
        if (iCounter2 > (OnePercentThreshold+0)):
              #arcpy.AddMessage("Processing Record "+str(iCounter) + " of "+ str(numRecords))
              iCounter2=0
    
        if (inGeom.length > alongDist):
              shapeList = splitShape(iInRow.shape,alongDist)

              for itmp in shapeList:
                    newRow = outputRows.newRow()
                    for ifield in fields:
                        if (ifield.editable):
                              newRow.setValue(ifield.name,iInRow.getValue(ifield.name))
                    newRow.shape = itmp
                    outputRows.insertRow(newRow)
        else:
              outputRows.insertRow(iInRow)

    del inputRows
    del outputRows


###END SPLIT LINE CODE ###

class ForestedBuffer(object):
    def __init__(self):
        self.label = "Forested Buffer"
        self.description = "Splits an input hydrology into 10m intervals. For each 10m portion, a polygon on the left/right side is added and clipped to the width of the input canopy cover map. The calculated length and presence of each side of the buffer is added is a field to the split hydrology feature output."
        self.canRunInBackground = False
        
    def getParameterInfo(self):
        out_folder = arcpy.Parameter(displayName="Intermediate Output Folder (Reuses files from here if available)", name="out_folder", datatype="DEFolder", parameterType="Required", direction="Input")
        out_folder.value = arcpy.env.scratchFolder

        in_hydro = arcpy.Parameter(displayName="Input Hydrological Feature (UTM Projection Required)", name="in_hydro", datatype="GPFeatureLayer", parameterType="Required", direction="Input")
        in_hydro.filter.list = ["Polyline"]

        #Need to check for integer type, with 1 = cover.        
        in_cover = arcpy.Parameter(displayName="Canopy Cover Map (UTM Projection Required)", name="in_cover", datatype="GPRasterLayer", parameterType="Required", direction="Input")
        
        transect_length = arcpy.Parameter(displayName="Buffer Length (m)", name="transect_length", datatype="GPLong", parameterType="Required", direction="Input")        
        transect_length.value = 100
        
        transect_width = arcpy.Parameter(displayName="Buffer Slice Width (m)", name="transect_width", datatype="GPLong", parameterType="Required", direction="Input")        
        transect_width.value = 10

        tolerance = arcpy.Parameter(displayName="Tolerated Distance from Hydro Feature (m)", name="tolerance", datatype="GPLong", parameterType="Required", direction="Input")        
        tolerance.value = 5
        
        out_hydro = arcpy.Parameter(displayName = "Output Hydrological Feature", name="out_hydro", datatype="GPFeatureLayer", parameterType="Required", direction="Output")
        out_hydro.value = ''.join([arcpy.env.scratchFolder,"\\frb_hydro.shp"])
        
        out_frb = arcpy.Parameter(displayName = "Output Buffer Feature", name="out_frb", datatype="GPFeatureLayer", parameterType="Required", direction="Output")
        out_frb.value = ''.join([arcpy.env.scratchFolder,"\\frb_swath.shp"])
        
        isDel = arcpy.Parameter(displayName = "Delete intermediate files", name="isDel", datatype="GPBoolean", parameterType="Required", direction="Input")
        isDel.value = 'false'
        
        reuseFiles = arcpy.Parameter(displayName = "Continue with existing files", name="reuseFiles", datatype="GPBoolean", parameterType="Required", direction="Input")
        reuseFiles.value = 'true'
        return [out_folder, in_hydro, in_cover, transect_length, transect_width, tolerance, out_hydro, out_frb, isDel, reuseFiles]
    
    def isLicensed(self):
        #Should check for spatial analyst and license for flat end buffers
        return True
    
    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
    def updateMessages(self, parameters):
        return
    
    def execute(self,parameters,messages):
        
        out_folder = parameters[0].valueAsText
        in_hydro = parameters[1].valueAsText
        in_cover = parameters[2].valueAsText
        transect_length = parameters[3].valueAsText
        transect_width = parameters[4].valueAsText
        tolerance = parameters[5].valueAsText
        out_hydro = parameters[6].valueAsText
        out_frb = parameters[7].valueAsText
        isDel = (parameters[8].valueAsText == 'true')
        reuseFiles = (parameters[9].valueAsText == 'true')
        
        #our workspace
        arcpy.env.overwriteOutput = True
        

        if not (os.path.isdir(''.join([out_folder,"\\frb_data"]))):
            arcpy.CreateFolder_management(out_folder,"frb_data")
        out_folder = ''.join([out_folder,"\\frb_data"])
                
        zonespath = ''.join([out_folder, "\\zones"])
        rivpath =''.join([out_folder, "\\riv_split"])
        transpath = ''.join([out_folder, "\\trans"])
        bufpath = ''.join([out_folder, "\\buf"])
        frbpath = ''.join([out_folder, "\\frb"])
        temppath = ''.join([out_folder, "\\temp"])

        if not (os.path.isdir(temppath)):
            arcpy.CreateFolder_management(out_folder,"temp")
            
        #1. Create fishnet and split rivers into 25 squares
        fishnet_fc = ''.join([temppath,"\\fishnet_fc.shp"])
        desc = arcpy.Describe(in_hydro)
        origin_coord = str(desc.extent.lowerLeft)
        y_axis_coord = ' '.join([str(desc.extent.XMin),str(desc.extent.YMax + 10)]).strip(" NaN")
        corner_coord = str(desc.extent.upperRight).strip(" NaN")
        # ~10km x 10km zones
        num_rows = math.ceil((desc.extent.XMax-desc.extent.XMin)/10000)
        num_cols = math.ceil((desc.extent.YMax-desc.extent.YMin)/10000)
        
        if not (os.path.isfile(fishnet_fc)):
            arcpy.AddMessage("Creating fishnet")
            #Errors here usually means you had spaces in the file name!!
            arcpy.CreateFishnet_management(fishnet_fc, origin_coord, y_axis_coord, 0, 0, num_rows, num_cols, corner_coord, "NO_LABELS", in_hydro, "POLYGON")
            arcpy.AddField_management(fishnet_fc,"Name","TEXT","#","#","#","#","NULLABLE","NON_REQUIRED","#")
            arcpy.CalculateField_management(fishnet_fc,"Name",""""zone_"+str(!FID!)""","PYTHON_9.3","#")
        
            if not os.path.isdir(zonespath):
                arcpy.CreateFolder_management(out_folder,"zones")
                arcpy.Split_analysis(in_hydro,fishnet_fc,"Name",os.path.join(out_folder,"zones"),"#")
        
        
        #2. Create rough frb perim
    
        #2.1 buffer buffer_length round of rivers
        buffer_fc = ''.join([temppath,"\\hydro_buffer.shp"])
        if not (os.path.isfile(buffer_fc)):
            arcpy.AddMessage("Create hydro_buffer")
            arcpy.Buffer_analysis(in_hydro,buffer_fc,transect_length+" Meters","FULL","ROUND","ALL","#")
        #2.2 extract mask to buffer
        cover_subset_fc = ''.join([temppath,"\\cover_subset.tif"])
        cover_subset_lyr = 'cover_subset'
        if not (os.path.isfile(cover_subset_fc)):
            arcpy.AddMessage("Extract mask to buffer")
            arcpy.CheckOutExtension('Spatial')
            arcpy.gp.ExtractByMask_sa(in_cover,buffer_fc,cover_subset_fc)
        
        #2.3. raster to polygon
        cover_poly_fc = ''.join([temppath, "\\cover_poly.shp"])
        if not (os.path.isfile(cover_poly_fc)):    
            arcpy.AddMessage("Raster to Polygon")  
            arcpy.MakeRasterLayer_management(cover_subset_fc, cover_subset_lyr) #needs to be in layer for raster-to-polygon tool
            arcpy.RasterToPolygon_conversion(cover_subset_lyr,cover_poly_fc,"SIMPLIFY","Value")
                                        
        #2.4 turn into layer, select 1, extract to shape
        perim_fc = ''.join([temppath,"\\rough_perim.shp"])
        perim_lyr = "rough_perim"
        if not (os.path.isfile(perim_fc)):
            arcpy.MakeFeatureLayer_management(cover_poly_fc,perim_lyr,""""GRIDCODE" =1""","#","FID FID VISIBLE NONE;Shape Shape VISIBLE NONE;ID ID VISIBLE NONE;GRIDCODE GRIDCODE VISIBLE NONE")
            arcpy.CopyFeatures_management(perim_lyr,perim_fc)

        perim_lyr = perim_fc
        
        #4. Create and save split-rivers with TID into folder in out_folder workspace.
        if not os.path.isdir(transpath):
            arcpy.CreateFolder_management(out_folder,"trans")
        if not os.path.isdir(rivpath):
            arcpy.CreateFolder_management(out_folder,"riv_split")
        if not os.path.isdir(bufpath):
            arcpy.CreateFolder_management(out_folder,"buf")
        if not os.path.isdir(frbpath):
            arcpy.CreateFolder_management(out_folder,"frb")

        riv_files = []
        frb_files = []
        for f in os.listdir(zonespath):
            fullpath= os.path.join(zonespath,f)
            if (os.path.isfile(fullpath) and f.endswith(".shp")):   
                id =  f[5:-4]
                if  (float(id) < -1 ): #change this when you need to skip some zones temporarily!
                    arcpy.AddMessage(''.join(["ZONE ",id, " (skipped)"]))
                    riv_files.append(''.join([rivpath,'\\rivsplit_',id,'.shp']))
                    frb_files.append(''.join([frbpath,"\\frb_",id,".shp"]))
                    continue
                arcpy.AddMessage(''.join(["ZONE ",id]))

                buf_fc = ''.join([bufpath,"\\buf_",id,".shp"])
                buf_clip_fc = ''.join([bufpath,"\\buf_clip_",id,".shp"])
                buf_sp_fc = ''.join([bufpath,"\\buf_sp_",id,".shp"])
                buf_sp_lyr = ''.join(["buf_sp_",id,"_lyr"])
                riv_fc = ''.join([rivpath,'\\rivsplit_',id,'.shp'])
                riv_lyr = ''.join(['rivsplit_',id])
                trans_fc = ''.join([transpath,"\\trans_",id,".shp"])
                frb_fc = ''.join([frbpath,"\\frb_",id,".shp"])
                
                #print "Building transects for zone " , id      
                if (not os.path.isfile(trans_fc)):
                    arcpy.AddMessage(''.join(["Building transects for zone " , id]))
                    LinearTransect_execute(fullpath,transect_length,float(transect_width),"METERS",trans_fc,riv_fc)

                #print "Initial buffer for zone " , id
                if (not os.path.isfile(buf_fc)):
                    arcpy.AddMessage(''.join(["Initial buffer for zone " , id]))
                    arcpy.Buffer_analysis(trans_fc,buf_fc,(str(float(transect_width)/2)+" Meters"),"FULL","FLAT","NONE","#")
                
                #print "Clipping initial buffer to perim for zone " , id
                if (not os.path.isfile(buf_clip_fc)):
                    arcpy.AddMessage(''.join(["Clipping initial buffer to perim for zone " , id]))
                    arcpy.Clip_analysis(buf_fc,perim_lyr,buf_clip_fc,"#")
                    
                if (not os.path.isfile(buf_sp_fc)):
                    arcpy.AddMessage(''.join(["Converting clipped buffer to singleparts for zone " , id]))
                    arcpy.MultipartToSinglepart_management(buf_clip_fc,buf_sp_fc)
                
                #print "Trimming buffer with a tolerated " , tolerance,  "m distance from input hydrology for zone " , id
                
                if (not os.path.isfile(frb_fc)):
                    arcpy.AddMessage(''.join(["Trimming buffer with a tolerated " , tolerance,  "m distance from input hydrology for zone " , id]))
                    arcpy.MakeFeatureLayer_management (buf_sp_fc, buf_sp_lyr)
                    arcpy.SelectLayerByLocation_management(buf_sp_lyr,"WITHIN_A_DISTANCE",in_hydro,tolerance+" Meters","NEW_SELECTION")                 
                    arcpy.CopyFeatures_management(buf_sp_lyr,frb_fc,"#","0","0","0")
                
                #print "Calculating FRBWidth field for zone " ,id
                    arcpy.AddMessage(''.join(["Calculating FRBWidth field for zone ", id]))
                    arcpy.AddField_management(frb_fc,"FRBWidth","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    exp = '!SHAPE.AREA@SQUAREMETERS!/'+transect_width
                    arcpy.CalculateField_management(frb_fc,"FRBWidth",exp,"PYTHON_9.3","#")
                
                    #print "Adding FRBWidth to River-Split for zone  " , id
                    arcpy.AddMessage(''.join(["Adding FRBWidth to River-split for zone " , id]))
                    arcpy.AddField_management(riv_fc, "FRBWidth1", "DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "FRBWidth2", "DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Num_FRB", "DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    
                    arcpy.AddField_management(riv_fc, "Over_7", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_15", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_30", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_45", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_60", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_75", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                    arcpy.AddField_management(riv_fc, "Over_90", "SHORT", "#","#","#","#","NULLABLE","NON_REQUIRED","#")
                
                    #src: http://gis.stackexchange.com/questions/95957/most-efficient-method-to-join-multiple-fields-in-arcgis
                    #Pull values from frb_fc
                    side1 = {}
                    side2 = {}
                    joinfields = ['TransectID','FRBWidth','Side']
                    with arcpy.da.SearchCursor(frb_fc,joinfields) as rows:
                        for row in rows:
                            tid = row[0]
                            if (row[2] == 1):
                                side1[tid] = row[1]
                            else:
                                side2[tid] = row[1]
                    del row, rows
                
                    #Put values into riv_split
                    global globTransectId
                    targfields = ['TransectID', 'FRBWidth1', 'FRBWidth2', 'Num_FRB','Over_7','Over_15','Over_30','Over_45','Over_60','Over_75','Over_90']
                    with arcpy.da.UpdateCursor(riv_fc, targfields) as recs:
                        for rec in recs:
                            over = {7:0,15:0,30:0,45:0,60:0,75:0,90:0}
                            globTransectId += 1
                            tid = rec[0]
                            numFRB = 0
                            if side1.has_key(tid):
                                rec[1] = side1[tid]
                                numFRB += 1
                                for threshVal in over:
                                    if (rec[1] > threshVal):
                                        over[threshVal] += 1
                                
                            if side2.has_key(tid):
                                rec[2] = side2[tid]
                                numFRB += 1
                                for threshVal in over:
                                    if (rec[2] > threshVal):
                                        over[threshVal] += 1
                            rec[3] = numFRB
                            rec[4] = over[7]
                            rec[5] = over[15]
                            rec[6] = over[30]
                            rec[7] = over[45]
                            rec[8] = over[60]
                            rec[9] = over[75]
                            rec[10] = over[90]
                            recs.updateRow(rec)
                        del rec, recs
                    
                    # Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
    # The following inputs are layers or table views: "al_frb_riv"
                    '''arcpy.AddMessage(''.join(["Counting key buffer widths for zone " , id]))

                    CodeBlock_frb="""def nums(f1,f2,w):
 num = 0
 if (f1>w):
  num+=1
 if (f2>w):
  num+=1
 return num"""
                    arcpy.CalculateField_management(riv_fc,"Over_7","nums(!FRBWidth1!,!FRBWidth2!,7)", "PYTHON_9.3", CodeBlock_frb)        
                    arcpy.CalculateField_management(riv_fc,"Over_15","nums(!FRBWidth1!,!FRBWidth2!,15)", "PYTHON_9.3", CodeBlock_frb)
                    arcpy.CalculateField_management(riv_fc,"Over_30","nums(!FRBWidth1!,!FRBWidth2!,30)", "PYTHON_9.3", CodeBlock_frb)
                    arcpy.CalculateField_management(riv_fc,"Over_45","nums(!FRBWidth1!,!FRBWidth2!,45)", "PYTHON_9.3", CodeBlock_frb)
                    arcpy.CalculateField_management(riv_fc,"Over_60","nums(!FRBWidth1!,!FRBWidth2!,60)", "PYTHON_9.3", CodeBlock_frb)
                    arcpy.CalculateField_management(riv_fc,"Over_75","nums(!FRBWidth1!,!FRBWidth2!,75)", "PYTHON_9.3", CodeBlock_frb)
                    arcpy.CalculateField_management(riv_fc,"Over_90","nums(!FRBWidth1!,!FRBWidth2!,90)", "PYTHON_9.3", CodeBlock_frb)'''
                
                riv_files.append(riv_fc)
                frb_files.append(frb_fc)

                #Clean up.
                arcpy.Delete_management(buf_sp_lyr)

                if isDel:
                    arcpy.Delete_management(buf_fc)
                    arcpy.Delete_management(buf_clip_fc)
                    arcpy.Delete_management(buf_sp_fc)
                    arcpy.Delete_management(trans_fc)
                
        #6. Recompile all the created rivers into one file.
        arcpy.Merge_management (riv_files, out_hydro)
        arcpy.Merge_management (frb_files, out_frb)
        #7. Recompile all the created buffer slices into one file(?)
        
        
