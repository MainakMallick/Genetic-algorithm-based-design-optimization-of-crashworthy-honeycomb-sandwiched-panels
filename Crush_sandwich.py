#Python Script for optimization of 2-D honeycomb structures for parameterized input from a text file
#with included Top and Bottom plates subjected to a uniformly distributed Pressure bending load.
#The parameters of drafting such as the total width, length, number of core cells and the wall thickness
#are provided through a text file.
#Developed at the Department of Aerospace Engineering and Applied Mechanics
#Indian Institute of Engineering Science and Technology, Shibpur.
from abaqus import*
from abaqusConstants import*
import regionToolset
import math
#independent quantities taken as user input
iteration_count=1
input_file=open('honeycomb_sandwich.txt')

for line in input_file:
    extracted_line=line
    extracted_word=extracted_line.split()
    #independent quantities taken as user input
    #from the textfile honeycomb_sandwich.txt
    #In order of precedence, these include the total length
    #of the honeycomb, the wall thickness of the core cells, the
    #number of core cells, the top layer thickness and the
    #bottom layer thickness respectively
    total_length_HC=float(extracted_word[0])
    wall_thickness_core_cells=float(extracted_word[1])
    No_of_core_cells=float(extracted_word[2])
    core_layer_thickness=float(extracted_word[3])
    topLayer_thickness=float(extracted_word[4])
    bottomLayer_thickness=float(extracted_word[5])
    # -----------------------------------------------------------------------------------------------------------------#
    # variable initialization: units=metres
    reportxy_name='\SandwichXYData'
    reportxy_path='G:\HC_structure_optimization'
    # the number of core cells are the
    # number of hexagons along each row and each column
    # they have to be the same
    # -----------------------------------------------------------------------------------------------------------------#
    # --------------------------------------INITIAL CALCULATIONS-------------------------------------------------------#

    # --------------------------------No Abaqus Statements appear here-------------------------------------------------#
    no_of_core_cells=math.floor(No_of_core_cells)

    # All the initial calculations are being done in the following section
    # All the scalar and initial quantities such as length, width, angles
    # population of the x-coordinate points and y-coordinate points into arrays
    # definition of points for section assignment, assembly and boundary conditions
    # would be done in the following section. All the ABAQUS Python commands and the
    # ABAQUS statements would follow in a later section.
    # This effort is made to reduce the computational/processing time.

    #-------------------------SPECIAL CALCULATIONS FOR DRAFTING OF CORE LAYER SECTION----------------------------------#

    # The Included angle of a regular hexagon is 120 degrees
    # Excluding the angle of perpendicularity that is 90 degrees
    # the inscribed angle is 30 degrees
    # and ABAQUS processes angles in radians so the angle eventually
    # will be converted into radians
    angle=30.0
    rangle=((math.pi)/180)*angle
    #In the process of parameterization by a script the outer edge will
    #be calculated on the basis of the total length of the honeycomb structure
    #from geometric considerations the distance in between the flats of a regular hexagon
    #is given by: 2*outeredge*cos(theta).Here theta is taken as rangle
    #If the above quantity multiplied by number of core cells, the total length of the honeycomb obtained.
    #outer edge is core_layer_edge_outer
    core_layer_edge_outer=(total_length_HC/(2*no_of_core_cells*math.cos(rangle)))
    #While Drafting the core layer section we are going to require certain quantities which
    #we will compute or calculate here itself
    #Our idea is to populate two empty arrays which would contain the values of the corresponding x-coordinates
    #and the y-coordinates. Firstly, our idea would be to create a blank structure as shown below for
    # 3 number of core cells
                  #    #    #
                #   #    #    #
                #   #    #    #
                  #    #    #
    # Calculations are as below
    # the expression to be used in a while loop for increment of x-coordinate for the above blank structure
    x_incr=(total_length_HC/(2*no_of_core_cells))
    # On reaching the limit of the x-increment with total length of honeycomb y must
    # be incremented depending on the the number of core cells. The increment would vary
    # with odd or even distribution of core cells
    if no_of_core_cells%2==0:
        y_incr=int((no_of_core_cells/2)+1)*(2*core_layer_edge_outer*math.sin(rangle)+core_layer_edge_outer)+\
               int(no_of_core_cells/2)*core_layer_edge_outer-\
               core_layer_edge_outer-2*core_layer_edge_outer*math.sin(rangle)
    else:
        y_incr=int(no_of_core_cells/2+1)*(2*core_layer_edge_outer*math.sin(rangle)+core_layer_edge_outer)+\
               int(no_of_core_cells/2)*core_layer_edge_outer-2*core_layer_edge_outer*math.sin(rangle)
    # For initially drawing the blank area of the core layer we would be using while loops for populating
    # the arrays of the x-coordinates and the y-coordinates. We would use to check initially whether or not
    # the current coordinate is less than that of the total length of the honeycomb. So if the total length
    # for instance is 50m and the last incremented value is 49.996m (which is approx 50m), so it
    # would be considered less than 50m and the as per the next increment the x-value would exceed 50m.
    # This is undesirable and hence we introduce a variable called clearcheck to negate it from the
    # total length of the honeycomb such that the above error is avoided. This would also be helpful
    # when we proceed from total_length_HC to zero. The clearcheck value is set as the wall thickness
    # of the core cells.
    clearcheck=wall_thickness_core_cells
    # Now we take two empty arrays for x-coordinate and y-coordinate and populate them
    # exactly as per our requirement of the blank
    vertices_x=[]
    vertices_y=[]
    # The x-coordinate would be at zero while the y-coordinate will be at a certain height
    xperiod=0.0
    yperiod=core_layer_edge_outer*math.sin(rangle)
    # A variable called count would be introduced to count the number of points of the outer blank
    # to be drafted. And another variable called const would be introduced. This is for arranging the
    # y-coordinates in a zig-zag way such that the outer boundaries of the blank exactly
    # match with that of the structure of a honeycomb.
    count=0
    const=1
    #The constant will be multiplied by -1 continnuously to vary y-coordinate in a zigzag fashion
    while xperiod<total_length_HC-clearcheck:
        vertices_x.append(xperiod)
        vertices_y.append(yperiod)
        count=count+1
        xperiod=xperiod+x_incr
        const=const*(-1)
        yperiod=yperiod+(const)*core_layer_edge_outer*math.sin(rangle)
    vertices_x.append(xperiod)
    vertices_y.append(yperiod)
    count=count+1
    #After the program encounters the end of the length of honeycomb
    #the y-coordinate is incremented by the y_incr value
    yperiod=yperiod+y_incr
    #Now incase of decrease from the total length of the honeycomb to zero
    #A variable called var is introduced which is similar to const
    #It varies with the odd number of cells or the even number of cells
    if no_of_core_cells%2==0:
        var=-1
    else:
        var=1
    while xperiod>clearcheck:
        vertices_x.append(xperiod)
        vertices_y.append(yperiod)
        count=count+1
        xperiod=xperiod-x_incr
        yperiod=yperiod+var*core_layer_edge_outer*math.sin(rangle)
        var=var*(-1)
    vertices_x.append(xperiod)
    vertices_y.append(yperiod)
    count=count+1
    #We will now compute the inner edge of the regular hexagon
    #This will be done by comparing with the length formed by the outer edges
    #The total wall thickness subtracted from the outer edge was found to be the
    #inner edge of the hexagon. This would be needed for doing Cut-Extrude of
    #hexagons inside the blank structure.
    edge_inner=core_layer_edge_outer-wall_thickness_core_cells
    #An important quantity is the total width of the honeycomb
    #This would be a dependent quantity and would be dependent on the
    #number of core cells
    if no_of_core_cells%2==0:
        total_width_HC=y_incr+core_layer_edge_outer*math.sin(rangle)
    else:
        total_width_HC=y_incr+2*core_layer_edge_outer*math.sin(rangle)
    #Another variable is the ydistmeasure which would be used for the
    #checking of the increment of the y-coordinate
    ydistmeasure=total_width_HC-core_layer_edge_outer*math.sin(rangle)
    #----------------------calculation for the top face point where the partition would be carried out-----------------#

    sketch_top_face_point=(total_length_HC/2,total_width_HC/2,core_layer_thickness)

    #--------Left edge is needed to be determined for making the sketch transform for partitioning---------------------#

    sketch_left_edge_point=(0.0,core_layer_edge_outer*math.sin(rangle) + \
                              no_of_core_cells*core_layer_edge_outer,core_layer_thickness)

    #--------------COORDINATE POINT ASSIGNMENT FOR ASSEMBLY OF THE HONEYCOMB SANDWICHED STRUCTURE----------------------#

    #-----------------------(i) the front face and the rear face points of the top layer ------------------------------#

    topLayer_front_face_point=(0.0,total_width_HC/2,topLayer_thickness/2)
    topLayer_back_face_point=(total_length_HC,total_width_HC/2,topLayer_thickness/2)

    #-----------------------------------(ii) the bottom face of the top layer------------------------------------------#

    topLayer_bottom_face_point=(total_length_HC/2,total_width_HC/2,0.0)

    #--------------------------------(ii(a)) top layer and bottom layer edge point ------------------------------------#

    topLayer_bottom_edge_point=(total_length_HC/2,total_width_HC,0.0)
    bottomLayer_top_edge_point=(total_length_HC/2,total_width_HC,bottomLayer_thickness)

    #------------------------(iii) the front face and the rear face points of the bottom layer-------------------------#

    bottomLayer_front_face_point=(0.0,total_width_HC/2,topLayer_thickness/2)
    bottomLayer_back_face_point=(total_length_HC,total_width_HC/2,topLayer_thickness/2)

    #----------------------------------(iv) the top face of the bottom layer ------------------------------------------#
    bottomLayer_top_face_point=(total_length_HC/2,total_width_HC/2,bottomLayer_thickness)

    #---------------------Important points for coincident assembly-----------------------------------------------------#

    #This is again dependent on odd or even number of core cells
    if no_of_core_cells%2==0:
        topLayer_bottom_partition_point=(total_length_HC-2*core_layer_edge_outer*math.cos(rangle),total_width_HC,0.0)
        bottomLayer_top_partition_point=(total_length_HC-2*core_layer_edge_outer*math.cos(rangle),total_width_HC,\
                                         bottomLayer_thickness)
        coreLayer_top_align_point=(total_length_HC-2*core_layer_edge_outer*math.cos(rangle),total_width_HC,\
                                   core_layer_thickness)
        coreLayer_bottom_align_point=(total_length_HC-2*core_layer_edge_outer*math.cos(rangle),total_width_HC,0.0)
    else:
        topLayer_bottom_partition_point=(total_length_HC-core_layer_edge_outer*math.cos(rangle),total_width_HC,0.0)
        bottomLayer_top_partition_point=(total_length_HC-core_layer_edge_outer* math.cos(rangle),total_width_HC, \
                                         bottomLayer_thickness)
        coreLayer_top_align_point=(total_length_HC-core_layer_edge_outer*math.cos(rangle),total_width_HC, \
                                   core_layer_thickness)
        coreLayer_bottom_align_point=(total_length_HC-core_layer_edge_outer*math.cos(rangle),total_width_HC,0.0)

    topLayer_bottom_partition_point_1 = (total_length_HC /2, total_width_HC, 0.0)
    #-------------------------------evaluate a clearance distance by the name of clearDist-----------------------------#

    clearDist=wall_thickness_core_cells*math.cos(rangle)

    # Point used to find the top surface of the core
    # We will be using the findAt command to find the top surface of the core
    # As an arguement we need to find the coordinates of a point on this surface
    # For an even number of cells, the exact center point of the top surface of
    # the core can be used
    # However for an odd number of cells, this point will lie over one of the holes
    # In that case we will pick a point offset a little way from it

    if no_of_core_cells%2==0:
        coreLayer_top_face_point=(total_length_HC-(edge_inner*math.cos(rangle)+clearDist),\
                                  total_width_HC-core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,\
                                  core_layer_thickness)
        coreLayer_bottom_face_point=(total_length_HC-(edge_inner*math.cos(rangle)+clearDist),\
                                     total_width_HC-core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,0.0)
        coreLayer_inside_coord=(total_length_HC-(edge_inner*math.cos(rangle)+clearDist),\
                                total_width_HC-core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,\
                                core_layer_thickness/2)
    else:
        coreLayer_top_face_point=(total_length_HC-(core_layer_edge_outer*math.cos(rangle)+\
                                  edge_inner*math.cos(rangle)+clearDist),total_width_HC-\
                                  core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,\
                                  core_layer_thickness)
        coreLayer_bottom_face_point=(total_length_HC-(core_layer_edge_outer*math.cos(rangle)+\
                                     edge_inner*math.cos(rangle)+clearDist),total_width_HC-\
                                     core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,0.0)
        coreLayer_inside_coord=(total_length_HC-(core_layer_edge_outer*math.cos(rangle)+\
                                edge_inner*math.cos(rangle)+clearDist),total_width_HC-\
                                core_layer_edge_outer*math.sin(rangle)-core_layer_edge_outer/2,core_layer_thickness/2)

    #--------------------------Inside Coords for Top Layer and Bottom Layer/for meshing--------------------------------#
    topLayer_inside_coord=(total_length_HC/2,total_length_HC/2,topLayer_thickness/2)
    bottomLayer_inside_coord=(total_length_HC/2,total_length_HC/2,bottomLayer_thickness/2)
    #-------------------------------Assign the top layer top surface point---------------------------------------------#

    topLayer_top_surface_point=(total_length_HC/2,total_width_HC/2,topLayer_thickness)

    #-----------------------------Assign the top layer bottom Surface Point--------------------------------------------#

    topLayer_bottom_surface_point=(total_length_HC/2,total_width_HC/2,0.0)

    #-----------------------------Assign the bottom layer top surface point--------------------------------------------#

    bottomLayer_top_surface_point=(total_length_HC/2,total_width_HC/2,bottomLayer_thickness)

    #---------------------------Assign the crusher bottom surface point------------------------------------------------#

    crusher_bottom_surface_point=(total_length_HC/2,total_width_HC/2,0.0)

    #------------COORDINATE POINT ASSIGNMENT FOR BOUNDARY CONDITIONS OF THE HONEYCOMB SANDWICHED STRUCTURE-------------#

    # ---------------------------Points for fixing the front and back face of the top layer----------------------------#

    topLayer_fix_front_face_point=(0.0,total_width_HC/2,topLayer_thickness/2)
    topLayer_fix_back_face_point=(total_length_HC,total_width_HC/2,topLayer_thickness/2)

    #--------------------------Points for fixing the front and back face of the bottom layer---------------------------#

    bottomLayer_fix_front_face_point=(0.0,total_width_HC/2,bottomLayer_thickness/2)
    bottomLayer_fix_back_face_point=(total_length_HC,total_width_HC/2,bottomLayer_thickness/2)

    vertex_coords_for_displacement_1 = (0.0,0.0,0.0)

    #------------------------Points for fixing the front and back face of the core layer-------------------------------#

    coreLayer_fix_front_face_points=[]
    coreLayer_fix_back_face_points=[]

    tworangle=2*rangle

    coord_x_front=0.0
    coord_x_back=total_length_HC
    wall_extra=2*wall_thickness_core_cells*math.cos(rangle)
    first_ydist=wall_thickness_core_cells*math.sin(rangle)+edge_inner+wall_extra*math.sin(tworangle)
    second_ydist=2*wall_extra*math.sin(tworangle)+edge_inner
    coord_y=core_layer_edge_outer*math.sin(rangle)+(first_ydist/2)
    coreLayer_fix_front_face_points.append((coord_x_front,coord_y,core_layer_thickness/2))
    coreLayer_fix_back_face_points.append((coord_x_back,coord_y,core_layer_thickness/2))

    coord_y=coord_y+(first_ydist/2)+2*edge_inner*math.sin(rangle)+edge_inner+(second_ydist/2)

    face_counter = 1
    while coord_y<total_width_HC:
        temp_point_front=(coord_x_front,coord_y,core_layer_thickness/2)
        temp_point_back=(coord_x_back,coord_y,core_layer_thickness/2)
        coreLayer_fix_front_face_points.append(temp_point_front)
        coreLayer_fix_back_face_points.append(temp_point_back)
        coord_y=coord_y+(second_ydist/2+2*edge_inner*math.sin(rangle)+edge_inner+second_ydist/2)
        face_counter=face_counter+1
   # -----------------------Determine the load top face point of the top layer-----------------------------------------#
    topLayer_load_top_face_point = (total_length_HC / 2, total_width_HC / 2, topLayer_thickness)
    #-----------------------Determine the bottom layer face point of the bottom layer----------------------------------#

    bottomLayer_bottom_face_point=(total_length_HC/2,total_width_HC/2,0.0)

    #-------------------------------COORDINATE POINTS FOR POST PROCESSING----------------------------------------------#

    if no_of_core_cells%2==0:
        xdist=total_length_HC/2-core_layer_edge_outer*math.cos(rangle)
    else:
        xdist=total_length_HC/2


    vertex_coords_for_displacement_2=(xdist,0.0,0.0)


    #--------------------------------------ABAQUS STATEMENTS START HERE------------------------------------------------#

    session.viewports['Viewport: 1'].setValues(displayedObject=None)

    #-------------------------------------------Create the Model-------------------------------------------------------#

    # Change the Model Name to HoneyComb Structure
    # Add the iteration count to the end of the Model's name
    modelname='HoneyComb Structure'+repr(iteration_count)
    #change the name of the model from Model-1 to HoneyComb Part
    #mdb.models.changeKey(fromName='Model-1',toName='HoneyComb Part')
    #set the object name hcmodel to HoneyComb Part
    hcModel=mdb.Model(name=modelname)

    #----------------------------------------------PART/SKETCH MODULE--------------------------------------------------#
    import part
    import sketch
    #--------------------------------------------Draft the Top Layer---------------------------------------------------#
    #Create a 3D deformable part with the name of "Top Layer"
    topLayerPart=hcModel.Part(name='Top Layer',dimensionality=THREE_D,type=DEFORMABLE_BODY)
    #Now create a sketch
    topLayerProfileSketch=hcModel.ConstrainedSketch(name='Top Layer Sketch', sheetSize=40.0)
    #Draw the geometry for the sketch
    topLayerProfileSketch.rectangle(point1=(0.0,0.0),point2=(total_length_HC,total_width_HC))
    #Extrude the above sketch
    topLayerPart.BaseSolidExtrude(sketch=topLayerProfileSketch,depth=topLayer_thickness)
    #Find the edge on which to create a partition
    topLayer_bottom_edge=topLayerPart.edges.findAt(topLayer_bottom_edge_point,)
    #Create partition of the edge by a point
    topLayerPart.PartitionEdgeByPoint(edge=topLayer_bottom_edge, point=topLayer_bottom_partition_point)

    reference_point = topLayerPart.ReferencePoint(point=(total_length_HC / 2, total_width_HC / 2, topLayer_thickness))
    topLayerPart.features.changeKey(fromName='RP', toName='Top Layer')


    #-------------------------------------------Draft the Core Layer---------------------------------------------------#

    #Create a 3D deformable part with the name of "HoneyComb Part"
    #For Drawing Hexagonal honeycombs we will draw a blank pointed structure first
    coreLayerPart=hcModel.Part(name='HoneyComb Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    #Now create a sketch for the blank structure
    hcSketch=hcModel.ConstrainedSketch(name='Sketch for HoneyComb part', sheetSize=100.0)
    #As the vertices are populated as per the lines .........
    #the sketch of the blank structure is drawn as follows
    for k in range(0, count-1):
        hcSketch.Line(point1=(vertices_x[k],vertices_y[k]), point2=(vertices_x[k+1], vertices_y[k+1]))
    hcSketch.Line(point1=(vertices_x[count-1], vertices_y[count-1]), point2=(vertices_x[0],vertices_y[0]))
    #Firstly use the BaseSolidExtrude() function and extrude the hc part
    coreLayerPart.BaseSolidExtrude(sketch=hcSketch, depth=core_layer_thickness)
    # After the blank structure is drafted we will create a sketch to extrude the draft
    # Just Cut out the through hole in the form of the honeyComb pattern
    # select the top surface of the created as per the point defined in line 157
    sketch_top_face=coreLayerPart.faces.findAt(sketch_top_face_point, )
    # select the left edge point in line 159
    sketch_left_edge=coreLayerPart.edges.findAt(sketch_left_edge_point, )
    # now for partitioning a drawing canvas would be brought up by sketch transform
    sketch_transform=coreLayerPart.MakeSketchTransform(sketchPlane=sketch_top_face, \
                                                       sketchUpEdge=sketch_left_edge, \
                                                       sketchPlaneSide=SIDE1, \
                                                       sketchOrientation=LEFT, \
                                                       origin=(0.0, 0.0, 0.0))
    # create an object insideSketch for drawing the inside honeycomb pattern
    insideSketch=hcModel.ConstrainedSketch(name='core layer cutout sketch', sheetSize=50.0, transform=sketch_transform)
    # Drafting the inside sketch and making the Cut Extrude
    nrow=1
    coor_y_start=vertices_y[0]
    while coor_y_start<ydistmeasure:
        if nrow%2==0:
            coor_x_start=vertices_x[0]-(core_layer_edge_outer)*math.cos(rangle)
        else:
            coor_x_start=vertices_x[0]
        coor_ini_x=coor_x_start+wall_thickness_core_cells*math.cos(rangle)
        coor_ini_y=coor_y_start+wall_thickness_core_cells*math.sin(rangle)
        while coor_ini_x<total_length_HC:
            coor_x_i=[]
            coor_y_i=[]
            # Populate the coordinates of inside x and y
            # Populate the coordinates of the x array
            coor_x_i.append(coor_ini_x)
            coor_x_i.append(coor_ini_x)
            coor_x_i.append(coor_ini_x+edge_inner*math.cos(rangle))
            coor_x_i.append(coor_ini_x+2*edge_inner*math.cos(rangle))
            coor_x_i.append(coor_ini_x+2*edge_inner*math.cos(rangle))
            coor_x_i.append(coor_ini_x+edge_inner*math.cos(rangle))
            # Populate the coordinates of the y array
            coor_y_i.append(coor_ini_y)
            coor_y_i.append(coor_ini_y+edge_inner)
            coor_y_i.append(coor_ini_y+edge_inner+edge_inner*math.sin(rangle))
            coor_y_i.append(coor_ini_y+edge_inner)
            coor_y_i.append(coor_ini_y)
            coor_y_i.append(coor_ini_y-edge_inner*math.sin(rangle))
            # Simply connect the vertices
            insideSketch.Line(point1=(coor_x_i[0], coor_y_i[0]), \
                              point2=(coor_x_i[1], coor_y_i[1]))
            insideSketch.Line(point1=(coor_x_i[1], coor_y_i[1]), \
                              point2=(coor_x_i[2], coor_y_i[2]))
            insideSketch.Line(point1=(coor_x_i[2], coor_y_i[2]), \
                              point2=(coor_x_i[3], coor_y_i[3]))
            insideSketch.Line(point1=(coor_x_i[3], coor_y_i[3]), \
                              point2=(coor_x_i[4], coor_y_i[4]))
            insideSketch.Line(point1=(coor_x_i[4], coor_y_i[4]), \
                              point2=(coor_x_i[5], coor_y_i[5]))
            insideSketch.Line(point1=(coor_x_i[5], coor_y_i[5]), \
                              point2=(coor_x_i[0], coor_y_i[0]))

            coor_ini_x=coor_ini_x+2*edge_inner*math.cos(rangle)+2*wall_thickness_core_cells*math.cos(rangle)

        coor_y_start=coor_y_start+(core_layer_edge_outer*math.sin(rangle))+core_layer_edge_outer
        nrow=nrow+1

    # Cut-Extrude Through Holes in a HoneyComb Pattern
    coreLayerPart.CutExtrude(sketchPlane=sketch_top_face, \
                             sketchUpEdge=sketch_left_edge, \
                             sketchPlaneSide=SIDE1, sketchOrientation=LEFT, \
                             sketch=insideSketch, \
                             flipExtrudeDirection=OFF)

    # ----------------------------------------Draft the Bottom Layer---------------------------------------------------#

    # Create a 3D deformable part with the name of "Bottom Layer"
    bottomLayerPart=hcModel.Part(name='Bottom Layer',dimensionality=THREE_D,type=DEFORMABLE_BODY)
    # Now create a sketch
    bottomLayerProfileSketch=hcModel.ConstrainedSketch(name='Bottom Layer Sketch',sheetSize=40.0)
    # Draw the geometry for the sketch
    bottomLayerProfileSketch.rectangle(point1=(0.0,0.0),point2=(total_length_HC,total_width_HC))
    # Extrude the above sketch
    bottomLayerPart.BaseSolidExtrude(sketch=bottomLayerProfileSketch, depth=bottomLayer_thickness)
    # Find the edge on which to create a partition
    bottomLayer_top_edge=bottomLayerPart.edges.findAt(bottomLayer_top_edge_point, )
    # Create partition of the edge by a point
    bottomLayerPart.PartitionEdgeByPoint(edge=bottomLayer_top_edge, point=bottomLayer_top_partition_point)
    # Create bottom Layer Bottom Edge
    bottomLayer_bottom_edge=bottomLayerPart.edges.findAt(vertex_coords_for_displacement_2,)
    # Create an edge partition for bottom layer
    bottomLayerPart.PartitionEdgeByPoint(edge=bottomLayer_bottom_edge, point=vertex_coords_for_displacement_2)

    crashLayerPart=hcModel.Part(name='Crashing Layer',dimensionality=THREE_D,type=DISCRETE_RIGID_SURFACE)
    #Now create a sketch
    crashLayerProfileSketch=hcModel.ConstrainedSketch(name='crash Layer Sketch', sheetSize=40.0)
    #Draw the geometry for the sketch
    crashLayerProfileSketch.rectangle(point1=(0.0,0.0),point2=(total_length_HC,total_width_HC))
    #Extrude the above sketch
    crashLayerPart.BaseShell(sketch=crashLayerProfileSketch)
    #Find the edge on which to create a partition
    reference_point=crashLayerPart.ReferencePoint(point=(total_length_HC/2,total_width_HC/2, 0.0))
    crashLayerPart.features.changeKey(fromName='RP', toName='Crash Layer')

    r = crashLayerPart.referencePoints
    refPoints = (r[2],)
    region = crashLayerPart.Set(referencePoints=refPoints, name='Set-1')
    crashLayerPart.engineeringFeatures.PointMassInertia(
        name='Inertia-1', region=region, mass=2000, alpha=0.0, composite=0.0)

    #------------------------------------------MATERIAL MODULE---------------------------------------------------------#

    import material

    # Create a material AlSi10Mg by assigning mass density, Young's Modulus and Poisson's ratio
    hcMaterial=hcModel.Material(name='AlSi10Mg')
    hcMaterial.Density(table=((2800,),))
    hcMaterial.Elastic(table=((70E9, 0.33),))
    hcMaterial.Plastic(hardening=JOHNSON_COOK,
                       table=((520, 477, .52, 1, 893, 293),))
    hcMaterial.plastic.RateDependent(
        type=JOHNSON_COOK, table=((.001, .0005),))

    #------------------------------------------SECTION MODULE----------------------------------------------------------#

    import section
    # Create the solid section to assign the crusher to it
    # Create the solid section and assign the beam to it
    # Create a section to assign to the top layer
    topLayerSection=hcModel.HomogeneousSolidSection(name='Top Layer Section', \
                                                    material='AlSi10Mg')
    # Assign the Top Layer to this section
    top_layer_region=(topLayerPart.cells,)
    topLayerPart.SectionAssignment(region=top_layer_region, \
                                   sectionName='Top Layer Section')
    # Create a section to assign to the bottom layer
    bottomLayerSection=hcModel.HomogeneousSolidSection(name='Bottom Layer Section', \
                                                       material='AlSi10Mg')
    # Assign the Bottom Layer to this section
    bottom_layer_region=(bottomLayerPart.cells,)
    bottomLayerPart.SectionAssignment(region=bottom_layer_region, \
                                      sectionName='Bottom Layer Section')
    # Create a section to assign to the core layer
    coreLayerSection=hcModel.HomogeneousSolidSection(name='Core Layer Section', \
                                                     material='AlSi10Mg')
    # Assign the Core Layer to this section
    core_layer_region=(coreLayerPart.cells,)
    coreLayerPart.SectionAssignment(region=core_layer_region, \
                                    sectionName='Core Layer Section')

    # -----------------------------------------ASSEMBLY MODULE---------------------------------------------------------#

    import assembly

    # This import assembly section is for creating instances of a part
    # Create the part instances
    hcAssembly=hcModel.rootAssembly
    topLayerInstance=hcAssembly.Instance(name='Top Layer Instance',part=topLayerPart,dependent=ON)
    bottomLayerInstance=hcAssembly.Instance(name='Bottom Layer Instance',part=bottomLayerPart,dependent=ON)
    coreLayerInstance=hcAssembly.Instance(name='Core Layer Instance',part=coreLayerPart,dependent=ON)
    crashLayerInstance = hcAssembly.Instance(name='Crash Layer Instance', part=crashLayerPart, dependent=ON)
    hcAssembly.translate(instanceList=('Crash Layer Instance',), vector=(0.0, 0.0, core_layer_thickness+topLayer_thickness))
    vertex_for_displacement_1 = topLayerInstance.vertices.findAt((vertex_coords_for_displacement_1,))
    hcAssembly.Set(vertices=vertex_for_displacement_1, name='displacement point set 1')
    prop=hcAssembly.getMassProperties(specifyDensity=True,density=2800.0,miAboutCenterOfMass=True)['mass']
    file_output_2=open(reportxy_path + '\honeycombstructure_output_3.txt', 'w')
    file_output_2.write(repr(prop) + "\n")
    file_output_2.close()
    # -----------------------------------Define the surfaces of the assembly first-------------------------------------#

    topLayer_bottom_surface = topLayerInstance.faces.findAt((topLayer_bottom_surface_point,))

    bottomLayer_top_surface = bottomLayerInstance.faces.findAt((bottomLayer_top_surface_point,))

    coreLayer_top_surface_point = coreLayer_top_face_point
    coreLayer_top_surface = coreLayerInstance.faces.findAt((coreLayer_top_surface_point,))

    coreLayer_bottom_surface_point = coreLayer_bottom_face_point
    coreLayer_bottom_surface = coreLayerInstance.faces.findAt((coreLayer_bottom_surface_point,))

    # -----------------------------Define the face where pressure load would be applied---------------------------------#

    topLayer_load_top_face = topLayerInstance.faces.findAt((topLayer_load_top_face_point,))

    #-----------------------------Define the face which would be fixed for BC------------------------------------------#

    bottomLayer_fix_face=bottomLayerInstance.faces.findAt((bottomLayer_bottom_face_point,))

    #----------------------------Define the faces where the encastre BCs would be applied------------------------------#

    #---------------------------------------------For Top Layer--------------------------------------------------------#

    topLayer_fix_front_face=topLayerInstance.faces.findAt((topLayer_fix_front_face_point,))
    topLayer_fix_front_region=regionToolset.Region(faces=(topLayer_fix_front_face))

    topLayer_fix_back_face=topLayerInstance.faces.findAt((topLayer_fix_back_face_point,))
    topLayer_fix_back_region=regionToolset.Region(faces=(topLayer_fix_back_face))

    #---------------------------------------------For Core Layer-------------------------------------------------------#

    #We will assign at the time of designating BCs

    #---------------------------------------------For Bottom Layer-----------------------------------------------------#

    bottomLayer_fix_front_face=bottomLayerInstance.faces.findAt((bottomLayer_fix_front_face_point,))
    bottomLayer_fix_front_region=regionToolset.Region(faces=(bottomLayer_fix_front_face))

    bottomLayer_fix_back_face = bottomLayerInstance.faces.findAt((bottomLayer_fix_back_face_point,))
    bottomLayer_fix_back_region = regionToolset.Region(faces=(bottomLayer_fix_back_face))

    #----------------------------------- and then the actual assembly operations --------------------------------------#

    core_top_point=coreLayerInstance.vertices.findAt(coreLayer_top_align_point, )
    core_down_point=coreLayerInstance.vertices.findAt(coreLayer_bottom_align_point, )
    top_down_point=topLayerInstance.vertices.findAt(topLayer_bottom_partition_point, )
    bottom_up_point=bottomLayerInstance.vertices.findAt(bottomLayer_top_partition_point, )

    # carry out assembly operation
    # Identify the point on the top layer to be placed with core layer

    hcAssembly.CoincidentPoint(movablePoint=bottom_up_point,fixedPoint=core_down_point)
    hcAssembly.CoincidentPoint(movablePoint=top_down_point,fixedPoint=core_top_point)
    #------------------------IMPORT STEP AND SET FIELD AND HISTORY OUTPUT REQUESTS-------------------------------------#

    # --- Apply Predefined Field-----------

    import step

    # Create the loading step
    hcModel.ExplicitDynamicsStep(name='Apply Load', previous='Initial', \
                                 description='Apply pressure load on top surface of the honeycomb structure',
                                               timePeriod=0.0003)

    hcModel.ContactProperty('IntProp-1')
    hcModel.interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
        table=((0.2,),), shearStressLimit=None, maximumElasticSlip=FRACTION,
        fraction=0.005, elasticSlipStiffness=None)
    hcModel.interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT)
    hcModel.ContactExp(name='Int-1', createStepName='Apply Load')
    hcModel.interactions['Int-1'].includedPairs.setValuesInStep(
        stepName='Apply Load', useAllstar=ON)
    hcModel.interactions['Int-1'].contactPropertyAssignments.appendInStep(
        stepName='Apply Load', assignments=((GLOBAL, SELF, 'IntProp-1'),))



    displacement_point_region_1 = hcAssembly.sets['displacement point set 1']

    hcModel.HistoryOutputRequest(name='Displacement output 1', \
                                 createStepName='Apply Load', \
                                 variables=('UT',), \
                                 region=displacement_point_region_1, \
                                 sectionPoints=DEFAULT, rebar=EXCLUDE)

    #------------------------------------------IMPORT INTERACTION -----------------------------------------------------#

    import interaction

    # import interaction properties

    hcAssembly.Surface(side1Faces=topLayer_bottom_surface, name='Top Layer Bottom')
    hcAssembly.Surface(side1Faces=bottomLayer_top_surface, name='Bottom Layer Top')
    hcAssembly.Surface(side1Faces=coreLayer_bottom_surface, name='Core Layer Bottom')
    hcAssembly.Surface(side1Faces=coreLayer_top_surface, name='Core Layer Top')

    regionDef=hcModel.rootAssembly.surfaces['Core Layer Bottom']
    hcModel.HistoryOutputRequest(name='H-Output-3',
        createStepName='Apply Load', variables=('CFTM', 'CFT1', 'CFT2',
        'CFT3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    region1 = hcAssembly.surfaces['Core Layer Top']
    region2 = hcAssembly.surfaces['Top Layer Bottom']

    hcModel.Tie(name='Constraint-1', master=region1, slave=region2, \
                adjust=ON, tieRotations=ON, \
                constraintEnforcement=SURFACE_TO_SURFACE)

    region1 = hcAssembly.surfaces['Core Layer Bottom']
    region2 = hcAssembly.surfaces['Bottom Layer Top']
    hcModel.Tie(name='Constraint-2', master=region1, slave=region2, \
                adjust=ON, tieRotations=ON, \
                constraintEnforcement=SURFACE_TO_SURFACE)


    # ---------------------------------IMPORT LOAD AND APPLY BOUNDARY CONDITIONS----------------------------------------#

    # --- Apply Predefined Field-----------

    # ---------------------------------IMPORT LOAD AND APPLY BOUNDARY CONDITIONS----------------------------------------#

    # ---TOP LAYER----


    # ---CORE LAYER----


    # --BOTTOM LAYER---

    hcModel.EncastreBC(name='Fix Bottom Layer Front', \
                       createStepName='Apply Load', \
                       region=bottomLayer_fix_front_region)

    hcModel.EncastreBC(name='Fix Bottom Layer Back', \
                       createStepName='Apply Load', \
                       region=bottomLayer_fix_back_region)

    # Apply Load

    bottomLayer_fix_face_region = regionToolset.Region(faces=(bottomLayer_fix_face))

    # ----fix bottom layer
    hcModel.EncastreBC(name='Fix Bottom Layer', \
                       createStepName='Apply Load', \
                       region=bottomLayer_fix_face_region)
    r1 = hcAssembly.instances['Crash Layer Instance'].referencePoints
    refPoints1=(r1[2], )
    region =hcAssembly.Set(referencePoints=refPoints1, name='Set-17')
    hcModel.Velocity(name='Predefined Field-1',
        region=region, field='', distributionType=MAGNITUDE, velocity1=0.0,
        velocity2=0.0, velocity3=-10.0, omega=0.0)
    r3 = hcAssembly.instances['Crash Layer Instance'].referencePoints
    refPoints2=(r3[2], )
    region = hcAssembly.Set(referencePoints=refPoints2, name='Set-18')
    hcModel.DisplacementBC(name='BC-17',
        createStepName='Apply Load', region=region, u1=SET, u2=SET, u3=UNSET,
        ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
        fieldName='', localCsys=None)
    # Apply Load

    topLayer_top_face_region = regionToolset.Region(side1Faces=topLayer_load_top_face)

    # ----loading sequence----

    # apply pressure load on this region


    # ----------------------------------------IMPORT MESH---------------------------------------------------------------#

    import mesh

    # ----------------------------------------------------------------------------------------------------
    # Mesh the Top Layer
    # We place a point somewhere inside it based on our knowledge of the geometry
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, \
                              kinematicSplit=AVERAGE_STRAIN, \
                              secondOrderAccuracy=OFF, \
                              hourglassControl=DEFAULT, distortionControl=DEFAULT)
    topLayerCells = topLayerPart.cells
    selectedTopLayerCells = topLayerCells.findAt(topLayer_inside_coord, )
    topLayerMeshRegion = (selectedTopLayerCells,)
    topLayerPart.setElementType(regions=topLayerMeshRegion, elemTypes=(elemType1,))

    topLayerPart.seedPart(size=0.04, deviationFactor=0.1)

    topLayerPart.generateMesh()
    # -------------------------------------------------------------------------------------------------------
    # Mesh the bottom Layer
    # We place a point somewhere inside it based on our knowledge of the geometry

    elemType2 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, \
                              kinematicSplit=AVERAGE_STRAIN, \
                              secondOrderAccuracy=OFF, \
                              hourglassControl=DEFAULT, distortionControl=DEFAULT)
    bottomLayerCells = bottomLayerPart.cells
    selectedBottomLayerCells = bottomLayerCells.findAt(bottomLayer_inside_coord, )
    bottomLayerMeshRegion = (selectedBottomLayerCells,)
    bottomLayerPart.setElementType(regions=bottomLayerMeshRegion, elemTypes=(elemType2,))

    bottomLayerPart.seedPart(size=0.04, deviationFactor=0.1)

    bottomLayerPart.generateMesh()
    # --------------------------------------------------------------------------------------------------------
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, \
                              kinematicSplit=AVERAGE_STRAIN, \
                              secondOrderAccuracy=OFF, \
                              hourglassControl=DEFAULT, distortionControl=DEFAULT)
    coreLayerCells = coreLayerPart.cells
    selectedCoreLayerCells = coreLayerCells.findAt(coreLayer_inside_coord, )
    coreLayerMeshRegion = (selectedCoreLayerCells,)
    coreLayerPart.setMeshControls(elemShape=TET, regions=coreLayerMeshRegion, technique=FREE, allowMapped=FALSE)
    coreLayerPart.setElementType(regions=coreLayerMeshRegion, elemTypes=(elemType3,))

    coreLayerPart.seedPart(size=0.01, deviationFactor=0.1)

    coreLayerPart.generateMesh()

    p = mdb.models['HoneyComb Structure1'].parts['Crashing Layer']
    p.seedPart(size=0.04, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['HoneyComb Structure1'].parts['Crashing Layer']
    p.generateMesh()
    # -----------------------------------------SUBMIT JOB SECTION------------------------------------------------------#

    import job

    # Create the job
    job_name = 'SandwichStructureJob' + repr(iteration_count)

    mdb.Job(name=job_name, model=modelname, type=ANALYSIS, \
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, \
            description='Run the contact simulation', \
            parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, \
            numDomains=1, userSubroutine='', numCpus=1, memory=50, \
            memoryUnits=PERCENTAGE, scratch='', echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)
    # Run the job
    mdb.jobs[job_name].submit(consistencyChecking=OFF)

    # Do not return control till job is finished running
    mdb.jobs[job_name].waitForCompletion()

    # -------------------------------------------POST PROCESSING SECTION------------------------------------------------#

    import odbAccess
    import visualization

    honeycomb_odb_path = job_name + '.odb'
    honeycomb_odb_object = session.openOdb(name=honeycomb_odb_path)

    # The main session viewport must be set to the odb object using
    # the following line. If you might receive an error message that states
    # The current viewport is not associated with the output data base file
    # Request operation cancelled
    session.viewports['Viewport: 1'].setValues(displayedObject=honeycomb_odb_object)

    keyarray = session.odbData[honeycomb_odb_path].historyVariables.keys()
    theoutputvariablename = []

    for x in keyarray:
        if (x.find('ETOTAL') > -1):
            theoutputvariablename.append(x)

    # You may enter an entire path if you wish to store the result in
    # a particular file location.
    # One way to do it is using the following syntax.
    reportxy_name_and_path = reportxy_path + reportxy_name + '.txt'
    # Note however that folder MyNewFolder must exist otherwise
    # you will likely get the following error.
    # "I/O Error: C;/MyNewFolder:Directory Not Found"
    # You must either first create the script in windows before running
    # the script. Or if you wish to create it using Python commands you
    # must use the os.makedir() or os.makedirs() function
    # os.makedirs() is preferable because you can create multiple nested
    # directories in one statement if you wish.
    # Note that this function returns an exception if this directory already exists.
    # hence it is a good idea to use a try block

    try:
        os.makedirs(reportxy_path)
    except:
        print("Directory exists and hence no need to recreate it. Move on to next statement")

    print(reportxy_name_and_path)

    session.XYDataFromHistory(name='honeycombXYData-1', odb=honeycomb_odb_object, \
                              outputVariableName=theoutputvariablename[0])
    honeycomb_xydata_object = session.xyDataObjects['honeycombXYData-1']
    session.xyReportOptions.setValues(totals=ON, minMax=ON)
    session.writeXYReport(fileName=reportxy_name_and_path, \
                          xyData=(honeycomb_xydata_object,), appendMode=OFF)

    honeycomb_odb_object.close()

    # -----------------------------------------------------------------------------------------------------------
    # Read the displacement from the report
    extracted_line = ''
    # Need a boolean variable to state whether we are reading the
    # correct section of the file
    file_xydata_section = 0

    f = open(reportxy_name_and_path)
    for line in f:
        str = line
        if 'honeycombXYData-1' in str:
            file_xydata_section = 1
        if 'MAXIMUM' in str and file_xydata_section == 1:
            extracted_line = str
            extracted_list = extracted_line.split()
            Energy_absorbed = extracted_list[2]
            print("The Energy_absorbed in the plate is" + \
                  repr(Energy_absorbed))

    f.close()
    # ----------------------------------------------------------------------------------------------------------
    # Write this value as well as inputs to the output file

    file_output = open(reportxy_path + '\honeycombstructure_output.txt', 'w')
    file_output.write(repr(float(Energy_absorbed)) + "\n")

    honeycomb_odb_path = job_name + '.odb'
    honeycomb_odb_object = session.openOdb(name=honeycomb_odb_path)

    # The main session viewport must be set to the odb object using
    # the following line. If you might receive an error message that states
    # The current viewport is not associated with the output data base file
    # Request operation cancelled
    session.viewports['Viewport: 1'].setValues(displayedObject=honeycomb_odb_object)
    keyarray = session.odbData[honeycomb_odb_path].historyVariables.keys()
    theoutputvariablename = []

    for x in keyarray:
        if (x.find('U3') > -1):
            theoutputvariablename.append(x)

    # You may enter an entire path if you wish to store the result in
    # a particular file location.
    # One way to do it is using the following syntax.
    reportxy_name_and_path = reportxy_path + reportxy_name + '.txt'
    # Note however that folder MyNewFolder must exist otherwise
    # you will likely get the following error.
    # "I/O Error: C;/MyNewFolder:Directory Not Found"
    # You must either first create the script in windows before running
    # the script. Or if you wish to create it using Python commands you
    # must use the os.makedir() or os.makedirs() function
    # os.makedirs() is preferable because you can create multiple nested
    # directories in one statement if you wish.
    # Note that this function returns an exception if this directory already exists.
    # hence it is a good idea to use a try block

    try:
        os.makedirs(reportxy_path)
    except:
        print("Directory exists and hence no need to recreate it. Move on to next statement")

    print(reportxy_name_and_path)

    session.XYDataFromHistory(name='honeycombXYData-2', odb=honeycomb_odb_object, \
                              outputVariableName=theoutputvariablename[0])
    honeycomb_xydata_object = session.xyDataObjects['honeycombXYData-2']
    session.xyReportOptions.setValues(totals=ON, minMax=ON)
    session.writeXYReport(fileName=reportxy_name_and_path, \
                          xyData=(honeycomb_xydata_object,), appendMode=OFF)

    honeycomb_odb_object.close()

    # -----------------------------------------------------------------------------------------------------------
    # Read the displacement from the report
    extracted_line = ''
    # Need a boolean variable to state whether we are reading the
    # correct section of the file
    file_xydata_section = 0

    f = open(reportxy_name_and_path)
    for line in f:
        str = line
        if 'honeycombXYData-2' in str:
            file_xydata_section = 1
        if 'MINIMUM' in str and file_xydata_section == 1:
            extracted_line = str
            extracted_list = extracted_line.split()
            max_displacement = extracted_list[2]
            print("The Displacement of of the node at end of beam is" + \
                  repr((-1) * float(max_displacement)))

    f.close()
    # ----------------------------------------------------------------------------------------------------------
    # Write this value as well as inputs to the output file

    file_output = open(reportxy_path + '\honeycombstructure_output_1.txt', 'w')
    file_output.write(repr((-1) * float(max_displacement)) + "\n")

    honeycomb_odb_path = job_name + '.odb'
    honeycomb_odb_object = session.openOdb(name=honeycomb_odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=honeycomb_odb_object)
    keyarray = session.odbData[honeycomb_odb_path].historyVariables.keys()
    theoutputvariablename = []

    for x in keyarray:
        if (x.find('CFT3') > -1):
            theoutputvariablename.append(x)

    # You may enter an entire path if you wish to store the result in
    # a particular file location.
    # One way to do it is using the following syntax.
    reportxy_name_and_path = reportxy_path + reportxy_name + '.txt'
    # Note however that folder MyNewFolder must exist otherwise
    # you will likely get the following error.
    # "I/O Error: C;/MyNewFolder:Directory Not Found"
    # You must either first create the script in windows before running
    # the script. Or if you wish to create it using Python commands you
    # must use the os.makedir() or os.makedirs() function
    # os.makedirs() is preferable because you can create multiple nested
    # directories in one statement if you wish.
    # Note that this function returns an exception if this directory already exists.
    # hence it is a good idea to use a try block

    try:
        os.makedirs(reportxy_path)
    except:
        print("Directory exists and hence no need to recreate it. Move on to next statement")

    print(reportxy_name_and_path)

    session.XYDataFromHistory(name='honeycombXYData-3', odb=honeycomb_odb_object, \
                              outputVariableName=theoutputvariablename[0])
    honeycomb_xydata_object = session.xyDataObjects['honeycombXYData-3']
    session.xyReportOptions.setValues(totals=ON, minMax=ON)
    session.writeXYReport(fileName=reportxy_name_and_path, \
                          xyData=(honeycomb_xydata_object,), appendMode=OFF)

    honeycomb_odb_object.close()

    # -----------------------------------------------------------------------------------------------------------
    # Read the displacement from the report
    extracted_line = ''
    # Need a boolean variable to state whether we are reading the
    # correct section of the file
    file_xydata_section = 0

    f = open(reportxy_name_and_path)
    for line in f:
        str = line
        if 'honeycombXYData-3' in str:
            file_xydata_section = 1
        if 'MAXIMUM' in str and file_xydata_section == 1:
            extracted_line = str
            extracted_list = extracted_line.split()
            peak_force = extracted_list[2]
            print("The Peak Force is" + \
                  repr(float(peak_force)))

    f.close()
    # ----------------------------------------------------------------------------------------------------------
    # Write this value as well as inputs to the output file

    file_output = open(reportxy_path + '\honeycombstructure_output_2.txt', 'w')
    file_output.write(repr(float(peak_force)) + "\n")
    file_output.close()

    iteration_count = iteration_count + 1
    # -----------------------------------------------------------------------------------------------------------

input_file.close()



















