import numpy as np
import xml.etree.ElementTree as ET
import math
from shapely.geometry import Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import cascaded_union
from shapely.geometry.polygon import orient as orient

np.set_printoptions(threshold='inf')


# define namespace

ET.register_namespace('core','http://www.opengis.net/citygml/1.0') 
ET.register_namespace('bldg', 'http://www.opengis.net/citygml/building/1.0')
ET.register_namespace('gml', 'http://www.opengis.net/gml')
ET.register_namespace('grp', "http://www.opengis.net/citygml/cityobjectgroup/1.0")
ET.register_namespace('app', "http://www.opengis.net/citygml/appearance/1.0")
ET.register_namespace('xAL', "urn:oasis:names:tc:ciq:xsdschema:xAL:2.0")
ET.register_namespace('gen', "http://www.opengis.net/citygml/generics/1.0")
ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")

ns = {'bldg': 'http://www.opengis.net/citygml/building/1.0', 'core':'http://www.opengis.net/citygml/1.0', 'gml':'http://www.opengis.net/gml'}




# INPUT PATH
inp = '2018-3611_LoD1_gesamt_GEB3D_307450_5679890_cgml.gml'
# OUTPUT PATH
outp = '2018-3611_LoD1_gesamt_GEB3D_307450_5679890_cgml_simple.gml'


################################
###         classes         ####
################################



# a building is defined by a horizontal ground surface and the roof z-coord
class build_simple:

	def __init__(self,ground, zground, zroof):
		self.ground = ground
		self.zroof = zroof
		self.zground = zground



################################
###         routines        ####
################################

# this routine transforms every LOD1 building consisting of building parts into a single building with simplified shape
def reconstr_build(build,build_simple):

	lod1solid = build.find('bldg:lod1Solid', ns)
	
	# remove envelope
	boundedby = build.find('gml:boundedBy',ns)
	build.remove(boundedby)
	# remove LOD1 solid
	build.remove(lod1solid)


	# reconstruct simplified gml-building
	wallsurfaces = []
        ground = build_simple.ground
	xvals = []
	yvals = []
	for coord in ground.boundary.coords[:]:
	        xvals.append(coord[0])
		yvals.append(coord[1])

	# print xvals
        zground = build_simple.zground
        zroof = build_simple.zroof
	
	# make ground surface
        boundedby_new = ET.SubElement(build,'bldg:boundedBy',ns) 
        groundsurface_new = ET.SubElement(boundedby_new,'bldg:GroundSurface',ns)
        lod2multisurface = ET.SubElement(groundsurface_new,'bldg:lod2MultiSurface',ns)
        multisurface = ET.SubElement(lod2multisurface,'gml:MultiSurface',ns)
	surfmember = ET.SubElement(multisurface,'gml:surfaceMember',ns)
	polygon  = ET.SubElement(surfmember,'gml:Polygon',ns)
	exterior_new  = ET.SubElement(polygon,'gml:exterior',ns)		
	linearring  = ET.SubElement(exterior_new,'gml:LinearRing',ns)
	poslist  = ET.SubElement(linearring,'gml:posList',ns, srsDimension="3")


	text = ""

	for i in xrange(0,len(xvals),1):
		text += str(xvals[i]) + " "
              	text += str(yvals[i]) + " "
               	text += str(zground) + " "
	poslist.text = text[:-1] 

	# make roof surface
       	boundedby_new = ET.SubElement(build,'bldg:boundedBy',ns)     
       	roofsurface_new = ET.SubElement(boundedby_new,'bldg:RoofSurface',ns)
       	lod2multisurface = ET.SubElement(roofsurface_new,'bldg:lod2MultiSurface',ns)
        multisurface = ET.SubElement(lod2multisurface,'gml:MultiSurface',ns)
        surfmember = ET.SubElement(multisurface,'gml:surfaceMember',ns)
        polygon  = ET.SubElement(surfmember,'gml:Polygon',ns)
        exterior_new  = ET.SubElement(polygon,'gml:exterior',ns)
        linearring  = ET.SubElement(exterior_new,'gml:LinearRing',ns)
        poslist  = ET.SubElement(linearring,'gml:posList',ns, srsDimension="3")
        text = ""

        for i in xrange(0,len(xvals)-1,1):
        	text += str(xvals[i]) + " "
                text += str(yvals[i]) + " "      
                text += str(zroof)+ " "
        poslist.text = text[:-1]

	# add the walls
	for i in xrange(0,len(xvals)-1,1):
                boundedby_new = ET.SubElement(build,'bldg:boundedBy',ns)
       	        wallsurface_new = ET.SubElement(boundedby_new,'bldg:WallSurface',ns)	                
                lod2multisurface = ET.SubElement(wallsurface_new,'bldg:lod2MultiSurface',ns)
                multisurface = ET.SubElement(lod2multisurface,'gml:MultiSurface',ns)
                surfmember = ET.SubElement(multisurface,'gml:surfaceMember',ns)
         	polygon  = ET.SubElement(surfmember,'gml:Polygon',ns)
                exterior_new  = ET.SubElement(polygon,'gml:exterior',ns)
                linearring  = ET.SubElement(exterior_new,'gml:LinearRing',ns)
                poslist  = ET.SubElement(linearring,'gml:posList',ns, srsDimension="3")
       	        text = ""
		text +=str(xvals[i]) +" "
                text +=str(yvals[i]) +" " 
                text +=str(zground) +" "
                text +=str(xvals[i+1]) +" "
                text +=str(yvals[i+1]) +" "
                text +=str(zground) +" "
                text +=str(xvals[i+1]) +" "
                text +=str(yvals[i+1]) +" "
                text +=str(zroof) +" "
                text +=str(xvals[i]) +" "
                text +=str(yvals[i]) +" "
                text +=str(zroof) +" "
                text +=str(xvals[i]) +" "
                text +=str(yvals[i]) +" "
                text +=str(zground)
		poslist.text = text

	# add new envelope
	boundedby_new = ET.SubElement(build,'gml:boundedBy', ns)
	envelope = ET.SubElement(boundedby_new,'gml:Envelope', ns, srsDimension="3", srsName="urn:ogc:def:crs,crs:EPSG:6.12:25833,crs:EPSG:6.12:5783")
	lowercorner = ET.SubElement(envelope, 'gml:lowerCorner', ns)
	lowercorner.text = str(min(xvals))+" "+str(min(yvals))+" "+str(zground)
	uppercorner = ET.SubElement(envelope, 'gml:upperCorner', ns)
        uppercorner.text = str(max(xvals))+" "+str(max(yvals))+" "+str(zroof)
        stripNs(build)


# return all horizontal surfaces of the original LOD1 buildingpart
def horsrf(buildpart):
        lod1solid = buildpart.find('bldg:lod1Solid', ns)
        if(not lod1solid ==None):	
                solid  = lod1solid.find('gml:Solid',ns)
                exterior= solid.find('gml:exterior',ns)
                compositesurface = exterior.find('gml:CompositeSurface',ns)
		heights = []
		surfaces= []
		for surfacemember in list(compositesurface.findall('gml:surfaceMember', ns)):
			# check if surface is planar and if yes return z coordinate
			hght = check_horplan(surfacemember)	
			if not hght==-999:
				surfaces.append(surfacemember)
		return surfaces
	else:
		return None

# checks if a surface is horizontal and planar
def check_horplan(surface):
	polygon = surface.find('gml:Polygon',ns)
	exterior = polygon.find('gml:exterior',ns)
	linearring= exterior.find('gml:LinearRing',ns)
	poslist = str(linearring[0].text)
        coords =  list(map(float, poslist.split(' ')))
	zvals = coords[2::3]
	size1 = np.unique(zvals).size	
	if size1==1:		
		return zvals[0]
	else:	
		return -999

# this routine gets all the horizontal surfaces of the building parts and reconstructs the common ground surface and the roof height
def ground_roof(surfacelist):
	polygons = []
	nsurf = len(surfacelist)
	zvals = []
	for surface in surfacelist:
		polygon = surface.find('gml:Polygon',ns)
        	exterior = polygon.find('gml:exterior',ns)
        	linearring= exterior.find('gml:LinearRing',ns)
        	poslist = str(linearring[0].text)
        	coords =  list(map(float, poslist.split(' ')))
		x  = coords[0::3]
		y = coords[1::3]

		tuples = []
		for i in xrange(0, len(x),1):
			tuples.append(tuple([x[i],y[i]]))

		#convert into shapely object
		shape_pol = Polygon(tuples)
		polygons.append(shape_pol)
		zvals.append(coords[2])

	# if there are more than 2 horizontal surfaces we have multiple building parts 
	if len(surfacelist)>2:
		# inspecting the data, we know, that every second horizontal surface is a ground surface
		groundsrfs  = (polygons[1::2])
		roofsrfs  = (polygons[0::2])
		zroofs = zvals[0::2]
		zgrounds = zvals[1::2]
		zground = min(zgrounds)
		# calculate roof height
		areas = []
		i = 0
		for roofsrf in list(roofsrfs):
			area = roofsrf.area
			# check if a groundsurface has the same height, then subtract it from the roof surface
			ind_roof = np.argwhere(np.array(zgrounds) == zroofs[i])


			if len(ind_roof):	
				area = max(area-groundsrfs[ind_roof[0][0]].area,0)
				zgrounds.remove(zgrounds[ind_roof[0][0]])
				groundsrfs.remove(groundsrfs[ind_roof[0][0]])

			areas.append(area)	
			i += 1	
		sum_areas= sum(np.array(areas))

		if not sum_areas==0:
			zroof=  np.dot(np.array(zroofs),np.array(areas))/sum_areas
		else:
			zroof= max(zroofs)	

		#use shapely function to merge all ground surfaces to  unified surface
		ground_uni = cascaded_union(polygons)
		# simplify polygon with a tolerance of 0.5m
		ground_uni = ground_uni.simplify(0.5)
                if not ground_uni.is_valid:
                        ground_uni = ground_uni.buffer(0)
		# in case not all surfaces could be merged, redo after simplification
		ground_uni = cascaded_union(ground_uni)

		# check if resulting building is larger than 125m^3 and higher than 5m   

		if ground_uni.area*(zroof-zground)>125 and zroof-zground>5:
                        coordslist = []
			# in case not all surfaces could be merged, take the largest as representative
			boundaries = ground_uni.boundary
        	        if isinstance(boundaries, MultiLineString):
				areas = []
				for member in boundaries: 
					areas.append(Polygon(member).area)
				ind = areas.index(max(areas))				
				# ensure ground surface is counter clockwise orientated
				new_ground = orient(Polygon(boundaries[ind]),1)     
				new_build = build_simple(new_ground, zground, zroof)

			else:	
				new_ground = orient(ground_uni,1)
                                new_build = build_simple(new_ground, zground, zroof)

			return new_build

		
		else:
			return None
	else:	
		# check if resulting building is larger than 125m^3 and higher than 5m	
                if polygons[0].area*(zvals[0]-zvals[1])>125 and zvals[0]-zvals[1]>5:
			new_ground = orient(polygons[0],1)
			new_build = build_simple(new_ground, zvals[1], zvals[0])

			return new_build

		else:
		
			return None




# routine to strip redundant namespaces
def stripNs(el):
        '''Recursively search this element tree, removing namespaces.'''
        el.attrib.pop('core',None)
        el.attrib.pop('bldg',None)
        el.attrib.pop('gml',None)
        for child in el:
                stripNs(child)


# routine to make a nice structure
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def make_model():
        k = 0

        tree = ET.parse(inp)

        root = tree.getroot()
        print "file parsed"
        for cityobjectmember in root.findall('core:cityObjectMember',ns):
		k +=1
                build = cityobjectmember.find('bldg:Building',ns)
		surfaces = []
		for consistsofbuildingparts in list(build.findall('bldg:consistsOfBuildingPart',ns)):
			buildpart = consistsofbuildingparts.find('bldg:BuildingPart', ns)
			surf = horsrf(buildpart)
	                if not surf == None:
				surfaces.extend(surf)
			build.remove(consistsofbuildingparts)

		surf = horsrf(build)
                if not surf == None:
			surfaces.extend(surf)
		print "Modify building {}".format(k)

		if len(surfaces):		
			new_build = ground_roof(surfaces) 

			if not new_build==None:
				reconstr_build(build, new_build)
			else:
				root.remove(cityobjectmember)
			del new_build
		else:
			root.remove(cityobjectmember)

	indent(root)
	new_tree = ET.ElementTree(root)
	new_tree.write(outp)



def main():

         make_model()




if __name__ =="__main__":
        main()

# END 

