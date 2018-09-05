#########  ###  ##########   ##  ###  ######  #########  #####   ##  ####  ####   ######  ################### ###### ########
#########  ###  ####  ####   ## ###  #######  ###    ##  ######  ##  ####  ###   #######  ###   ############# ###    ###  ###
###                          #####  ###       ###    ##  ###  ## ##   ###  ##   ##        ###  ####    ###   ######  ########
###  ##    ###     ###  ##         ###        ###    ##  ###  #####   ###  #   #######    ##           ###  #######  #######
###  ##    ###     ###  ##  ###                                                ######     ##   ###     ###  ###      ##   ##
###        ###     ###     ###    ###         ###    ##  ###  ######   #####   ###             ###     ###  ###      ##   ##
#########  ###     ###    ###     #### ####   ###    ##  ###   #####   #####   ###        ##    ###    ###  #######  ##    ##
#########  ###     ###   ###  ##  #### #####  #########  ###    ####   #####   ########  ##      ###   ###  #######  ##     ##

import sys
sys.path.append('/vols/fs1/work/weger/PROG/PYTHON/lib/lib/python2.7/site-packages/')
import numpy as np
import xml.etree.ElementTree as ET
from rotkoord import Rotgrid
import rotkoord
import math
from pyproj import Proj, transform
from shapely.geometry import Polygon, LineString, MultiLineString, MultiPolygon
import shapely.ops as sp 
from netCDF4 import Dataset
import time

np.set_printoptions(threshold='inf')

#define namespace to read city gml data

ns = {'bldg': 'http://www.opengis.net/citygml/building/1.0', 'core':'http://www.opengis.net/citygml/1.0', 'gml':'http://www.opengis.net/gml'}


####################################################
########          !!!!settings!!!!            ######
####################################################


citygml_path = '2018-3611_LoD1_gesamt_GEB3D_307450_5679890_cgml_simple.gml'
impsurf_path = '/vols/fs1/work/weger/URBAN/CITYGML/Leipzig_in/imp_surf_lpz_rot0d0018.nc'
output_path= '/vols/fs1/work/weger/URBAN/CITYGML/Leipzig_out/leipzig_stats.nc'

#define grid: 
xpole = -170
ypole = 40
xfirst = 1.41
yfirst = 1.30
nx = 80
ny = 75
xinc = 0.0018
yinc = 0.0018

# dont change this!
lats = np.linspace(yfirst, yfirst+ (ny-1)*yinc, ny , dtype=float)
lons = np.linspace(xfirst, xfirst+ (nx-1)*xinc, nx , dtype=float)



# set projection of the city-gml data
Proj1 = Proj(init='epsg:25833')

# dont change this!
Proj2 = Proj(init='epsg:4326')

# set urban half levels
hhlu = [0,5,10,15,20,25,30,40,70,120,170]


#set street directions
street_dir = [-45,0,45,90]

# set radius within walls are considered for street-width calculation 
radius_check =  150
# set minimum distance for walls 
dist_min = 2

##############################################
######           classes               #######
##############################################


# a grid cell is approximated by a rectangular polygon
class cell:
        def __init__(self, polygon):
                self.polygon = polygon



# a wall surface is rectangular and vertical, it has several attributes:
# id_bldg -> id of building containing the wall
# id_wall -> id of wall in the list of walls of a building
# pos -> geometric centroid
# surfinter -> ground surface intersection 
# area -> surface area
# height -> roof height
# normal -> orientation vector perpendicular to wall and pointing outwards
# proj -> index of street direction, for which projection of normal vector is maximal 
# distance -> mean distance to other visible wall surfaces 
 
class wallsurface:

# at initialization, properties are calculated
	def __init__(self,points, id_bldg, id_wall):
		self.id_bldg = id_bldg
		self.id_wall = id_wall
		surfpt1 = points[0]
		surfpt2 = points[1]
		self.pos = [(surfpt1[0]+surfpt2[0])/2., (surfpt1[1]+surfpt2[1])/2.]
		self.surfinter = surfpt1[0:2],surfpt2[0:2]
		self.height = points[3][2]-points[0][2]
		self.area = LineString([surfpt1[0:2],surfpt2[0:2]]).length*self.height 	
		self.dist = 0 
		dotprod = []
                surf_normal = [surfpt2[1]-surfpt1[1],surfpt1[0]-surfpt2[0]]
		surf_normal_lnght = np.sqrt(np.power(surfpt2[1]-surfpt1[1],2)+np.power(surfpt1[0]-surfpt2[0],2))
		surf_normal = surf_normal[0]/surf_normal_lnght,  surf_normal[1]/surf_normal_lnght
		self.normal = surf_normal	
                counter = 0
                for direction in street_dir:	
                        vector = np.array([np.cos(np.radians(direction)), np.sin(np.radians(direction))])
			dotprod.append(abs(np.dot(vector, surf_normal)))                      
		self.proj = dotprod.index(max(dotprod))

	# function to modify area in case of overlapping walls
	def set_area(self,area_new):
		self.area = area_new	

	# function to set distance after initialization
	def set_distance(self,dist):
		self.dist = dist

# a building consists of all  surfaces, furthermore it has following attributes:

# id_bldg -> id of building in the city model
# ground -> ground surface as a polygon
# center -> centroid projected on the ground
# next_blgs -> a sorted list of neighboring building indices, beginning with the closest

class buildg:
	def __init__(self,polygon,walls, id_bldg):
		self.id_bldg = id_bldg
		self.walls = list(walls)
		self.next_blgs = []
		self.ground = polygon
		self.center = polygon.centroid.coords[0][:]
	# remove a wall surface from a building in case of maximal overlap		
	def remove(self, id_wall):
		walls = self.walls
		for wall in walls:
			if wall.id_wall == id_wall:
				walls.remove(wall)
				break
		self.walls = walls 

	# get a wall by its id
	def get(self, id_wall):
		walls = self.walls
                for wall in walls:

                        if wall.id_wall == id_wall:
                             	return wall 
                           

	# add index of a nearby building
	def addnext(self, index):
		next_blgs = self.next_blgs
		next_blgs.append(index)		
		self.next_blgs = next_blgs



############################################
###             routines                ####
############################################
			
# open citygml file and extract all buildings, and initialize building and wall class objects from  above	
def extract_buildings(citygml_path):
        k = 0
	print "parse file"
	tree = ET.parse(citygml_path)
        root = tree.getroot()
        print "file parsed"
	bldg_lst = []
	counter = 0
        for cityobjectmember in list(root.findall('core:cityObjectMember',ns))[0:500]:
		print "extract building " + str(counter)
                gmlbuilding = cityobjectmember.find('bldg:Building',ns)
		gmlwallsurfaces = gmlbuilding.findall('bldg:boundedBy/bldg:WallSurface',ns)
		walls = []
		areas = []
		id_wall = 0
		for gmlwallsurface in gmlwallsurfaces:
			lod2multisurface = gmlwallsurface.find('bldg:lod2MultiSurface',ns)
			multisurface = lod2multisurface.find('gml:MultiSurface',ns)
			surfmember = multisurface.find('gml:surfaceMember',ns)
			polygon = surfmember.find('gml:Polygon',ns)
			exterior = polygon.find('gml:exterior',ns)
			linearring = exterior.find('gml:LinearRing',ns)
			poslist = linearring.find('gml:posList',ns).text
			coords = list(map(float,poslist.split(' ')[:-1]))
			xvals = coords[0::3]
			yvals = coords[1::3]
			zvals = coords[2::3]
			zvals.append(zvals[0])		
			points = []	
       		        for i in xrange(0,len(xvals),1):
                       		points.append(tuple([xvals[i],yvals[i], zvals[i]]))
               		walls.append(wallsurface(points, counter, id_wall))
			id_wall +=1
		# sort walls based on surface area
		for wall in walls:
			areas.append(wall.area)
		indices = np.argsort(np.array(areas), kind= 'mergesort')
		total_area = sum(areas)
		walls_sorted = []
		area_accum = 0
		for index in reversed(indices):
			walls_sorted.append(walls[index])
			area_accum+=walls[index].area
			# remove very small walls such that total surface remains >= 99%
			if area_accum/total_area>0.99:
				break	

		gmlgroundsurface = gmlbuilding.find('bldg:boundedBy/bldg:GroundSurface',ns)
		lod2multisurface = gmlgroundsurface.find('bldg:lod2MultiSurface',ns)
		multisurface = lod2multisurface.find('gml:MultiSurface',ns)
		surfmember = multisurface.find('gml:surfaceMember',ns)
		polygon = surfmember.find('gml:Polygon',ns)
		exterior = polygon.find('gml:exterior',ns)
		linearring = exterior.find('gml:LinearRing',ns)
		poslist = linearring.find('gml:posList',ns).text
		coords = list(map(float,poslist.split(' ')[:-1]))
		xvals = coords[0::3]
		yvals = coords[1::3]
		points = []
		for i in xrange(0,len(xvals),1):
			points.append(tuple([xvals[i],yvals[i]]))
		ground = Polygon(points)
		bldg_lst.append(buildg(ground,list(walls_sorted), counter))
		counter +=1
	return bldg_lst


# find all neighboring buildings around a building
def connect_buildings(bldg_lst):
        nbuild = len(bldg_lst)
        print nbuild
        counter = 0
	# find all buildings within check_radius

        for building  in bldg_lst:
                print "connect building " + str(counter)
                center1 = building.center
                dist = []
                args = []
		# find closest buildings

                i = counter
		# the list is already sorted, so it is enough to consider only 1000 buildings in the list around the index
                for i in xrange(max(counter-500,0),min(counter+500,nbuild-1),1):
                        building2 = bldg_lst[i]
                        args.append(i)
			# calculate the distance to all buildings
                        dist.append(np.sqrt(pow(building2.center[0]-center1[0],2)+pow(building2.center[1]-center1[1],2)))
                dist = np.array(dist)
		# sort the list
                args_sorted = np.argsort(dist, kind='mergesort')
		# save the first 49 closest building indices 
                for arg in args_sorted[1:50]:	
                        building.addnext(args[arg])
                counter +=1


#calculate the street width
def calc_dist(bldg_lst):
	counter = 0
	# just a variable to save the mimimum distance in the urban calculation
	min_dist = 1E6  
	remove_walls = []
	for building in bldg_lst:
		print "calculate distance for walls of building " + str(counter)
		#get all neighboring buildings
		neighbors = building.next_blgs
		walls = building.walls
		for wall in walls:
			remove_wall_status = False
			# get all potential visible wall surfaces, which are orientated towards the wall 
			# (i.e. dot product of neighbor surface normal vector and position vector must be negative)
			# and located in the positive half space (i.e. dot product of normal vector and position vector must be postive)
			pot_viswalls = []
			dists = []
			pos_vects = []
			for index in neighbors:
				neighbor = bldg_lst[index]
				for neighborwall in neighbor.walls:
					# the position vector points at the neighboring wall
					pos_vec = [neighborwall.pos[0]-wall.pos[0], neighborwall.pos[1]-wall.pos[1]]
					if (np.dot(pos_vec, neighborwall.normal)<0 and np.dot(wall.normal, pos_vec)>0):
						pot_viswalls.append(neighborwall)
						# effective distance is length of pos_vec projected on normal vector of wall
	        	                        dists.append(np.dot(wall.normal, pos_vec))	
						pos_vects.append(pos_vec)	
			# sort all distances to start with the shortest distance
			args = np.argsort(np.array(dists), kind = 'mergesort')
			# construct a MultiLineString consisting of all pot_viswalls
			lines = []
			if len(pot_viswalls)>0:
				for pot_viswall in pot_viswalls:
					lines.append(pot_viswall.surfinter)
				multiline = MultiLineString(lines)
				areasum_eff = 0
				dist_wght = 0
				for i in xrange(0,len(pot_viswalls),1):
					pot_viswall = pot_viswalls[args[i]]
					# for every potential visible wall, construct a connecting line
					#add 0.01 of the normal vectors on both sides in order to avoid intersection with the two target walls
					st_point = tuple([pot_viswall.pos[0]+pot_viswall.normal[0]*0.01,pot_viswall.pos[1]+wall.normal[1]*0.01])
					end_point = tuple([wall.pos[0]+wall.normal[0]*0.01,wall.pos[1]+wall.normal[1]*0.01])
					connection = LineString([st_point,end_point])
					if not connection.intersects(multiline):
						# we have a wall which is "visible"		
						vis_wall = pot_viswall

						# if distance < dist_min one wall is obscured by the other wall. 
						# In this case, the smaller wall is subtracted from the larger one and the smaller one removed completely
						if dists[args[i]] < dist_min:
							area1 = wall.area
							area2 = pot_viswall.area
							if area1 > area2:
								wall.set_area(area1-area2)
								id_bldg = vis_wall.id_bldg
								id_wall = vis_wall.id_wall	
								bldg_lst[id_bldg].remove(id_wall)	
								vis_wall.set_area(0)	
								continue
							else:		
								pot_viswall.set_area(area2-area1)
								# the wall is completely obscured, in this case exit and continue with next wall
								id_wall = vis_wall.id_wall
								id_bldg = vis_wall.id_bldg
								bldg_lst[id_bldg].get(id_wall).set_area(area2-area1)
								id_wall = wall.id_wall
                                                                id_bldg = wall.id_bldg	
								remove_walls.append([id_bldg,id_wall])
								remove_wall_status = True
								break

						# this is a valid distance and total accumulated surface area of visible walls 
						# has not exceeded yet wall surface area
						if  dists[args[i]] < radius_check:
							if areasum_eff< wall.area:
								fac = (abs(np.dot(pos_vects[args[i]],vis_wall.normal))*
									abs(np.dot(pos_vects[args[i]],wall.normal))/
									pow(np.linalg.norm(pos_vects[args[i]]),2))
								area_eff =  vis_wall.area*fac
								dist_wght += dists[args[i]]* area_eff
								areasum_eff  += area_eff
							else:
							# exit loop and set the distance, if area of visible walls exceeds area of considered wall 	
								break	
				# if considered wall has no visible walls within check radius, then mark wall for removal
				if remove_wall_status:	
					continue
                                if areasum_eff< 1E-3:
                                	id_wall = wall.id_wall
                                        id_bldg = wall.id_bldg		
					remove_walls.append([id_bldg,id_wall])
					continue
				else:
					wall.set_distance(dist_wght/areasum_eff)
					min_dist = min(min_dist, dist_wght/areasum_eff)                              
					continue 
			# if considered wall has no potential visible walls within check radius, then mark wall for removal
			else:
	                        id_wall = wall.id_wall
                                id_bldg = wall.id_bldg
				remove_walls.append([id_bldg,id_wall])
		counter +=1	
	# remove all marked walls
	for ids in remove_walls:
		id_bldg = ids[0]
		id_wall = ids[1]
		bldg_lst[id_bldg].remove(id_wall)
	# return the minumum distance of all walls in the city model (should be hardly above dist_min)
	return min_dist


# transform rotated latlon grid to projection of city-gml data 
def transform1(point):
        mapping = Rotgrid(xpole, ypole,polerotate=10)
        point= mapping.transform(point[0], point[1], inverse=True)
        point = transform(Proj2, Proj1,point[0],point[1])
        return point



# transform projection of city-gml data to rotated latlon
def transform2(point):
        point = transform(Proj1, Proj2,point[0],point[1])
        mapping = Rotgrid(xpole, ypole,polerotate=10)
        point= mapping.transform(point[0], point[1], inverse=False)
        return point



# transform the predefined latlon grid to projection of city-gml data using rectangular polygons 
def transform_grid():
        cells = np.empty([ny,nx],dtype=object)
        vcell = np.vectorize(cell)
        for i  in xrange(0,ny,1):
                for j in xrange(0,nx,1):
                        lowleft = transform1([lons[j]-xinc/2.,lats[i]-yinc/2.])
                        lowright = transform1([lons[j]+xinc/2.,lats[i]-yinc/2.])
                        topleft = transform1([lons[j]-xinc/2.,lats[i]+yinc/2.])
                        topright = transform1([lons[j]+xinc/2.,lats[i]+yinc/2.])
                        polygon =  Polygon([lowleft,lowright,topright,topleft, lowleft])
                        cells[i,j]= polygon
        return vcell(cells)


			
# project on grid and write output
def calc_urb_param(bldg_lst, grid):
	fr_street = np.empty([len(street_dir),ny,nx])
	street_w = np.empty([len(street_dir),ny,nx])
	buildprop = np.empty([len(street_dir),len(hhlu),ny,nx])
	fr_street.fill(0.)
	buildprop.fill(0.)
	street_w.fill(0.)
        building_fr = np.empty([ny,nx])
        building_fr.fill(0.0)
	wall_counter = 0
	zero_wall_counter = 0
	counter = 0
	for building in bldg_lst:
		print "write parameters for building " + str(counter)

		# street direction, building height probability and street width are accumulated "wall wise"
		# and not "building wise" for better results at high grid resolutions 
		for wall in building.walls:
			wall_counter += 1
			surfinter = wall.surfinter
			angle = wall.proj
			height = np.argmin(abs(np.array(hhlu)-wall.height))
			area = wall.area
			if wall.dist < 1E-3:
				zero_wall_counter += 1
			dist = wall.dist
			# project values on latlon grid			
			point1 =  transform2(surfinter[0])
			point2 =  transform2(surfinter[1])
			lat_min_ind = np.argmin(abs(lats-min(point1[1],point2[1])))
			lat_max_ind = np.argmin(abs(lats-max(point1[1],point2[1])))
			lon_min_ind = np.argmin(abs(lons-min(point1[0],point2[0])))
			lon_max_ind = np.argmin(abs(lons-max(point1[0],point2[0])))	
			for i in xrange(lat_min_ind, lat_max_ind+1,1):
				for j in xrange(lon_min_ind, lon_max_ind+1,1):
					# frac is the fraction of the wall, thats contained in the grid cell
					frac = (grid[i,j].polygon.intersection(LineString(surfinter))).length/LineString(surfinter).length	
					fr_street[angle,i,j] += area*frac
					buildprop[angle,height,i,j] += area*frac
					street_w[angle,i,j] += dist*area*frac
			area_tot =  fr_street

		# aggregate ground surfaces to calculate building fraction
		ground = building.ground
		# get the envelope containing the ground surface  and select all grid cells containing this envelope 
		envelope = ground.envelope.boundary.coords[:]
                lowleft =  transform2(envelope[0])
                lowright = transform2(envelope[1])
                topright = transform2(envelope[2])
                topleft =  transform2(envelope[3])
                lon_min = min(lowleft[0],topleft[0])
                lon_max = max(lowright[0],topright[0])
                lat_min = min(lowleft[1],lowright[1])
                lat_max = max(topleft[1],topright[1])
                lat_min_ind = np.argmin(abs(lats-lat_min))
                lat_max_ind = np.argmin(abs(lats-lat_max))
                lon_min_ind = np.argmin(abs(lons-lon_min))
                lon_max_ind = np.argmin(abs(lons-lon_max))

                for i in xrange(lat_min_ind,lat_max_ind+1,1):
                        for j in xrange(lon_min_ind,lon_max_ind+1,1):
				# calculate the intersection of each groundsurface with each grid cell in the selection
                                building_fr[i,j]+= (ground.intersection(grid[i,j].polygon).area)

		counter +=1
	#calculate normalized and area averaged quantities
	for i in xrange(0,ny,1):
		for j in xrange(0,nx,1):
			building_fr[i,j] = building_fr[i,j]/grid[i,j].polygon.area
			for k in xrange(0,len(street_dir),1):
				if fr_street[k,i,j] >0 :
					street_w[k,i,j] = street_w[k,i,j]/area_tot[k,i,j]
				norm = np.sum(buildprop[k,:,i,j])
                                if norm > 0:
                                        buildprop[k,:,i,j] = buildprop[k,:,i,j]/norm
					for l in xrange(0, len(hhlu),1):
						if buildprop[k,l,i,j] ==1:
							buildprop[k,l,i,j] -= 1E-10
						if buildprop[k,l,i,j] ==0:
							buildprop[k,l,i,j] += 1E-10
		        norm = np.sum(fr_street[:,i,j])
                        if norm>0:
                                fr_street[:,i,j] = fr_street[:,i,j]/norm

	# get urban fraction 
        file_impsurf = Dataset(impsurf_path,  format='NETCDF4') 
	imp_surf = file_impsurf.variables['Imp_surf'][:,:]
	imp_surf = imp_surf/100.
        build_w = np.empty([len(street_dir),ny,nx])
	build_w.fill(0.)	
	urbancl_fr = np.empty([ny,nx])
	urbancl_fr.fill(0.)
	
	# do some corrections and calculate building width
	for i in xrange(0,ny,1):
		for j in xrange(0,nx,1):
			for k in xrange(0, len(street_dir),1):
				if street_w[k,i,j] == 0:
					fr_street[k,i,j] = 0
			if np.sum(street_w[:,i,j]) ==0:
				building_fr[i,j] = 0
			if building_fr[i,j] > 0:
				urbancl_fr[i,j] = 1
			if imp_surf[i,j] < building_fr[i,j]:
				imp_surf[i,j] = min(building_fr[i,j]+0.01,1)
			if urbancl_fr[i,j]==0:
				imp_surf[i,j] = 0
			if building_fr[i,j]>0:
				build_w[:,i,j] = street_w[:,i,j]*building_fr[i,j]/(imp_surf[i,j]-building_fr[i,j])




	# write output to netcdf

        dataset_new = Dataset(output_path, 'w', format='NETCDF4')

        ncdim_rlat = dataset_new.createDimension('rlat', lats.size)
        ncdim_rlon = dataset_new.createDimension('rlon', lons.size)
	ncdim_nuc = dataset_new.createDimension('nuc', 1)
	ncdim_strdir = dataset_new.createDimension('streetdir', len(street_dir))
	ncdim_hghs = dataset_new.createDimension('uheight1', len(hhlu))

        ncvar_lat = dataset_new.createVariable('rlat', np.float, 'rlat')
        ncvar_lon = dataset_new.createVariable('rlon', np.float, 'rlon')
	ncvar_hghs = dataset_new.createVariable('uheight1', np.float, ('uheight1'))
	ncvar_strdir = dataset_new.createVariable('streetdir', np.float, ('streetdir'))
	ncvar_frstr = dataset_new.createVariable('FR_STREETD', np.float, ('nuc','streetdir','rlat', 'rlon'))
	ncvar_streetw = dataset_new.createVariable('STREET_W', np.float, ('nuc','streetdir','rlat', 'rlon'))
        ncvar_buildw = dataset_new.createVariable('BUILD_W', np.float, ('nuc','streetdir','rlat', 'rlon'))
        ncvar_bprop = dataset_new.createVariable('BUILD_PROP', np.float, ('nuc','streetdir','uheight1','rlat', 'rlon'))
        ncvar_frurban = dataset_new.createVariable('FR_URBAN', np.float, ('rlat', 'rlon'))
	ncvar_frurbcl  = dataset_new.createVariable('FR_URBANCL', np.float, ('nuc','rlat','rlon'))
	ncvar_frbuild = dataset_new.createVariable('FR_BUILD', np.float, ('nuc','rlat','rlon'))
	
        dataset_new.description = 'Leipzig building statistics'
        ncvar_lat.units = 'degree_north'
        ncvar_lon.units = 'degree_east'
        ncvar_lat[:]   = lats[:]
        ncvar_lon[:]  = lons[:]
	ncvar_hghs[:] = hhlu[:]
        ncvar_strdir[:] = street_dir[:]
        ncvar_frstr[0,:,:,:] = fr_street[:,:,:]
	ncvar_streetw[0,:,:,:]  = street_w[:,:,:] 
	ncvar_buildw[0,:,:,:] = build_w[:,:,:]
        ncvar_bprop[0,:,:,:,:] = buildprop[:,:,:,:]
	ncvar_frurban[:,:] = imp_surf[:,:]
	ncvar_frurbcl[0,:,:] = urbancl_fr[:,:]
	ncvar_frbuild[0,:,:] = building_fr[:,:]

	
	ncvar_lat.grid_mapping = "rotated_pole"
	ncvar_lon.grid_mapping = "rotated_pole"
	ncvar_hghs.grid_mapping = "rotated_pole"
	ncvar_strdir.grid_mapping = "rotated_pole"
        ncvar_frstr.grid_mapping = "rotated_pole"
        ncvar_streetw.grid_mapping = "rotated_pole"
        ncvar_buildw.grid_mapping = "rotated_pole"
        ncvar_bprop.grid_mapping = "rotated_pole"
        ncvar_frurban.grid_mapping = "rotated_pole"
        ncvar_frurbcl.grid_mapping = "rotated_pole"
	ncvar_frbuild.grid_mapping = "rotated_pole"

        dataset_new.close()
	


def main():
        start_time = time.time() 
	bldg_lst = extract_buildings(citygml_path)
        print "buildings extracted"
	connect_buildings(bldg_lst)	
        print "buildings connected"
	min_dist = calc_dist(bldg_lst)
	print "street width calculated"
	print "make grid"
        grid = transform_grid()
	urb_param = calc_urb_param(bldg_lst, grid)
	print "mimimum wall distance: " + str(min_dist)
        print "Urban calculation took {} minutes".format(int((time.time() - start_time)/60))
if __name__ =="__main__":
        main()
