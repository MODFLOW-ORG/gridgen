/*! \file MeshReader.h
\brief Tri mesh reader. Not implemented yet.

\author <a href="masc.cs.gmu.edu/">MASC group</a>, George Mason University,
        <a href="profile.usgs.gov/langevin/">Christian Langevin</a>, USGS
\bug    No known bugs.
*/

#pragma once

#include <string>
using namespace std;

#include <Point>
using namespace mathtool;


namespace cusg
{

/*! \brief an exp tri mesh reader

	Not implemented yet...
*/

class Mesh
{
public:

	/// Open the file, read the number of vertices and number of cells and
	/// allocate arrays.
	Mesh(const string & filename)
	{
	
	}
	
	
	void create_connectivity()
	{
	}
	
	double get_cell_area(int nodeid)
	{
	
	}
	
	Point2d get_cell_centroid(int nodeid)
	{
		Point2d o;
		return o;
	}
	
	Point2d get_cell_circumcenter(int nodeid)
	{
		Point2d o;
		return o;
	}
	
	void draw()
	{
	
	}
	
	void get_usg_csr_data(bool storediagonal=false)
	{
	
	}
};

}//end namespace cusg

/*
class Mesh(object):
    def __init__(self, filename):
        
        #Open the file, read the number of vertices and number of cells and
        #allocate arrays.
        f = open(filename, 'r')
        line = f.readline()
        linelist = line.split()
        ncells, nvertices = linelist[0:2]
        try:
            numcellparameters, numvertparameters = linelist[2:4]
        except:
            numcellparameters = 0
            numvertparameters = 0
        self.nvertices = int(nvertices)
        self.ncells = int(ncells)
        self.numcellparameters = int(numcellparameters)
        self.numvertparameters = int(numvertparameters)
        self.nodes = self.ncells
        print self.ncells, self.nvertices, self.numcellparameters, self.numvertparameters
        self.vertices = numpy.zeros((self.nvertices, 3), dtype=object)
        self.cellvertices = numpy.zeros( (self.ncells), dtype=object)
        if numcellparameters > 0:
            self.cellparameters = numpy.empty( (self.ncells, self.numcellparameters),
                dtype=float)
        if numvertparameters > 0:
            self.vertparameters = numpy.empty( (self.nvertices, self.numvertparameters),
                dtype=float)

        #read the vertices
        f.readline()
        for ivert in range(self.nvertices):
            line = f.readline()
            linelist = line.split()
            c, iv, x, y = linelist[0:4]
            self.vertices[ivert, :] = numpy.array([float(x), float(y), 0])
            if self.numvertparameters > 0:
                self.vertparameters[ivert, :] = numpy.array(linelist[4:])

        #read the cell information
        for icell in range(self.ncells):
            line = f.readline()
            linelist = line.split()
            iverts = []
            for ic in linelist[2:5]:
                iverts.append(int(ic) - 1) #convert to zero index
            self.cellvertices[icell] = iverts
            if self.numcellparameters > 0:
                self.cellparameters[icell, :] = numpy.array(linelist[5:])

        #close file and create connectivity       
        f.close()
        self._create_connectivity()
        
        return        

    def _create_connectivity(self):
        import node
        self.nodegroup = node.NodeGroup(self.ncells, chunksize=1)

        #create a face dictionary that maps the face to the connected cells        
        facedict = {}
        self.connectionlist = numpy.empty( (self.ncells), dtype=object)
        for i in range(self.ncells):
            self.connectionlist[i] = []
        
        for nodeid in range(self.ncells):
            nverts = len(self.cellvertices[nodeid])
            for iface in range(nverts):
                #create the face, which is a list of two nodes (in 2d)
                idx = iface + 1
                if iface == nverts - 1:
                    idx = 0
                face = [ self.cellvertices[nodeid][iface], self.cellvertices[nodeid][idx] ]
                face.sort()
                face = tuple(face)
                #store this face and its nodeid to facedict
                if facedict.has_key(face):
                    cell_list = facedict[face]
                    cell_list.append(nodeid)
                    facedict[face] = cell_list
                    #if isinstance(self.connectionlist[nodeid], list):
                    #    self.connectionlist[nodeid].append( [cell_list[0], face] )
                else:
                    facedict[face] = [nodeid]
                    #self.connectionlist[nodeid] = []
            
        self.nja = 0
        for face, c in facedict.iteritems():
            if len(c) > 1:
                self.connectionlist[c[0]].append( [c[1], face] )
                self.connectionlist[c[1]].append( [c[0], face] )
                self.nja += 2
            else:
                #this is a boundary face
                pass
           #meshcell = MeshCell(nodeid, connections)
        return

    def get_cell_area(self, nodeid):
        nverts = len(self.cellvertices[nodeid])
        a = 0.
        #create a copy of the cell vertices and close the polygon
        cellvertices = list(self.cellvertices[nodeid])
        cellvertices.append(self.cellvertices[nodeid][0])
        for i in range(nverts):
            iv = cellvertices[i]
            iv1 = cellvertices[i + 1]
            x = self.vertices[iv, 0]
            y = self.vertices[iv, 1]
            xp1 = self.vertices[iv1, 0]
            yp1 = self.vertices[iv1, 1]
            a += (x * yp1 - xp1 * y)
        a = a * 0.5
        return a
        
    def get_cell_centroid(self, nodeid):
        nverts = len(self.cellvertices[nodeid])
        cx = 0.
        cy = 0.
        #create a copy of the cell vertices and close the polygon
        cellvertices = list(self.cellvertices[nodeid])
        cellvertices.append(self.cellvertices[nodeid][0])
        for i in range(nverts):
            iv = cellvertices[i]
            iv1 = cellvertices[i + 1]
            x = self.vertices[iv, 0]
            y = self.vertices[iv, 1]
            xp1 = self.vertices[iv1, 0]
            yp1 = self.vertices[iv1, 1]
            cx += (x + xp1) * (x * yp1 - xp1 * y)
            cy += (y + yp1) * (x * yp1 - xp1 * y)
        a = self.get_cell_area(nodeid)
        cx = cx * 1./ 6. / a
        cy = cy * 1./ 6. / a
        return (cx, cy)        

    def get_cell_circumcenter(self, nodeid):
        '''
        Calculate cell circumcenter.  Cell must be a triangle.
        '''
        assert len(self.cellvertices[nodeid]) == 3
        ia = self.cellvertices[nodeid][0]
        ib = self.cellvertices[nodeid][1]
        ic = self.cellvertices[nodeid][2]
        ax = self.vertices[ia, 0]
        ay = self.vertices[ia, 1]
        bx = self.vertices[ib, 0]
        by = self.vertices[ib, 1]
        cx = self.vertices[ic, 0]
        cy = self.vertices[ic, 1]
        
        d = 2. * ( ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        Cx = ((ax ** 2 + ay ** 2) * (by - cy) + 
            (bx ** 2 + by ** 2) * (cy - ay) +
            (cx ** 2 + cy ** 2) * (ay - by)) / d
        Cy = ((ax ** 2 + ay ** 2) * (cx - bx) +
            (bx ** 2 + by ** 2) * (ax - cx) +
            (cx ** 2 + cy ** 2) * (bx - ax)) / d
        return (Cx, Cy)        

    def draw(self, ax, arr=None):
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        patches = []
        for icell in range(self.ncells):
            vertices = []
            for ivert in self.cellvertices[icell]:
                vertices.append( (self.vertices[ivert,0], self.vertices[ivert,1]) )
            poly = Polygon(vertices, fill=False)
            patches.append(poly)
        p = PatchCollection(patches, match_original=True, 
                            facecolors=None,
                            linestyle='-')
        ax.add_collection(p)

        if arr is None:
            return
            
        import matplotlib.cm as cm
        import matplotlib.colors as colors
        norm = colors.normalize(arr.min(), arr.max())
        colors = []
        for val in arr:
            color = cm.jet(norm(val))
            colors.append(val)
        p.set_array(numpy.array(colors))
        p.set_norm(norm)
        return

    def get_usg_csr_data(self, storediagonal=False):
        '''
        
        '''
        #note that self.nja does not include the diagonal
        nja = self.nja
        if storediagonal: nja = nja + self.ncells
        ia = numpy.empty( (self.ncells + 1), dtype=int)
        ja = numpy.empty( (nja), dtype=int)
        area = numpy.empty( (self.ncells), dtype=float)
        cl1 = numpy.empty( (nja), dtype=numpy.float )
        cl2 = numpy.empty( (nja), dtype=numpy.float )
        fahl = numpy.empty( (nja), dtype=numpy.float )
        fldr = numpy.empty( (nja), dtype=numpy.int )

        #go through each cell and build data
        japos = 0
        for nodeid in range(self.ncells):
            
            area[nodeid] = self.get_cell_area(nodeid)
            ia[nodeid] = japos
            if storediagonal:
                #diagonal position
                ja[japos] = -nodeid
                cl1[japos] = 0.
                cl2[japos] = 0.
                fahl[japos] = 0.
                fldr[japos] = 0
                japos += 1
            
            #go through each connection
            for tonodeid, face in self.connectionlist[nodeid]:
                ja[japos] = tonodeid
                fldr[japos] = 1 #fldir
                cx, cy = self.get_cell_centroid(nodeid)
                xf0 = self.vertices[face[0], 0]
                yf0 = self.vertices[face[0], 1]
                xf1 = self.vertices[face[1], 0]
                yf1 = self.vertices[face[1], 1]
                xm = 0.5 * (xf0 + xf1)
                ym = 0.5 * (yf0 + yf1)
                cl1[japos] = ( (xm - cx) ** 2 + (ym - cy) ** 2) ** 0.5
                #do for tonodeid
                cx, cy = self.get_cell_centroid(tonodeid)
                cl2[japos] = ( (xm - cx) ** 2 + (ym - cy) ** 2) ** 0.5
                fahl[japos] = ( (xf0 - xf1) ** 2 + (yf0 - yf1) ** 2) ** 0.5
                japos += 1
        ia[nodeid + 1] = japos
        return sparse.UsgCsrData(self.ncells, self.nja, area, ia, ja, cl1, cl2, fahl, fldr)
        
*/
