# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import time
import sys
import numpy as np
import json
import cv2
import pickle
from scipy.sparse import lil_matrix, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, bicgstab
from scipy.signal import convolve2d

class Chrono:
    
    indentation = 0
    
    def __init__(self, name):
        self.name = name
        sys.stdout.write( '\n%s %s: ' % ('--'*Chrono.indentation, name))
        sys.stdout.flush()
        Chrono.indentation += 1
    def __enter__(self):
        self.t = time.time()
    def __exit__(self, *args, **kwargs):
        # print('%s: %.2f sec' % (self.name, time.time() - self.t))
        sys.stdout.write( '%.2f sec' % (time.time() - self.t))
        sys.stdout.flush()
        Chrono.indentation -= 1
        
        
        

class Edge:
    # types
    VESSEL = 0
    CAPILLARY = 1   
    
    #properties
    # length = 60
    
    def __init__(self, nodes, edgeType=VESSEL):
        self.radius = 0
        #self.length = self.latticeConstant
        self.index = None
        self.nodes = nodes
        self.edgeType = edgeType
   
    def isTip(self):
        return self.nodes[0].isTip() or self.nodes[1].isTip()
        

    ### Properties:

    @property
    def length(self):
        return Vascularization.latticeConstant

    @property
    def viscosity( self):
        return 4 * 1e-6 * (6*np.exp(-0.17*self.radius) + 3.1 - 2.44*np.exp(-0.06*(2*self.radius)**0.645))
    

class Node:

    ROOT = 0
    INTERNAL = 1    
    TIP = 2
    
    def __init__(self, x,y,nodeType=INTERNAL, pressure=0):
        self.x = x
        self.y = y
        self.index = None
        self.nodeType = nodeType
        self.pressure = pressure
        self.neighbors = []
        self.edges     = []
        
    def __str__(self):
        return '(%d,%d)' % (self.x,self.y)        
        
        
    def _connect(self, node, edge ):
        if( node in self.neighbors):
            print('already connected')
        else:
            self.neighbors.append(node)
            self.edges.append(edge)
            
    def _disconnect(self, node, edge ):
        if( node not in self.neighbors):
            print('not connected')
        else:
            self.neighbors.remove(node)
            self.edges.remove(edge)
            
    @staticmethod
    def connect(node1, node2 ):
        edge = Edge([node1,node2])
        node1._connect(node2, edge)
        node2._connect(node1, edge)
        return edge
        
    @staticmethod
    def disconnect( edge):
        edge.nodes[0]._disconnect(edge.nodes[1], edge)
        edge.nodes[1]._disconnect(edge.nodes[0], edge)
        
    def isTip(self):
        #return len(self.neighbors)==1       
        return self.nodeType != Node.ROOT and np.sum([ edge.edgeType != Edge.CAPILLARY for edge in self.edges]) == 1       
        
    def isNeighbor(self, node):
        return node in self.neighbors
    #def update(self):
        

    ### Properties:
        
    @property
    def volume(self):
        V = 0
        for edge in self.edges:
            Vij = np.pi * edge.radius**2 * edge.length
            V += Vij / 2.
        return V

    @property
    def surface(self):
        A = 0
        for edge in self.edges:
            Aij = 2. * np.pi * edge.radius * edge.length
            A += Aij / 2.
        return A
    
class Region:
    def __init__(self, center, radius, **properties):
        self.radius = radius
        self.x = center[0]
        self.y = center[1]
        self.properties = properties
        
    def containes(self, pos):
        # if( isinstance( node_or_edge, Node)):
        #     node = node_or_edge
        #     return (node.x-self.x)**2 + (node.y-self.y)**2 < self.radius**2
        # if( isinstance( node_or_edge, Edge)):
        #     edge = node_or_edge
        #     node1 = edge.nodes[0]
        #     node2 = edge.nodes[1]
        #     return (node1.x-self.x)**2 + (node1.y-self.y)**2 < self.radius**2 or (node2.x-self.x)**2 + (node2.y-self.y)**2 < self.radius**2
        # if( isinstance( pos, [tuple])):
        x = pos[0]
        y = pos[1]
        return (x-self.x)**2 + (y-self.y)**2 < self.radius**2
            
    def getProperty(self, pos, name):
        # print( name, self.properties)
        if( name in self.properties and self.containes( pos)):
            return self.properties[ name]
        else:
            return None
        
        
        
        
class Vascularization:
    
    latticeConstant = 60.
    
    def __init__(self, width, height, latticeConstant=60.):
        ''' Returns the reversed String.
         Parameters:
            width (int): The string which is to be reversed.
            height (int): The string which is to be reversed.
            latticeConstant (float): The string which is to be reversed.
         Returns:
            reverse(str1): The string which gets reversed.  
        '''
        # properties
        self.width = width
        self.height = height

        Vascularization.latticeConstant = latticeConstant
        
        self.nodes = []
        self.edges = []
        self.regions = []
        self.grid = np.empty((width, height), dtype=object)
    
    def __str__(self):
        return '(nodes: %d, edges: %d)' %(len(self.nodes), len(self.edges))
    
    # def save(self, filename):
    #     print('build json')
        
    def save(self, filename):
        pickle.dump( self, open( filename, 'wb'))
    
    @staticmethod
    def load( filename):
        vascularization = pickle.load( open( filename, 'rb'))
        return vascularization
    
    def addNode(self, node):
        node.index = len(self.nodes)
        self.nodes.append(node)
        self.grid[node.x][node.y] = node
    
    def removeNode(self, node):  
        self.grid[node.x][node.y] = None
        i = node.index
        j = len(self.nodes) - 1
        self.nodes[i], self.nodes[j] = self.nodes[j], self.nodes[i]
        self.nodes[i].index = i
        self.nodes.pop()
    
    def addEdge(self, edge):
        edge.index = len(self.edges)
        self.edges.append(edge)
    
    def addRegion(self, region):
        self.regions.insert( 0, region)
    
    def _posFree( self,x,y):
        return self.grid[x][y] is None
    def _posInDomain( self,x,y):
        return x>=0 and y>=0 and x<self.width and y<self.height
    
    def _getProperty(self, pos, name, default):
        # print( name, pos)
        for r in self.regions:
            value = r.getProperty( pos, name)
            # print(value)
            if( value is not None):
                # print('got region value: %s=%f' % (name, value))
                return value
        return default
        
    
    def addCapillaries(self):
        for i, node in enumerate(self.nodes):
            if( node.isTip()):
                #print('found tip')
                DX = [[ 0,-1],
                      [ 0,+1],
                      [-1, 0],
                      [+1, 0]]
                for dx in DX: 
                    x,y = node.x+dx[0],node.y+dx[1]
                    if( self._posInDomain( x,y) and not self._posFree( x,y)):
                        #print('found neighbor')
                        neighbor = self.grid[x,y]
                        if( neighbor.isTip() and not node.isNeighbor(neighbor)):
                            new_edge = Node.connect( node,neighbor)
                            new_edge.edgeType = Edge.CAPILLARY
                            new_edge.radius   = 4
                            self.addEdge( new_edge)
    
    def removeCapillaries(self):
        for edge in reversed(self.edges):
            if(edge.edgeType == Edge.CAPILLARY):
                #print('Remove capillary')
                self.edges.remove( edge)
                Node.disconnect( edge)
    
    def grow(self, iterations=1, sproutingProbability=1.0):
        for it in range(iterations):
            n = len(self.nodes)
            I = np.random.permutation(n)
            #for node in self.nodes[:n]:
            
            count_sprouts = 0
            
            for i in I:
                node = self.nodes[i]
                #print( str(node))
                #dx = np.random.randint(2)*2-1
                DX = [[0,-1],
                      [0,+1],
                      [-1,0],
                      [+1,0]]
                
                baseProbability = self._getProperty( (node.x,node.y), 'sproutingProbability', sproutingProbability)
                
                for dx in DX:    
                    if( self._posInDomain( node.x+dx[0],node.y+dx[1]) and self._posFree( node.x+dx[0],node.y+dx[1]) and np.random.rand() < baseProbability/4.):
                        new_node = Node( node.x+dx[0],node.y+dx[1])
                        new_edge = Node.connect( node,new_node)
                        self.addNode( new_node)
                        self.addEdge( new_edge)
                        count_sprouts += 1
                        #no
                    #print( )
            # print('%d sprouts in iteration %d' % (count_sprouts, it))

                
   # def grow(self, sproutingProbability=0.5):
        

    def toEPS(self, filename, plotRoots=False, plotTips=False):
        scale = .1
        with open( filename, 'w+') as fp:
            fp.write('%!PS-Adobe-3.1 EPSF-3.0\n')
            fp.write('%%%%BoundingBox: 0 0 %d %d\n' % (self.width*self.latticeConstant*scale,self.height*self.latticeConstant*scale))
            for node in self.nodes:
                if( plotRoots and node.nodeType == Node.ROOT ):
                    fp.write('%d %d %f 0 3self.latticeConstant arc\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale, self.latticeConstant*scale))
                if( plotTips and node.isTip()):
                    fp.write('%d %d %f 0 3self.latticeConstant arc\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale, self.latticeConstant/3*scale))
                    
                for i, neighbor in enumerate(node.neighbors):
                    pressure = (node.pressure + neighbor.pressure) / 2.
                    rel = (pressure - 2.) / (12. - 2.)
                    fp.write('%f %f %f setrgbcolor\n' % (1-rel,0,rel))
                    fp.write('%f setlinewidth\n' % (node.edges[i].radius*scale))
                    fp.write('%d %d moveto\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale))
                    fp.write('%d %d lineto\n' % (neighbor.x*self.latticeConstant*scale,neighbor.y*self.latticeConstant*scale))
                    fp.write('stroke\n')    
                    
    def toImage(self, filename, plotRoots=False, plotTips=False, plotRegions=True):
        scale = int(.1 * self.latticeConstant)
        
        img_h = int( self.height*scale)
        img_w = int( self.width*scale)
        img = np.ones( (img_h, img_w, 3))*255
        
        print( 'Image: '+str(img.shape))
        
        if( plotRegions):
            img_r = img.copy()
            colors = [[255,0,0], [0,255,0], [0,0,255]]
            for i, r in enumerate(self.regions):
                cv2.circle( img_r, (r.x*scale,r.y*scale), r.radius*scale, colors[i], -1)
            alpha = 0.3
            img = cv2.addWeighted( img_r, alpha, img, 1-alpha, 0)


        for node in self.nodes:
            # if( plotRoots and node.nodeType == Node.ROOT ):
            #     fp.write('%d %d %f 0 3self.latticeConstant arc\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale, self.latticeConstant.*scale))
            # if( plotTips and node.isTip()):
            #     fp.write('%d %d %f 0 3self.latticeConstant arc\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale, 20.*scale))
                
            for i, neighbor in enumerate(node.neighbors):
                pressure = (node.pressure + neighbor.pressure) / 2.
                rel = (pressure - 2.) / (12. - 2.)
                color = (256-int(256*rel), 0, int(256*rel))
                pt1 = (int(node.x*scale), int(node.y*scale))
                pt2 = (int(neighbor.x*scale), int(neighbor.y*scale))
                thickness = int( np.ceil(node.edges[i].radius/self.latticeConstant*scale))
                cv2.line( img, pt1,pt2, color, thickness)
                # fp.write('%f %f %f setrgbcolor\n' % (1-rel,0,rel))
                # fp.write('%f setlinewidth\n' % (node.edges[i].radius*scale))
                # fp.write('%d %d moveto\n' % (node.x*self.latticeConstant*scale,node.y*self.latticeConstant*scale))
                # fp.write('%d %d lineto\n' % (neighbor.x*self.latticeConstant*scale,neighbor.y*self.latticeConstant*scale))
                # fp.write('stroke\n')
                
                
        cv2.imwrite( filename, img)
    
                    
    def toEPSshear(self, filename):
        
        shearStress = [ edge.shearStress for edge in self.edges ]     
        shearStress_min = np.min(shearStress)
        shearStress_max = np.max(shearStress)
        
        scale = .1
        with open( filename, 'w+') as fp:
            fp.write('%% BoundingBox: 0 0 %d %d\n' % (self.width*self.latticeConstant*scale,self.height*self.latticeConstant*scale))
            for edge in self.edges:
                rel = (edge.shearStress - shearStress_min) / (shearStress_max - shearStress_min)
                fp.write('%f %f %f setrgbcolor\n' % (rel,0,1-rel))
                fp.write('%f setlinewidth\n' % (edge.radius*scale))
                fp.write('%d %d moveto\n' % (edge.nodes[0].x*self.latticeConstant*scale,edge.nodes[0].y*self.latticeConstant*scale))
                fp.write('%d %d lineto\n' % (edge.nodes[1].x*self.latticeConstant*scale,edge.nodes[1].y*self.latticeConstant*scale))
                fp.write('stroke\n')
                    
    def calculateRadii(self):
        # print('init radii')plotRoots
        #todo = []
        count_radii_set= 0        
        
        for edge in self.edges:
            if( edge.isTip()):
                edge.radius = 1
                count_radii_set += 1
            else:
                edge.radius = 0
                #todo.append( edge)
                
        #todo = list( self.nodes)
        #n=len(todo)
        
        
        for it in range(1000):
            #print(n)
            count_not_set_nodes = 0
            for node in self.nodes:#self.nodes:
                count_not_set = 0
                edge_not_set = None
                for (edge,neighbor) in zip(node.edges,node.neighbors):
                    if(edge.radius == 0):
                        count_not_set += 1
                        edge_not_set = edge
                        node_not_set = neighbor
                        
                already_set = False
                
                if( count_not_set==0):
                    already_set = True
                
                if( count_not_set==1 and node.nodeType != Node.ROOT ):
                    edge_not_set.radius = min( self.latticeConstant, np.sum([ edge.radius**2 for edge in node.edges ])**0.5)
                    count_radii_set += 1
                    already_set = True
                    
                if( not already_set):
                    count_not_set_nodes += 1
                    
            #print(count_not_set_nodes)
            if( count_not_set_nodes == 0):
                # print( 'set %d / %d radii ' % (count_radii_set, len(self.edges)))
                return
        #print('')
                
    def calculatePressure(self):
        n = len(self.nodes)
        # print('create sparse matrix: %dx%d' % (n,n))
        # A = csc_matrix((n,n))
        A = lil_matrix((n,n))
        b = [0] * n
        
        
        with Chrono('build A'):
            for i, node in enumerate(self.nodes):
                assert( i == node.index)
                if(node.nodeType == Node.ROOT):
                    A[i,i] = 1
                    b[i] = node.pressure
                else:
                    A[i,i] = 0
                    #b[i] = 0
                    for neighbor, edge in zip(node.neighbors, node.edges):
                        j = neighbor.index
                        Gij = np.pi / 8. * edge.radius**4 / (edge.length * edge.viscosity)
                        #print('i:%d, j:%d' % (i,j))
                        A[i,i] += Gij
                        #print(j)
                        A[i,j]  = -Gij

        with Chrono('solve x=b/A'):
            x = spsolve(A,b)
        #x = bicgstab(A,b)[0]
        # print(x)              
        with Chrono('update pressure'):
            for i, node in enumerate(self.nodes):  
                node.pressure = x[i] 
    
    def calculateShearStress(self):   
        for edge in self.edges:
            edge.shearStress = np.fabs( edge.nodes[0].pressure - edge.nodes[1].pressure) * edge.radius / (2 * edge.length)
            
    def collapse(self, iterations=1, collapseProbability=1.0):

        shearStress = [ edge.shearStress for edge in self.edges ]     
        shearStress_min = np.min(shearStress)
        shearStress_max = np.max(shearStress)
        # print( 'shear stress: min=%f, max=%f, mean=%f' % (shearStress_min, shearStress_max, np.mean(shearStress)))
        
        for it in range(iterations):
            count_collapses = 0        
            
            #np.random.shuffle( self.edges)
            #tips = [ edge for edge in reversed(self.edges) if edge.isTip()]
            for edge in reversed(self.edges):
            #for edge in tips:
        
                baseProbability = self._getProperty( (edge.nodes[0].x,edge.nodes[0].y), 'collapseProbability', collapseProbability)
                proba = baseProbability * (shearStress_max - edge.shearStress) / (shearStress_max-shearStress_min)            
                # proba = collapseProbability * (shearStress_max - edge.shearStress) / (shearStress_max-shearStress_min)            
                
                
                if( edge.isTip() and np.random.rand() < proba):
                    #print('collapse')
                    self.edges.remove(edge)
                    Node.disconnect( edge)
                    if( len(edge.nodes[0].neighbors) == 0):
                        self.removeNode(edge.nodes[0])
                    if( len(edge.nodes[1].neighbors) == 0):
                        self.removeNode(edge.nodes[1])
                    count_collapses += 1
            # print('%d collapses' % count_collapses)
            
    def calculatePerfusion(self, starttime, endtime, timestep, AIF, Kps=1, diffusionCoef=0):
        timestep_plot = 0.1 # sec
        
        n = len(self.nodes)
        # print('create sparse matrix: %dx%d' % (n,n))
        # A = csc_matrix((n,n))
        # A = lil_matrix((n,n))

        # precalculate constant properties
        v = [ node.volume for node in self.nodes]
        s = [ node.surface for node in self.nodes]
        
        b = np.zeros( (n))
        x = np.zeros( (n))
        concP = np.zeros((self.width,self.height))
        
        concI0 = np.zeros((self.width,self.height))
        concI1 = np.zeros((self.width,self.height))

        conc = np.zeros((self.width,self.height))
        
        # t = 0
        for t in np.arange(starttime, endtime, timestep):
            with Chrono('explicit Euler (t=%f sec)' % (t)):
                
                concI1[:,:] = concI0[:,:]
                
                for i, node in enumerate(self.nodes):
                    assert( i == node.index)
                    if(node.nodeType == Node.ROOT):
                        # A[i,i] = 1
                        b[i] = AIF(t)
                    else:
                        # A[i,i] = 1
                        b[i] = x[i] # current conecntration
                        for neighbor, edge in zip(node.neighbors, node.edges):
                            j = neighbor.index
                            Gij = np.pi / 8. * edge.radius**4 / (edge.length * edge.viscosity)
                            Fij = Gij * (neighbor.pressure - node.pressure) * timestep
                            if( Fij > 0):
                                # influx
                                b[i] += x[j] * Fij / v[i]
                            else:
                                # outflux
                                b[i] += x[i] * Fij / v[i]
                    perm = self._getProperty( [node.x, node.y], 'Kps', Kps)
                    
                    maxVolFrac = 0.9
                    if( self.latticeConstant**3*maxVolFrac < v[i]):
                        # print('vessel volume too large')
                        concI1[node.x, node.y] += perm*s[i] * (x[i] - concI0[node.x, node.y]) / ((1.-maxVolFrac)*self.latticeConstant**3) * timestep
                    else:
                        concI1[node.x, node.y] += perm*s[i] * (x[i] - concI0[node.x, node.y]) / (self.latticeConstant**3 - v[i]) * timestep
                                
            # Diffusion
            # DIFFUSION_COEFFICIENT = 0.001
            if( (b < 0).any() ):
                print('negative conc (plasma)')
            if( (concI1 < 0).any() ):
                print('negative conc (interstitial)')
            concI1 += diffusionCoef * convolve2d( concI0, [[0, 1, 0], [1, -4, 1], [0, 1, 0]], mode='same', boundary='wrap') * timestep
    
            # with Chrono('solve x=b/A'):
                # x = spsolve(A,b)
                
            x[:] = b[:]
            concI0[:,:] = concI1[:,:]
            
            # plot
            if( int(t/timestep_plot) != int((t+timestep)/timestep_plot) ):
                conc[:,:] = concI0[:,:]
                for i, node in enumerate(self.nodes):
                    concP[ node.x, node.y] = x[i] * v[i] / (self.latticeConstant**3)
                    
                    conc[ node.x, node.y] = concP[ node.x, node.y] + concI0[ node.x, node.y]*(1 - v[i]/self.latticeConstant**3)

                # print( )
                it = int(t/timestep_plot)

                cv2.imwrite( 'concP_%03d.png' % (it), concP.T/1.*256)
                cv2.imwrite( 'concI_%03d.png' % (it), concI0.T/1*256)
                cv2.imwrite( 'conc_%03d.png' % (it), conc.T/1*256)
                    
                # print( np.mean(x))
                
            t += timestep
        #x = bicgstab(A,b)[0]
        # print(x)              
        # with Chrono('update pressure'):
        #     for i, node in enumerate(self.nodes):  
        #         node.pressure = x[i] 
        
        
def AIF(t):
    t = t / 60 # sec to min
    
    A = [0.809, 0.330]
    T = [0.171, 0.365]
    S = [0.056, 0.132]
    a = 1.050
    b = 0.169
    s = 38.078
    tau = 0.483
    
    return np.sum([ np.exp(-(t-T[i])**2/(2*S[i])**2) + a*np.exp(-b*t) / (1 + np.exp(s*(tau-t))) for i in [0,1]]) 
    
# if( __name__=='__main__'):
    
#     import pickle
#     vascularization = pickle.load( open('vask.pkl', 'rb'))
#     # vascularization = Vascularization( 100, 100)
    
#     # vascularization.addNode( Node( 50,50, Node.ROOT, 12))

#     # vascularization.addNode( Node(  0,50, Node.ROOT, 2))
#     # vascularization.addNode( Node( 80, 0, Node.ROOT, 2))
#     # vascularization.addNode( Node( 80,99, Node.ROOT, 2))

#     # vascularization.grow( 70)
     
    
#     with Chrono('grow'):
#         vascularization.grow(10, sproutingProbability=0.25)
#     with Chrono('calculateRadii'):
#         vascularization.calculateRadii()
#     with Chrono('addCapillaries'):
#         vascularization.addCapillaries() 
#     with Chrono('calculatePressure'):
#         vascularization.calculatePressure()
            
#     vascularization.toEPS('before_perfusion.eps')
#     vascularization.calculatePerfusion( 0.1, AIF)
    
# else:
#     for i in range(200):
#         with Chrono('grow'):
#             vascularization.grow(10, sproutingProbability=0.25)
#         with Chrono('calculateRadii'):
#             vascularization.calculateRadii()

#         # print( 'before: '+str(vascularization))
#         #vascularization.toEPS('bla_before.eps')
#         with Chrono('addCapillaries'):
#             vascularization.addCapillaries() 
#         with Chrono('calculatePressure'):
#             vascularization.calculatePressure()
#         with Chrono('calculateShearStress'):
#             vascularization.calculateShearStress()
#         vascularization.toEPS('bla_between.eps')
#         vascularization.toEPS('bla_%i.eps' % (i))
#         #vascularization.toEPSshear('shear_between.eps')
#         # print( 'between: '+str(vascularization))
#         with Chrono('removeCapillaries'):
#             vascularization.removeCapillaries()
#         with Chrono('collapse'):
#             vascularization.collapse( 1, collapseProbability=1)
#         # print( 'after: '+str(vascularization))
#         #vascularization.calculatePressure()        
#         #vascularization.toEPS('bla_after.eps')
        
        
#         pickle.dump(vascularization, open('vask.pkl', 'wb'))
    