# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import time
import sys
import numpy as np
import json
from scipy.sparse import lil_matrix, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, bicgstab

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
    
    VESSEL = 0
    CAPILLARY = 1    
    
    def __init__(self, nodes, edgeType=VESSEL):
        self.radius = 0
        self.length = 60
        self.index = None
        self.nodes = nodes
        self.edgeType = edgeType
   
    def isTip(self):
        return self.nodes[0].isTip() or self.nodes[1].isTip()
        
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
        
        
    def connect_(self, node, edge ):
        if( node in self.neighbors):
            print('already connected')
        else:
            self.neighbors.append(node)
            self.edges.append(edge)
            
    def disconnect_(self, node, edge ):
        if( node not in self.neighbors):
            print('not connected')
        else:
            self.neighbors.remove(node)
            self.edges.remove(edge)
            
    @staticmethod
    def connect(node1, node2 ):
        edge = Edge([node1,node2])
        node1.connect_(node2, edge)
        node2.connect_(node1, edge)
        return edge
        
    @staticmethod
    def disconnect( edge):
        edge.nodes[0].disconnect_(edge.nodes[1], edge)
        edge.nodes[1].disconnect_(edge.nodes[0], edge)
        
    def isTip(self):
        #return len(self.neighbors)==1       
        return self.nodeType != Node.ROOT and np.sum([ edge.edgeType != Edge.CAPILLARY for edge in self.edges]) == 1       
        
    def isNeighbor(self, node):
        return node in self.neighbors
    #def update(self):
        
    def volume(self):
        V = 0
        for edge in self.edges:
            Vij = np.pi * edge.radius**2 * edge.length
            V += Vij / 2.
        return V
    def surface(self):
        A = 0
        for edge in self.edges:
            Aij = 2. * np.pi * edge.radius * edge.length
            A += Aij / 2.
        return A
    

class Vascularization:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.nodes = []
        self.edges = []
        self.grid = np.empty((width, height), dtype=object)
    
    def __str__(self):
        return '(nodes: %d, edges: %d)' %(len(self.nodes), len(self.edges))
    
    # def save(self, filename):
    #     print('build json')
        
    
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
    
    
    def posFree( self,x,y):
        return self.grid[x][y] is None
    def posInDomain( self,x,y):
        return x>=0 and y>=0 and x<self.width and y<self.height
    
    
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
                    if( self.posInDomain( x,y) and not self.posFree( x,y)):
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
    
    def grow(self, iterations=1, sproutingProbability=0.5):
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
                for dx in DX:    
                    if( self.posInDomain( node.x+dx[0],node.y+dx[1]) and self.posFree( node.x+dx[0],node.y+dx[1]) and np.random.rand() < sproutingProbability):
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
            fp.write('%%%%BoundingBox: 0 0 %d %d\n' % (self.width*60*scale,self.height*60*scale))
            for node in self.nodes:
                if( plotRoots and node.nodeType == Node.ROOT ):
                    fp.write('%d %d %f 0 360 arc\n' % (node.x*60*scale,node.y*60*scale, 60.*scale))
                if( plotTips and node.isTip()):
                    fp.write('%d %d %f 0 360 arc\n' % (node.x*60*scale,node.y*60*scale, 20.*scale))
                    
                for i, neighbor in enumerate(node.neighbors):
                    pressure = (node.pressure + neighbor.pressure) / 2.
                    rel = (pressure - 2.) / (12. - 2.)
                    fp.write('%f %f %f setrgbcolor\n' % (1-rel,0,rel))
                    fp.write('%f setlinewidth\n' % (node.edges[i].radius*scale))
                    fp.write('%d %d moveto\n' % (node.x*60*scale,node.y*60*scale))
                    fp.write('%d %d lineto\n' % (neighbor.x*60*scale,neighbor.y*60*scale))
                    fp.write('stroke\n')
                    
    def toEPSshear(self, filename):
        
        shearStress = [ edge.shearStress for edge in self.edges ]     
        shearStress_min = np.min(shearStress)
        shearStress_max = np.max(shearStress)
        
        scale = .1
        with open( filename, 'w+') as fp:
            fp.write('%% BoundingBox: 0 0 %d %d\n' % (self.width*60*scale,self.height*60*scale))
            for edge in self.edges:
                rel = (edge.shearStress - shearStress_min) / (shearStress_max - shearStress_min)
                fp.write('%f %f %f setrgbcolor\n' % (rel,0,1-rel))
                fp.write('%f setlinewidth\n' % (edge.radius*scale))
                fp.write('%d %d moveto\n' % (edge.nodes[0].x*60*scale,edge.nodes[0].y*60*scale))
                fp.write('%d %d lineto\n' % (edge.nodes[1].x*60*scale,edge.nodes[1].y*60*scale))
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
                    edge_not_set.radius = min( 60, np.sum([ edge.radius**2 for edge in node.edges ])**0.5)
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
                        Gij = np.pi / 8. * edge.radius**4 / (edge.length * edge.viscosity())
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
            
    def collapse(self, iterations=1, collapseProbability=0.2):

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
        
                proba = collapseProbability * (shearStress_max - edge.shearStress) / (shearStress_max-shearStress_min)            
                
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
            
    def calculatePerfusion(self, timestep, AIF, Kps=1):
        n = len(self.nodes)
        # print('create sparse matrix: %dx%d' % (n,n))
        # A = csc_matrix((n,n))
        # A = lil_matrix((n,n))
        v = [ node.volume() for node in self.nodes]
        s = [ node.surface() for node in self.nodes]
        b = [0] * n
        x = [0] * n
        concP = np.zeros((100,100))
        
        concI0 = np.zeros((100,100))
        concI1 = np.zeros((100,100))
        
        t = 0
        for it in range(10000):
            with Chrono('build A'):
                
                for i in range(100):
                    for j in range(100):
                        concI1[i][j] = concI0[i][j]
                
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
                            Gij = np.pi / 8. * edge.radius**4 / (edge.length * edge.viscosity())
                            # print(Gij)
                            # Gij = 1
                            Fij = Gij * (neighbor.pressure - node.pressure) * timestep
                            # print(Fij)
                            if( Fij > 0):
                                # influx
                                b[i] += x[j] * Fij / v[i]
                            else:
                                # outflux
                                b[i] += x[i] * Fij / v[i]
                    concI1[node.x, node.y] += Kps*s[i] * (x[i] - concI0[node.x, node.y]) / (60**3 - v[i]) * timestep
                                
    
            # with Chrono('solve x=b/A'):
                # x = spsolve(A,b)
                
            x[:] = b[:]
            concI0[:,:] = concI1[:,:]
            if( it%10==0):
                for i, node in enumerate(self.nodes):
                    concP[ node.x, node.y] = x[i] * v[i] / (60**3)
                # import matplotlib.pyplot as plt
                import cv2
                # plt.imshow( conc)
                # plt.draw()
                # plt.savefig('conc.png')
                cv2.imwrite( 'concP_%03d.png' % (it/10), concP/1*256)
                cv2.imwrite( 'concI_%03d.png' % (it/10), concI0/6*256)
                cv2.imwrite( 'conc_%03d.png' % (it/10), (concI0 + concP)/6.*256)
                    
                print( np.mean(x))
                
            t += timestep
        #x = bicgstab(A,b)[0]
        # print(x)              
        # with Chrono('update pressure'):
        #     for i, node in enumerate(self.nodes):  
        #         node.pressure = x[i] 
        
        
def AIF(t):
    t = t / 60. # sec to min
    
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
    