#import math
#import PDBTools
#import Vector
class Node:
    def __init__(self,loc,leftChild,rightChild,axis):
        self.location = loc
        self.leftChild = leftChild
        self.rightChild = rightChild
        self.leaf = False
        self.axis = axis
class Leaf:
    def __init__(self,atom):
        self.location = atom.location
        self.data = atom
        self.leaf = True
        
class KDTree:
    DIMENSIONS = 3
    def __init__(me):
        me.root = None
    @staticmethod
    def loadAtomArray(ary):
        me = KDTree()
        if (len(ary) == 0):
            me.root = None
        else:
            me.root = makeTree(ary)
        me.myList = ary
        return me
    @staticmethod
    def load2DPoints(ary):
        me = KDTree()
        me.root = makeTree(ary)
        return me
    def radiusSearch(me,center,radius,returnWithDist=False,condition = None):
        def searchNode(node):
            if node is None:
                return
            if node.leaf:
                dist = node.location.distance(center)
                if (dist < radius):
                    if (condition is None) or (condition(node.data)):
                        if returnWithDist:
                            masterList.append((node.data,dist))
                        else:
                            masterList.append(node.data)
                return
            else:
                axis = node.axis       
                median = node.location.dims[axis]
                low = center.dims[axis] - radius
                high = center.dims[axis] + radius
                if (low <= median):
                    searchNode(node.leftChild)
                if (high >= median):
                    searchNode(node.rightChild)
                return
        #radiusSQ = radius ** 2
        if (me.root is None):
            return []
        masterList = []
        searchNode(me.root)
        return masterList
    def nearestNeighbor(me,center,returnWithDist=False,condition=None):
        result = me.nearestNeighborWithDist(center,condition)
        if returnWithDist:
            return result
        return result[0]
    def nearestNeighborWithDist(me,center,condition=None):
        if (me.root is None):
            return (None,None)
        rootDist = center.distance(me.root.location)
        tup = NNS(center,me.root,None,rootDist,condition)
        return tup[0].data,tup[1]
def NNS(point,node,bestPoint,bestDist,condition=None):
    if node is None:
        #print("none")
        return
    #callID = KDTree.depth
    callID = 5
    #print("\nCall ",callID)
    #KDTree.depth+=1
    if (node.leaf):
        if (condition is None) or (condition(node.data)):
            newDist = point.distance(node.location)
            if (bestPoint is None) or (newDist <= bestDist):
                #print("Improve ",node.leaf)
                bestDist = newDist
                bestPoint = node
        #print("\tLeaf  ",bestDist)
    else:
        #print("not a leaf")
        axis = node.axis
        value = node.location.dims[axis]
        #print(axis,"  ",value)
        if (point.dims[axis] <= node.location.dims[axis]):
            searchLeftFirst = True
        else:
            searchLeftFirst = False
        #print("Left ",searchLeftFirst)
        #print(bestPoint.location.dims[axis], bestDist, value,bestPoint.location.dims[axis] - bestDist, bestPoint.location.dims[axis] + bestDist )
        if (searchLeftFirst):
            #print("Left top method")
            if (point.dims[axis] - bestDist <= value):
                #print("Search left ",callID)
                bestPoint,bestDist = NNS(point,node.leftChild,bestPoint,bestDist,condition)
            if (point.dims[axis] + bestDist >= value):
                #print("Search right",callID)
                bestPoint,bestDist = NNS(point,node.rightChild,bestPoint,bestDist,condition)
        else:
            if (point.dims[axis] + bestDist >= value):
                #print("Search right",callID)
                bestPoint,bestDist = NNS(point,node.rightChild,bestPoint,bestDist,condition)            
            if (point.dims[axis] - bestDist <= value):
                #print("Search left",callID)
                bestPoint,bestDist = NNS(point,node.leftChild,bestPoint,bestDist,condition)
    return bestPoint,bestDist
            
def makeTree(ary,depth=0):
    if (len(ary) == 1):
        return Leaf(ary[0])
    axis = depth % KDTree.DIMENSIONS
 
    # Sort point list and choose median as pivot element
    ary  = sorted(ary, key=lambda x:x.location.dims[axis])
    median = len(ary) // 2 # choose median
 
    # Create node and construct subtrees
    return Node(ary[median].location,makeTree(ary[:median], depth + 1),makeTree(ary[median:], depth + 1),axis)

