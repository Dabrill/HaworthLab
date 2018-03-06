import Permutations
CUTOFF = 200
class Node:
    def __init__(self,value,leftChild,rightChild,depth):
        self.value = value
        self.leftChild = leftChild
        self.rightChild = rightChild
        self.leaf = False
        self.depth = depth
class Leaf:
    def __init__(self,poseNum,depth):
        self.data = poseNum
        self.leaf = True
        self.depth = depth
class Tree:
    def __init__(me,data):
        me.root = Tree.makeTree(data,list(range(data.size)),0)
        me.data = data
    @staticmethod
    def makeTree(data,ary,depth):
        if (len(ary) == 1):
            return Leaf(ary[0],depth)
        atom = depth//3
        dim = depth%3
        ary = sorted(ary,key = lambda x:data.getCoord(x,atom,dim))
        medianIndex = len(ary) // 2
        medianPoseNumber = ary[medianIndex]
        medianVal = data.getCoord(medianPoseNumber,atom,dim)
        return Node(medianVal,Tree.makeTree(data,ary[:medianIndex], depth + 1),Tree.makeTree(data,ary[medianIndex:], depth + 1),depth)
    def radiusSearch(me,poseNum,radius,returnWithDist=False,condition = None):
        def searchNode(node):
            if node is None:
                return
            if node.leaf:
                if (node.data != poseNum):
                    if me.remainderMatches(poseNum,node.data,node.depth):
                        masterList.append(node.data)
                return
            else:
                atom = node.depth//3
                axis = node.depth%3
                searchValue = me.data.getCoord(poseNum,atom,axis)
                low = searchValue - radius
                high = searchValue + radius
                if (low <= node.value):
                    searchNode(node.leftChild)
                if (high >= node.value):
                    searchNode(node.rightChild)
                return
        if (me.root is None):
            return []
        masterList = []
        searchNode(me.root)
        return masterList
    def remainderMatches(me,aNum,bNum,depth):
        #print(aNum," to ",bNum," at depth ",depth)
        #for i in range(depth,me.data.atomNum*3):
        for i in range(me.data.atomNum*3):
            atom = i//3
            axis = i%3
            A = me.data.getCoord(aNum,atom,axis)
            B = me.data.getCoord(bNum,atom,axis)
            if (abs(A-B) > CUTOFF):
                #print("\t",atom,axis,A,B)
                return False
        return True

def findOverlays(data):
    def EValue(x):
        if data.fixvar.poses[x-1].energy is None:
            return 0
        return data.fixvar.poses[x-1].energy
        
    #outf = open("OverlayResultsCull.out",'w')
    #data.size = 10000
    dataTree = Tree(data)
    print("Tree constructed")

    count = 0
    for poseNumA in range(1,data.size+1):
        if (data.poseGood[poseNumA-1]):
            if (count % 100 == 0):
                print(poseNumA)
            #print(poseNumA)
            group = dataTree.radiusSearch(poseNumA,CUTOFF)
            group.append(poseNumA)
            group = sorted(group,key = EValue)
            '''
            for m in group:
                print("\t",m,data.fixvar.poses[m-1].energy)
            '''
            for j in range(1,len(group)):
                data.bad(group[j])
            count+=1
    #outf.close()
