class Node:
    def __init__(me,ID):
        me.ID = ID
        me.edges = []
        me.rEdges = []
    def toString(me):
        s = "{"+me.ID+"}\n->\n"
        for e in me.edges:
            s+=e.toString()+"\n"
        s+="<-\n"
        for e in me.rEdges:
            s+=e.toString()+"\n"
        return s
    def __lt__(self,other):
        if (self.__cmp__(other) < 0):
            return True
        return False
    def __gt__(self,other):
        if (self.__cmp__(other) > 0):
            return True
        return False
    def __cmp__(self,other):
        if (self.cost is None):
            if (other.cost is None): return 0
            return 1
        else:
            if (other.cost is None): return -1
        if (self.cost < other.cost):    return -1
        if (self.cost == other.cost):   return 0
        if (self.cost > other.cost):    return 1
        return None
class Edge:
    def __init__(me,source,sink,cost,num):
        me.source = source
        me.sink = sink
        me.cost = cost
        me.num = num
        me.include = False
    def invert(me):
        me.sink.edges.append(me)
        me.sink.rEdges.remove(me)
        me.source.edges.remove(me)
        me.source.rEdges.append(me)
        temp = me.sink
        me.sink = me.source
        me.source = temp
        me.cost = 0 - me.cost
        me.include = not me.include
        return me
    def wholeCost(me):
        return me.cost  - me.sink.price + me.source.price
    def toString(me):
        return str(me.source.ID)+"=["+str(me.cost)+"]=>"+str(me.sink.ID)
class Heap:
    def __init__(me,n):
        me.array = [None] * n
        me.dict = {}
        me.size = 0
    def openPaths(me):
        if (me.size < 1):
            return False
        return me.array[0].cost is not None
    def parent(me,i):
        return int((i-1)/2)
    def insert(me,item):
        me.dict[item.ID] = me.size
        me.array[me.size] = item
        me.heapifyUp(me.size)
        me.size+=1
    def delete(me,i):
        val = me.array[i].ID
        me.array[i] = me.array[me.size-1]
        me.size-=1
        me.heapifyDown(i)   
        me.dict[val] = None
    def changeKey(me,ID,val):
        i = me.dict[ID]
        if (i == None):
            return None
        me.array[i].cost = val
        me.heapifyUp(i)
        me.heapifyDown(i)
        return i
    def swap(me,i,j):
        me.dict[me.array[i].ID] = j
        temp = me.array[i]
        me.array[i] = me.array[j]
        me.array[j] = temp
        me.dict[me.array[i].ID] = i
    def heapifyUp(me,i):
        if (i > 0):
            j = me.parent(i)
            if (me.array[i] < me.array[j]):
                me.swap(i,j)
                me.heapifyUp(j)
    def heapifyDown(me,i):
        top = 2 * i + 2
        if (top > me.size):
            return
        if (top < me.size):
            left = 2 * i + 1
            right = 2* i + 2
            if (me.array[left] < me.array[right]):
                j = left
            else:
                j = right
        elif (top == me.size):
            j = top - 1
        if (me.array[i] > me.array[j]):
            me.swap(i,j)
            me.heapifyDown(j)
    def getMin(me):
        return me.array[0]
    def extractMin(me):
        ret = me.array[0]
        me.delete(0)
        return ret
class Network:
    def __init__(me):
        me.nodes = {}
        me.edgeList = []
        me.cutoff = 3.0
    def addNode(me,ID):
        if not ID in me.nodes:
            me.nodes[ID] = Node(ID)
    def addNodes(me,IDs):
        for ID in IDs:
            me.addNode(ID)
    def addEdgesTo(me, u, vs, c):
        for v in vs:
            me.addEdge(u,v,c)
    def addEdgesFrom(me, us, v, c):
        for u in us:
            me.addEdge(u,v,c)
    def addEdgeLoad(me, u, v, c):
        if not u in me.nodes:
            me.addNode(u)
        if not v in me.nodes:
            me.addNode(v)
        if u == v:
            raise ValueError("can't add edge to self")
        edge = Edge(me.nodes[u],me.nodes[v],c,len(me.edgeList))
        me.nodes[u].edges.append(edge)
        me.nodes[v].rEdges.append(edge)
        me.edgeList.append(edge)
    def addMutualEdge(me,u,v,c):
        me.addEdge(u,v,c)
        me.addEdge(v,u,c)
    def addEdge(me, u, v, c):
        if u == v:
            raise ValueError("can't add edge to self")
        if me.getEdge(u,v) is not None:
            return
        edge = Edge(me.nodes[u],me.nodes[v],c,len(me.edgeList))
        me.nodes[u].edges.append(edge)
        me.nodes[v].rEdges.append(edge)
        me.edgeList.append(edge)
    def getEdge(me,s,t):
        for edge in me.nodes[s].edges:
            if (edge.sink.ID == t):
                return edge
        return None
    def mapDepth(me,start="S"):
        from collections import deque
        for nID in me.nodes:
            me.nodes[nID].depth = None
        explored = set(start)
        queue = deque([(start,0)])
        while (len(queue) > 0):
            at,depth = queue.popleft()
            curNode = me.nodes[at]      
            curNode.depth = depth
            for edge in curNode.edges:
                toID = edge.sink.ID
                if toID not in explored:
                    queue.append((toID ,depth+1))
                    explored.add(toID)
        
    def shortestPath(me,start,end):
  
        H = Heap(len(me.nodes))
        prior = {}
        for node in me.nodes:
            me.nodes[node].cost = None
            prior[node] = None
            H.insert(me.nodes[node])
            me.nodes[node].exp = False
        H.changeKey(start,0)
        Done = False         
        while (H.openPaths()):
            u = H.extractMin()
            u.exp = True
            #print("%s(%.2f)"%(u.ID,u.cost),end=",")
            if (u.ID == end):
                break
            for edge in u.edges:
                v = edge.sink
                if not (v.exp):
                    dist = u.cost + edge.wholeCost()
                    if (v.cost is None) or (dist < v.cost):
                        H.changeKey(v.ID,dist)
                        prior[v.ID] = u.ID       
        path = [end]
        cur = end
        total = 0
        #print("")
        while not (cur == start):
            if (prior[cur] == None):
                return
            rear = me.getEdge(prior[cur],cur)
            if not (prior[cur] == "S"):
                total += rear.cost - rear.sink.price + rear.source.price
            if not (cur == "T"):
                me.nodes[cur].price += me.nodes[prior[cur]].cost
            cur = prior[cur]
            path.append(cur)
        path.reverse()
        '''
        print(path)
        for nodeID in me.nodes:
            print("\t",nodeID,me.nodes[nodeID].price)
        '''
        
        if (total <= me.cutoff):
            return path
        return None             
    def leastCostMatching(me,s="S",t="T"):
        me.findPrice(s,t)
        path = me.shortestPath(s,t)
        while not (path == None):
            for i in range(1,len(path)):
                e = me.getEdge(path[i-1],path[i]).invert()
            path = me.shortestPath(s,t)
        r = []
        for e in me.edgeList:
            if (e.include):
                if not ((e.source.ID == t) or (e.sink.ID == s)):
                    line = [e.sink.ID,e.source.ID]
                    r.append(line)
        return r
    def findPrice(me,s,t):
        #print("\nCall")
        for nodeID in me.nodes:
            me.findPriceNode(s,t,nodeID)
            #print("\t",nodeID,me.nodes[nodeID].price)
        me.nodes["S"].price = 0
        me.nodes["T"].price = 0
    def findPriceNode(me,s,t,nodeID):
        least = 0
        if (len(me.nodes[nodeID].rEdges) > 0):
            least = me.nodes[nodeID].rEdges[0].cost
            for i in range(1,len(me.nodes[nodeID].rEdges)):
                edge = me.nodes[nodeID].rEdges[i]
                if not ((edge.sink.ID == t) or (edge.source.ID == s)):
                    if (edge.cost < least):
                        least = edge.cost
        me.nodes[nodeID].price = least
    def load(me,file):
        inf = open(file,'r')
        for l in inf:
            line = l.replace("\n","")
            a = line.split("\t")
            me.addEdgeLoad(a[0],a[1],float(a[2]))
    def dump(me,file):
        outf = open(file,'w')
        for line in me.edgeList:
            if (line.include):
                outf.write(line.sink.ID+"\t"+line.source.ID+"\t"+str(line.cost)+"\n")
            else:
                outf.write(line.source.ID+"\t"+line.sink.ID+"\t"+str(line.cost)+"\n")
        outf.close()
