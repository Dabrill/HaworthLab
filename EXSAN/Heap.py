class Heap:
    def __init__(me,n,key = None):
        me.array = [None] * n
        me.size = 0
        me.dict = {}
        me.keyFunc = key
    def applyKey(me,element):
        return me.keyFunc(element)
    def parent(me,i):
        return int((i-1)/2)
    def insert(me,item):
        if (me.size == len(me.array)):
            me.array.append(None)
        me.array[me.size] = item
        me.heapifyUp(me.size)
        me.size+=1
    def delete(me,i):
        me.array[i] = me.array[me.size-1]
        me.size-=1
        me.heapifyDown(i)   
    def changeKey(me,ID,val):
        raise NotImplemented("Cannot yet change key")
        i = me.dict[ID]
        if (i == None):
            return None
        me.array[i].cost = val
        me.heapifyUp(i)
        me.heapifyDown(i)
        return i
    def swap(me,i,j):
        temp = me.array[i]
        me.array[i] = me.array[j]
        me.array[j] = temp
        
    def heapifyUp(me,i):
        if (i > 0):
            j = me.parent(i)
            if (me.applyKey(me.array[i]) < me.applyKey(me.array[j])):
                me.swap(i,j)
                me.heapifyUp(j)
    def heapifyDown(me,i):
        top = 2 * i + 2
        if (top > me.size):
            return
        if (top < me.size):
            left = 2 * i + 1
            right = 2* i + 2
            if (me.applyKey(me.array[left]) < me.applyKey(me.array[right])):
                j = left
            else:
                j = right
        elif (top == me.size):
            j = top - 1
        if (me.applyKey(me.array[i]) > me.applyKey(me.array[j])):
            me.swap(i,j)
            me.heapifyDown(j)
    def getMin(me):
        return me.array[0]
    def extractMin(me):
        ret = me.array[0]
        me.delete(0)
        return ret
'''
h = Heap(10,key = lambda x :x[1])
h.insert(("Four",4,"4X"))
h.insert(("Eight",8,"X8"))
h.insert(("Two",2,"tuuu"))
h.insert(("Seven",7,"Xseven"))
h.insert(("Three",3,"XXX"))
print("\n\n\n\n")
print(h.extractMin())
print(h.extractMin())
print(h.extractMin())
print(h.extractMin())
print(h.extractMin())
'''
