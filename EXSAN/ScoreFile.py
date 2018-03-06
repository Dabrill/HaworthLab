from Heap import Heap
class ScoreFile:
    def __init__(me,allPoses,filename):
        me.outf = open(filename,'w')
        me.allPoses = allPoses
        me.at = 0
        me.received = Heap(50,key = lambda x:x[0])
    def feedPose(me,data):
        if (me.at >= len(me.allPoses)):
            print("Too many poses")
        endBracket = data.find("]")
        number = int(data[1:endBracket])
        if (number == me.allPoses[me.at]):
            me.outf.write(data)
            me.at+=1
            minPose = me.myMin()
            while (me.at < len(me.allPoses)) and (minPose == me.allPoses[me.at]):
                pose = me.received.extractMin()
                me.outf.write(pose[1])
                me.at+=1
                minPose = me.myMin()
        else:
            me.received.insert((number,data))
    def myMin(me):
        minPose = me.received.getMin()
        if (minPose is not None):
            return minPose[0]
        return None
    def end(me):
        me.outf.write("[END]")
        me.outf.close()
'''
sf = ScoreFile([1,2,3,4,5],"D:/Trump/BE9/Delta_900_fix/test.txt")
sf.feedPose("[5]\nE\nE\n")
sf.feedPose("[3]\nC\nC\n")
sf.feedPose("[1]\nA\nA\n")
sf.feedPose("[4]\nD\nD\n")
sf.feedPose("[2]\nB\nB\n")
sf.end()
'''
