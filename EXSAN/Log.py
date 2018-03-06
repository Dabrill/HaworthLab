import time
import math
import os
import shutil
class Log:
    def __init__(self,file = None, permLogLoc = None):
        if file is None:
            self.stem = "log_"+dateString()
            self.ending = ".txt"
            file = self.stem+self.ending
        else:
            self.stem, self.ending = os.path.splitext(file)
        self.filename = file
        self.out = open(file,'w')
        self.open = True
        self.permLogLoc = permLogLoc
        self.top = os.path.abspath('')
    def setPermLogLoc(self,loc):
        self.permLogLoc = loc
    def write(self,string):
        self.out.write(string)
    def close(self):
        self.out.close()
        try:
            os.remove(self.top+"/"+self.stem+"_current"+self.ending)
        except:
            pass
        self.open = False
    def copyLog(self):
        self.out.close()
        copyFilename = self.top+"/"+self.stem+"_current"+self.ending
        bup = open(copyFilename,'w')
        soFar = open(self.top+"/"+self.filename,'r')
        bup.write(timeString()+"\n")
        for line in soFar:
            bup.write(line)
        self.out = open(self.top+"/"+self.filename,'a')
        bup.close()
        soFar.close()
        if self.permLogLoc is not None:
            shutil.copy(copyFilename,permLogLoc+"/"+self.stem+self.ending)
class QuantityLog(Log):
    DIST_CULL = "distCull"
    EXTEND_CULL = "extendCull"
    TMD = "TMD"
    WAT_CULL = "waterCull"
    ATM_CULL = "atomCull"
    BB_CULL = "backboneCull"
    HBSEEK = "HB Seek"
    HEADER = "Seq\tAction\tCount\tDate\tRemark\n"
    class Entry:
        def __init__(obj,s,a,c,o,e):
            obj.seq = s
            obj.action = a
            obj.count = c
            obj.time = dateString()
            obj.order = o
            obj.extra = e
        def toLine(obj):
            ret = obj.seq+"\t"+obj.action+"\t"+str(obj.count)+"\t"+obj.time+"\t"
            if (obj.extra is not None):
                ret+=obj.extra
            return ret+"\n"
    def __init__(me,prior = None, seq = None, permLogLoc = None):
        me.table = []    
        me.curSeq = seq
        #raise NotImplementedError
        me.filename = "PoseQuantity.tsv"
        try:
            read = open("PoseQuantity.tsv",'r')
        except:
            read = None
        #print(read.readline())
        if (prior is not None) and (read is not None):
            pastHeader = False
            for line in read:
                if (pastHeader):
                    ary = line.split("\t")
                    extra = None
                    if (len(ary) > 4):
                        extra = ary[4].replace("\n","")
                    E = me.Entry(ary[0],ary[1],ary[2],len(me.table),extra)
                    E.time = ary[3]
                    me.table.append(E)
                pastHeader = True
            me.alignStep(prior)
            read.close()
            me.remakeFile()     
        else:
            Log.__init__(me,"PoseQuantity.tsv")
            me.open = True
            #me.write(timeString()+"\n")
            me.write(QuantityLog.HEADER)         
            me.curSeq = None
        me.open = True
        me.stem = "PoseQuantity"
        me.ending = ".tsv"
        me.permLogLoc = permLogLoc
    def remakeFile(me):
        me.out = open("PoseQuantity.tsv",'w')
        me.write(QuantityLog.HEADER)
        for e in me.table:
            me.out.write(e.toLine())
    def alignStep(me,step):
        rightSeq = False
        for i in range(len(me.table)):
            e = me.table[i]
            if (e.seq == me.curSeq):
                rightSeq = True
                command = step["type"]
                if (command == e.action):
                    #print("@Aligned Quant File by exact match")
                    #print(command+" "+e.action)
                    #print(e.seq+" "+me.curSeq)
                    me.table = me.table[:i]
                    return
            elif(rightSeq):
                #print("@Subsequent line found")
                me.table = me.table[:i]
                return
        #print("@No match")
    def setSequence(me,seq):
        me.curSeq = seq
    def recordStep(me,action,count,extra = None):
        
        E = me.Entry(me.curSeq,action,count,len(me.table),extra)
        me.table.append(E)
        me.write(E.toLine())
        me.copyLog()
    def lastNumberPoses(me):
        latest = me.table[len(me.table)-1]
        return latest.count
    def dump(me):
        s = ""
        for e in me.table:
            s+=e.toLine()
        return s
'''
        #print(E.toLine())
        g = time.time()
        while (time.time() < g + 2):
            pass
'''        
        
            
        
startTime = time.time()
def dateString():
    return time.strftime("%m.%d.%y_%H.%M")
def timeString():
    l = time.strftime("%H:%M:%S")
    timeDifference = int(time.time()-startTime)
    s = secondsToHMS(timeDifference)
    t = l+" ("+s+")"
    return t
def secondsToHMS(secs):
    return "%i:%02i:%02i"%(secs/3600,(secs%3600)/60,secs%60)
