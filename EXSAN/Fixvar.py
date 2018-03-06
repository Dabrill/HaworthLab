PERFORM_DUMP = True
def padString(s,numChars,front = True):
    r = s
    while (len(r) < numChars):
        if (front):
            r = " "+r
        else:
            r = r+" "
    return r
def splitLine(string):
    a = []
    for i in range(0,int(len(string) / 4)):
        a.append(int(string[i*4:(i+1)*4]))
    return a
class FixvarPoses:
    def __init__(me,a,num):
        me.val = a #Array of values in pose line
        me.num = num #Pose Number
        me.retain = True
        me.backboneCount = 0
        me.energy = None
    def toLine(me):
        s = ""
        for v in me.val:
            s+=padString(str(v),4)
        return s
def makePerfectCandCA(zmat,model):
    zFirst = zmat.zAtoms[0]
    xFirst = model.getSpecificAtom[0][zFirst.atomType]
    dist1 = xFirst.distance(zmat.dist)
    ang1 = xFirst.angle(zmat.dist,zmat.ang)
    tor1 = xFirst.torsion(zmat.dist,zmat.ang,zmat.tor)
    zSecond = zmat.zAtoms[1]
    xSecond = model.getSpecificAtom[0][zSecond.atomType]
    ang2 = xFirst.angle(zmat.dist,zmat.ang)
    tor2 = xFirst.torsion(zmat.dist,zmat.ang,zmat.tor)

    ret = Fixvar(None)
    ret.atom = [4,4,4,5,5,6,8]
    ret.kind = [2,1,0,1,0,0,0]
    ret.header = [padString(str(len(new.atom)),8)+" "+zmat.filename]
    
    
    i = 1
    for a in range(0,359,30):
        for b in range(0,359,60):
            me.poses.append([dist1,ang1,tor1,ang2,tor2,a,b],i)
            i+=1
    ret.save("Fixvar.out")
class Fixvar():
    def __init__(me,file = "fixvar.out"):
        me.header = []
        me.poses = []
        if file is None:
            me.atom = []
            me.kind = []
        else:
            fv = open(file,'r')
            me.header.append(fv.readline().replace("\n","").split("\t")[0])
            me.header.append(fv.readline().replace("\n",""))
            me.header.append(fv.readline().replace("\n",""))
            me.atom = splitLine(me.header[1])
            me.kind = splitLine(me.header[2])
            line = fv.readline().replace("\n","")
            while line:
                me.poses.append(FixvarPoses(splitLine(line),len(me.poses)+1))
                line = fv.readline().replace("\n","")
            fv.close()
            me.indices = {}
            for a in range(0,len(me.atom)):
                if (me.kind[a] == 0):
                    me.indices[me.atom[a]] = a
            me.current = 0
        me.energyInfo = False
        me.dumps = 1
        me.bestPose = None
    def makeChildFixvar(old,atoms,zmatName):
        new = Fixvar(None)
        for i in range(len(old.atom)):
            new.atom.append(old.atom[i])
            new.kind.append(old.kind[i])
        for a in atoms:
            new.atom.append(a)
            new.kind.append(0)
        new.header = [padString(str(len(new.atom)),8)+" "+zmatName.filename]
        return new
    def sortByAngles(self):
        for i in reversed(range(len(self.poses[0].val))):
            self.poses = sorted(self.poses,key = lambda x: x.val[i])
    def incorporate(self,other):
        if (len(self.atom) != len(other.atom)):
            print(self.atom)
            print(other.atom)
            raise Exception("Header size mismatch")
        for i in range(len(self.atom)):
            '''
            if (self.atom[i] != other.atom[i]):
                print(self.atom[i])
                print(other.atom[i])
                raise Exception("Header atom mismatch")
            '''
            if (self.kind[i] != other.kind[i]):
                raise Exception("Header kind mismatch")
        for pose in other.poses:
            self.poses.append(pose)
        if (len(self.poses) < 1):
            return
        for i in reversed(range(len(self.poses[0].val))):
            self.poses = sorted(self.poses,key = lambda x: x.val[i])

        replace = [self.poses[0]]
        for i in range(1,len(self.poses)):
            allSame = True
            for j in reversed(range(len(self.poses[i].val))):
                if (self.poses[i].val[j] != self.poses[i-1].val[j]):
                    allSame = False
                    break
            if not allSame:
                replace.append(self.poses[i])
            '''
            else:
                print(self.poses[i].val,"\n",self.poses[i-1].val)
                input("")
            '''

        for i in range(len(replace)):
            replace[i].num = i+1
        self.poses = replace
    def setBB(self,zmat):
        self.bbAtoms = []
        for i in range(len(self.atom)):
            a = self.atom[i]
            #print("At ",a,"  ",zmat.zAtoms[a-4].atomType)
            if (zmat.zAtoms[a-4].atomType in ["C","CA","N"]):
                self.bbAtoms.append(i)
                #print("\tBB ",i)
    def bbStr(self,pose):
        s = ""
        indx = pose-1
        for bb in self.bbAtoms:
            s+=str(self.poses[indx].val[bb])+"."
        return s
    def getLength(self):
        return len(self.poses)
    def getPose(self,num):
        if (num == self.poses[num-1].num):
            return self.poses[num-1]
        else:
            PRINT("FIXVAR OUT OF ORDER")
            raise Exception("FIXVAR OUT OF ORDER")
    def getEnergy(self,num):
        p = self.getPose(num)
        return p.energy
    def goodPoseList(fxvr):
        poseList = []
        includeAll = True
        for i in range(0,len(fxvr.poses)):
            if (fxvr.poses[i].retain):
                poseList.append(fxvr.poses[i].num)
            else:
                includeAll = False
        return poseList,includeAll
    def numGood(me):
        count = 0
        for p in me.poses:
            if (p.retain):
                count+=1
        return count
    def nextGood(me):
        me.current+=1
        if (me.current == len(me.poses)):
            return None
        while not (me.poses[me.current].retain):
            me.current+=1
            if (me.current == len(me.poses)):
                return None
        return me.poses[me.current].num
    def good(self,num):
        self.current+=1
    def bad(self,num):
        #print("FX "+str(num-1))
        if (num == self.poses[num-1].num):
            self.poses[num-1].retain = False
        else:
            PRINT("FIXVAR OUT OF ORDER")
            raise Exception("FIXVAR OUT OF ORDER")
        self.current+=1
    def arrange(me):
        me.poses = sorted(me.poses, key = lambda p:p.num)
        me.current = 0
        while not (me.poses[me.current].retain):
            me.current+=1
            if (me.current == len(me.poses)):
                return None
    def currentPose(me):
        return me.poses[me.current]
    def whitelist(me,keep,override = False):
        cur = 0
        me.poses = sorted(me.poses, key = lambda p:p.num)
        keep.sort()
        for i in range(0,len(me.poses)):
            #PRINT("Pose "+str(me.poses[i].num))
            if (cur < len(keep)):
                if (keep[cur] == me.poses[i].num):
                    #PRINT("Pose: "+str(keep[cur])+" "+str(me.poses[i].retain))
                    if (override):
                        me.poses[i].retain = True
                    else:
                        me.poses[i].retain = me.poses[i].retain
                    cur+=1
                else:
                    me.poses[i].retain = False
            else:
                me.poses[i].retain = False
    def createApprovedFixvar(self,filename):
        outFile = open(filename, 'w')
        kindStr = ""
        atomStr = ""
        for i in range(len(self.atom)):
            atomStr+=padString(str(self.atom[i]),4)
            kindStr+=padString(str(self.kind[i]),4)
        outFile.write(self.header[0]+"\n")
        outFile.write(atomStr+"\n")
        outFile.write(kindStr+"\n")
        for p in self.poses:
            if (p.retain):
                outFile.write(p.toLine()+"\n")
        outFile.close()
    def backboneCluster(me,atomNumbersVary):
        import random
        #me.bbAtoms is by the index number of the self.atoms array.
        #This conversion is to find out the fixvar atom number of each entry
        atomNumBB = []
        for b in me.bbAtoms:
            atomNumBB.append(me.atom[b])
        #print("BB: ",atomNumBB)

        outf = open("BackBoneDump.txt",'w')
        groupAtomsByNum = []
        groupAtomsByIndex = []
        for i in range(len(atomNumBB)):
            b = atomNumBB[i]
            if b not in atomNumbersVary:
                groupAtomsByNum.append(b)
                
        #print("By Number: ",groupAtomsByNum)
        for i in range(len(me.atom)):
            atomNum = me.atom[i]
            if atomNum in groupAtomsByNum:
                groupAtomsByIndex.append(i)            
        #print("By Index: ",groupAtomsByIndex)
        #input("")
        groups = me.generateGroups(groupAtomsByIndex)

        varyIndices = []
        for i in range(len(me.atom)):
            a = me.atom[i]
            if a in atomNumbersVary:
                varyIndices.append(i)

        for i in range(len(groups)):
            g = groups[i]
            me.cullBackboneGroup(g,varyIndices,60)
            if (outf is not None):
                groupPosenumberStart = min(g, key=lambda x: x.num)
                groupPosenumberEnd = max(g, key=lambda x: x.num)
                outf.write("\nGroup starts at pose %i and ends at pose %i\nKeep:\n"%(groupPosenumberStart.num,groupPosenumberEnd.num))  
                for p in g:
                    #masterEntry = me.poses[p.num-1]
                    if (p.retain):
                        outf.write("%i\n"%p.num)  
        '''
        outf.write("Master List:\n")
        for i in range(len(me.poses)):
            outf.write("%i %s\n"%(me.poses[i].num,str(me.poses[i].retain)))
        '''
        outf.close()
        '''
        switch = False
        if switch:
            for i in range(len(groups)):
                g = groups[i]
                back = me.cullBackboneGroup(g,varyIndices,180)
                back2 = me.cullBackboneGroup(g,varyIndices,0)
                print(i,len(g),len(back),len(back2))
        else:
            back = me.cullBackboneGroup(groups[21],varyIndices,10)
            print("Done: ",len(back))
        '''
    def cullBackboneGroup(me,group,atoms,bound = 60):
        def angularDifference(ang1,ang2):
            a = (ang1-ang2)%360
            b = (ang2-ang1)%360
            return min(a,b)
        if (len(group) == 0):
            return []
        #print("G: ",len(group))
        group = sorted(group,key=lambda p: p.energy)
        best = group[0]
        subgroup = []
        #print(best.val)
        for j in range(1,len(group)):
            pose = group[j]
            #print(pose.val)
            good = False
            for a in atoms:
                ang1 = pose.val[a]
                ang2 = best.val[a]
                good = good or (angularDifference(ang1,ang2) >= bound)
                #print(a,ang1,ang2,good)
            if good:
                subgroup.append(pose)
            else:
                #print("Cull ",me.numGood())
                culledPoseNum = pose.num
                if (me.poseData is not None):
                    me.poseData.bad(culledPoseNum)
                else:
                    me.bad(culledPoseNum)
        return [best] + me.cullBackboneGroup(subgroup,atoms)
    def generateGroups(me,indexes):
        poseList = [p for p in me.poses]
        for a in reversed(indexes):
            poseList = sorted(poseList,key=lambda p: p.val[a])

    
        groups = []
        thisGroup = []
        prototype = None
        for thisPose in poseList:
            if thisPose.retain:
                if prototype is None:
                    prototype = thisPose
                    thisGroup.append(thisPose)
                    groups.append(thisGroup)
                else:
                    same = True
                    for i in indexes:
                        same = same and (prototype.val[i] == thisPose.val[i])
                    if (same):
                       thisGroup.append(thisPose)
                    else:
                        prototype = thisPose
                        thisGroup = [thisPose]
                        groups.append(thisGroup)
        '''
        for i in range(len(groups)):
            g = groups[i]
            print(i,len(g))


        for j in range(len(groups[21])):
            l = groups[21][j]
            print(l.val)               
         '''
        return groups   
            
        
        
        
            
    def atomFreeMult(me,atomNum):
        outf = open("SideChainDump"+str(me.dumps)+".txt",'w')
        outf.write("Only Consider Free Rotation\n")
        me.dumps+=1
        import random
        indices = []
        for x in atomNum:
            indices.append(me.indices[x])
        values = []
        valuesNum = 1
        for x in atomNum:
            varities = me.possible(x)
            values.append(varities)
            outf.write(str(varities)+"\t"+str(len(varities))+"\n")
            valuesNum *= len(varities)
        print ("VN "+str(valuesNum))
        outf.write("Atoms "+str(atomNum)+"\n")
        outf.write("Possible "+str(valuesNum)+"\n")
        outf.write("  "+me.header[1]+"\n")
        if me.bestPose is not None:
            energyCutoff =  me.bestPose.energy * .7
            outf.write("Cutoff: "+str(energyCutoff))
        poseList = [p for p in me.poses]
        for a in range(0,len(me.atom)):
            i = len(me.atom) - 1 - a
            if not (i in indices):
                poseList = sorted(poseList,key=lambda p: p.val[i])           
        for i in range(0,len(poseList)- valuesNum+1):
            if (poseList[i].retain):
                first = poseList[i]
                seen = []
                allSame = True
                for j in range(0,valuesNum):
                    if(allSame):
                        current = poseList[i+j]
                        same = current.retain
                        for k in range(0,len(me.atom)):
                            if not (k in indices):
                                same = same and (first.val[k] == current.val[k])
                        if (same):
                            entry = []
                            for a in indices:
                                entry.append(current.val[a])
                            seen.append(entry)
                    allSame = same and allSame
                if (allSame) and (len(seen) == valuesNum):
                    outf.write("\n\nGroup starts at "+str(poseList[i].num)+"\n")
                    if (me.energyInfo):
                        keep = 0
                        keepE = poseList[i].energy
                        for j in range(1,valuesNum):
                            if (poseList[i+j].energy < keepE):
                                keepE = poseList[i+j].energy
                                keep = j
                    else:
                        keep = random.randint(0,valuesNum-1)
                    outf.write("Random choices "+str(keep)+"\n")
                    for j in range(0,valuesNum):
                        immune = False
                        if (me.energyInfo):
                            if (poseList[i+j].energy < energyCutoff):
                                outf.write("I")
                                immune = True
                        if not (j == keep) and not immune:
                            outf.write("X:"+poseList[i+j].toLine()+"  E:"+str(poseList[i+j].energy)+"\n")
                            thisnum = poseList[i+j].num
                            if (me.poseData is not None):
                                me.poseData.bad(thisnum)
                            else:
                                me.bad(thisnum)
                        else:
                            outf.write(">:"+poseList[i+j].toLine()+"  E:"+str(poseList[i+j].energy)+"\n")
    def atomFreeMultB(me,atomNum):
        #if (PERFORM_DUMP):
        outf = open("Dump"+str(me.dumps)+".txt",'w')
        outf.write("Consider All\n")
        me.dumps+=1
        import random
        indices = []
        for x in atomNum:
            indices.append(me.indices[x])
        values = []
        valuesNum = 1
        for x in atomNum:
            varities = me.possible(x)
            values.append(varities)
            valuesNum *= len(varities)
        print ("VN "+str(valuesNum))
        outf.write("Atoms "+str(atomNum)+"\n")
        outf.write("Possible "+str(valuesNum)+"\n")
        outf.write("  "+me.header[1]+"\n")
        if me.bestPose is not None:
            energyCutoff =  me.bestPose.energy * .85
            outf.write("Cutoff: "+str(energyCutoff)+"\n")
        poseList = [p for p in me.poses]
        for i in reversed(range(len(me.atom))):
            if not (i in indices):
                poseList = sorted(poseList,key=lambda p: p.val[i])           
        for i in range(0,len(poseList)- valuesNum+1):
            #outf.write(str(i)+"\n")
            if (poseList[i].retain):
                first = poseList[i]
                seen = []
                allSame = True
                streak = 0
                oneNotCulled = False
                for j in range(0,valuesNum):
                    #outf.write("\t"+str(i+j)+"\n")
                    if(allSame):
                        current = poseList[i+j]
                        oneNotCulled = oneNotCulled or current.retain
                        same = True
                        for k in range(0,len(me.atom)):
                            if not (k in indices):
                                same = same and (first.val[k] == current.val[k])
                        if (same):
                            entry = []
                            for a in indices:
                                entry.append(current.val[a])
                            seen.append(entry)
                            streak+=1
                    allSame = same and allSame
                #if (allSame) and (len(seen) == valuesNum):
                if (streak > 1):
                    outf.write("\n\nGroup starts at "+str(poseList[i].num)+"\n")
                    keep = None
                    keepE = None
                    for j in range(0,streak):
                        if poseList[i+j].energy is not None:
                            if keep is None:
                                keep = j
                                keepE = poseList[i+j].energy
                            elif poseList[i+j].retain and (poseList[i+j].energy < keepE):
                                keepE = poseList[i+j].energy
                                keep = j
                    if keep is not None:
                        outf.write("Best Choice "+str(keep)+"\n")
                    else:
                        keep = random.randint(0,streak-1)
                        outf.write("Random choices "+str(keep)+"\n")
                    for j in range(0,streak):
                        immune = False
                        if poseList[i+j].energy is not None:
                            if (poseList[i+j].energy < energyCutoff):
                                outf.write("I")
                                immune = True
                        if not (j == keep) and not immune:
                            outf.write("X"+str(i+j)+":"+poseList[i+j].toLine()+"  E:"+str(poseList[i+j].energy)+"\t"+str(poseList[i+j].retain)+"\n")
                            thisnum = poseList[i+j].num
                            if (me.poseData is not None):
                                me.poseData.bad(thisnum)
                            else:
                                me.bad(thisnum)
                        else:
                            outf.write(">"+str(i+j)+":"+poseList[i+j].toLine()+"  E:"+str(poseList[i+j].energy)+"\t"+str(poseList[i+j].retain)+"\n")
    def possible(me,atomNum):
        index = me.indices[atomNum]
        possible = []
        for p in me.poses:
            if not p.val[index] in possible:
                possible.append(p.val[index])
        return possible
    def setEnergy(self,pose,val):
        if (pose == self.poses[pose-1].num):
            self.poses[pose-1].energy = val
            if self.bestPose is None:
                self.bestPose = self.poses[pose-1]
            if (val < self.bestPose.energy):
                self.bestPose = self.poses[pose-1]         
        else:
            PRINT("FIXVAR OUT OF ORDER")
            raise Exception("FIXVAR OUT OF ORDER")
    def confirmAllEnergySet(me):
        for p in me.poses:
            if (p.energy is None):
                me.energyInfo = False
                return False
        me.energyInfo = True
        return True
