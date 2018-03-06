from PDBTools import Atom,Peptide
import math
import copy
import Vector
class WaterPermutations:
    def __init__(me,size):
        me.table = [None] * size
        me.size = size
    @staticmethod
    def readCWF(fileIn = "Water.cwf"):
        inf = open(fileIn,'rb')
        line = inf.readline().decode('utf-8')
        pSig = PeptidePermutations.POSE_NUM_SIGNAL
        poseCount = int(line[len(pSig):])
        database = WaterPermutations(poseCount)
        for pose in range(poseCount):
            precon = inf.read(2)
            waterCount = int.from_bytes(precon, byteorder='big', signed=False)
            #print("WC: ",waterCount)
            byteCount = waterCount * 27
            database.table[pose] = inf.read(byteCount)
        return database
    def saveFile(me,fileOut):                  
        outf = open(fileOut,'wb')
        sizeStr = PeptidePermutations.POSE_NUM_SIGNAL+str(me.size)+"\n"
        outf.write(sizeStr.encode('utf-8'))
        for i in range(me.size):
            if me.table[i] is None:
                watCount = 0
            else:
                watCount = int(len(me.table[i]) / 27)
            outf.write(watCount.to_bytes(2, byteorder='big', signed=False))
            if me.table[i] is not None:
                outf.write(me.table[i])
        outf.close()
    def feedPose(me,pose):
        seek = pose[2] - 1
        try:
            me.table[seek] = pose[1]
        except:
            print("PoseNum ",pose[2])
            print("Tab ",len(me.table))
    def getCoord(me,pose,atom,coord):
        b = me.coordBit(pose,atom,coord)
        #print(b)
        val = int.from_bytes(b, byteorder='big', signed=True)
        #print(val)
        return val
    def getAtomAry(me,pose,atom):
        return[me.getCoord(pose,atom,0),me.getCoord(pose,atom,1),me.getCoord(pose,atom,2)]
    def coordBit(me,pose,atom,coord):
        indx = pose - 1
        poseBytes = me.table[indx]
        startByte = 9 * atom + 3 * coord
        return poseBytes[startByte:startByte+3]
    def makeWaterAtom(me,pose,atom,hydrogens):
        a = Atom("")
        if (hydrogens):
            a.atomNumber = atom+1
            if (atom%3 == 0):
                a.atomType = "O"
            elif (atom%3 == 1):
                a.atomType = "H1"
            elif (atom%3 == 2):
                a.atomType = "H2"
            a.chain = "I"
        else:
            a.atomNumber = int(atom/3)+1  
            a.atomType = "O"
            a.chain = "W"            
        a.residueType = "WAT"
        a.residueNumber = int(atom/3)+1
        temp = [None] * 3
        for coord in range(3):
            val = me.getCoord(pose,atom,coord)
            f = PeptidePermutations.intToFloat(val)
            temp[coord] = f
        a.location = Vector.Vector3D(temp)
        #a.location = Vector.Vector3D(me.getAtomAry(pose,atom))
        return a
    def grabPose(me,pose,hydrogens=True):
        indx = pose - 1
        ret = []
        data = me.table[indx]
        atomCount = int(len(data)/9)
        for atom in range(atomCount):
            if hydrogens or (atom % 3 == 0):
                a = me.makeWaterAtom(pose,atom,hydrogens)
                ret.append(a)
        return ret
    def grabPoseString(me,pose,hydrogens=True,array=False):
        if (array):
            ret = []
        else:
            ret = ""
        indx = pose - 1
        #print("Schwifty-Five")
        data = me.table[indx]
        atomCount = int(len(data)/9)
        for atom in range(atomCount):
            if hydrogens or (atom % 3 == 0):
                a = me.makeWaterAtom(pose,atom,hydrogens)
                if (array):
                    ret.append(a.toPDBLine()+"\n")
                else:
                    ret+=a.toPDBLine()+"\n"
        return ret
    def grabPoseOandH(me,pose):
        indx = pose - 1
        retO = []
        retH = []
        data = me.table[indx]
        atomCount = int(len(data)/9)
        for atom in range(atomCount):
            a = me.makeWaterAtom(pose,atom,True)
            if (atom % 3 == 0):             
                retO.append(a)
            else:
                retH.append(a)
        return retO,retH
class PeptidePermutations:
    POSE_NUM_SIGNAL = "Number of Poses:"
    def __init__(me,example,size):
        me.size = size
        me.template = example
        me.atomNum = len(example.atoms)
        me.linked = False
        me.table = bytearray(me.atomNum * 9 * size)
        me.poseGood = [True] * size
        me.bytePerPose = 9 * me.atomNum
        me.max = [me.floatToIntCoords(example.atoms[0].location.dims[0]),me.floatToIntCoords(example.atoms[0].location.dims[1]),me.floatToIntCoords(example.atoms[0].location.dims[2])]
        me.min = [me.floatToIntCoords(example.atoms[0].location.dims[0]),me.floatToIntCoords(example.atoms[0].location.dims[1]),me.floatToIntCoords(example.atoms[0].location.dims[2])]
    @staticmethod
    def floatToIntCoords(f):
        s = "%6.3f"%(f)
        return int(s.replace(".",""))
    @staticmethod
    def pdbToByteCoords(line):
        if (len(line) < 3) or ("TER" in line):
            return None            
        try:
            return (int(line[30:38].replace(".","")).to_bytes(3, byteorder='big', signed=True)) + (int(line[38:46].replace(".","")).to_bytes(3, byteorder='big', signed=True)) + (int(line[46:54].replace(".","")).to_bytes(3, byteorder='big', signed=True))
        except:
            print(line)
            raise Exception("Couldn't convert")
    @staticmethod
    def readPDB(fileIn,rejectH = True):
        inf = open(fileIn,'r')
        line = inf.readline()
        templateAry = []
        while not (line[:3] == "TER"):
            if (line[13:14] != "H") or not rejectH:
                templateAry.append(Atom(line))
            line = inf.readline()
            if (len(line) < 3):
                break
        poseCount = 1
        line = inf.readline()
        while not (line == ""):
            if (line[:3] == "TER"):
                poseCount+=1
            line = inf.readline()
        template = Peptide(templateAry,1)
        template.makeDictionary()
        inf.close()

        ATOM_COUNT = len(template.atoms)
        
        database = PeptidePermutations(template,poseCount)
        pos = 0
        inf = open(fileIn,'r')
        line = inf.readline()
        while not (line == ""):
            if (line[:3] != "TER"):
                if (line[13:14] != "H") or not rejectH:
                    ba = PeptidePermutations.pdbToByteCoords(line)
                    for i in range(9):
                        pose = int(pos / database.bytePerPose)
                        atomNum = int((pos % database.bytePerPose) / 9)
                        #print(pos,pose,atomNum,i)
                        database.table[pos] = ba[i]
                        pos+=1 
                        #print(pos,pose,atomNum,i)
            line = inf.readline()                   
        return database
    @staticmethod
    def readCPF(fileIn="all.cpf"):
        inf = open(fileIn,'rb')
        line = inf.readline().decode('utf-8')
        templateAry = []
        pSig = PeptidePermutations.POSE_NUM_SIGNAL
        while not (line[:len(pSig)] == pSig):
            templateAry.append(Atom(line))
            line = inf.readline().decode('utf-8')
        template = Peptide(templateAry,1)
        template.makeDictionary()
        poseCount = int(line[len(pSig):])
        database = PeptidePermutations(template,poseCount)
        database.table = inf.read()
        return database
    def saveFile(me,fileOut):                  
        outf = open(fileOut,'wb')
        for a in me.template.atoms:
            atmStr = a.toPDBLine()+"\n"
            outf.write(atmStr.encode('utf-8'))
        sizeStr = PeptidePermutations.POSE_NUM_SIGNAL+str(me.size)+"\n"
        outf.write(sizeStr.encode('utf-8'))     
        outf.write(me.table)
        outf.close()
    def savePDB(me,fileOut):
        goods = me.goodPoseList()[0]
        outf = open(fileOut,'w')
        for i in goods:
            #print(i)
            outf.write(me.getPosePDB(i))
        outf.close()
    def feedPose(me,pose):
        indx = pose[2] - 1
        pos = me.bytePerPose * indx
        coords = pose[0]
        if (me.bytePerPose != len(coords)):
            raise Exception("Wrong number of bytes at pose "+str(pose[2]))
        for i in range(me.bytePerPose):
            me.table[pos+i] = coords[i]
    def feedPoseIntArray(me,poseData,poseNum):
        indx = poseNum - 1
        pos = me.bytePerPose * indx
        offset = 0
        for atom in range(len(poseData)):
            for coord in range(3):
                position = poseData[atom][coord]
                if (position > me.max[coord]):
                    me.max[coord] = position
                if (position < me.min[coord]):
                    me.min[coord] = position
                value = position.to_bytes(3, byteorder='big', signed=True)
                for byte in range(3):
                    me.table[pos+offset] = value[byte]
                    offset+=1
    def goodPoseList(me,useSelfRegardless = False):
        if (me.linked) and (not useSelfRegardless):
            #print("Linked!")
            return me.fixvar.goodPoseList()
        else:
            #print("By self")
            poseList = []
            includeAll = True
            for i in range(0,len(me.poseGood)):
                if (me.poseGood[i]):
                    poseList.append(i+1)
                else:
                    includeAll = False
            return poseList,includeAll
    def calcMinMax(me):
        print("Calc")
        for i in range(me.size):
            for j in range(me.atomNum):
                atom = me.getAtomAry(i+1,j)
                for d in range(3):
                    if (atom[d] > me.max[d]):
                        me.max[d] = atom[d]
                    if (atom[d] < me.min[d]):
                        me.min[d] = atom[d]
        print("Done")
    def getCoord(me,pose,atom,coord):
        return int.from_bytes(me.coordBit(pose,atom,coord), byteorder='big', signed=True)
    def getAtomAry(me,pose,atom):
        return[me.getCoord(pose,atom,0),me.getCoord(pose,atom,1),me.getCoord(pose,atom,2)]
    def getPoseAry(me,pose):
        ret = [None] * me.atomNum
        for i in range(me.atomNum):
            ret[i] = me.getAtomAry(pose,i)
        return ret       
    def coordBit(me,pose,atom,coord):
        indx = pose - 1
        startByte = me.bytePerPose * indx
        atomStart = startByte + 9 * atom
        coordStart = atomStart + 3 * coord
        return me.table[coordStart:coordStart+3]
    def atomBit(me,pose,atom):
        indx = pose - 1
        startByte = me.bytePerPose *indx
        atomStart = startByte + 9 * atom
        return me.table[atomStart:atomStart+9]
    def poseBit(me,pose):
        indx = pose - 1
        startByte = me.bytePerPose * indx
        return me.table[startByte:startByte+me.bytePerPose] 
    @staticmethod
    def intToFloat(val):
        s = str(abs(val))
        while (len(s) < 6):
            s="0"+s
        s = s[:-3]+"."+s[-3:]
        if (val < 0):
            s = "-"+s
        f = float(s)
        return f
    @staticmethod
    def applyBytesToPeptide(peptide,byt):
        i = 0
        byteNum = 0
        while (byteNum + 5 < len(byt)):
            atomNum = i//3
            byteNum = 3*i
            #print(i,atomNum,byteNum)
            atom = peptide.atoms[atomNum]

            byteSlice = byt[byteNum:byteNum+3]
            val = int.from_bytes(byteSlice, byteorder='big', signed=True)
            f = val / 1000.0
            atom.location.dims[i%3] = f 
            i+=1
            
        return peptide
    def getPeptide(me,pose):
        for atom in range(len(me.template.atoms)):
            for coord in range(3):
                val = me.getCoord(pose,atom,coord)
                #f = me.intToFloat(val)
                f = val / 1000.0
                me.template.atoms[atom].location.dims[coord] = f
                me.template.indexNum = pose
        return me.template
    def getCopiedPeptide(me,pose):
        ret = Peptide([],pose)
        for atom in range(len(me.template.atoms)):
            newAtm = me.template.atoms[atom].copy()
            for coord in range(3):
                val = me.getCoord(pose,atom,coord)
                #f = me.intToFloat(val)
                f = val / 1000.0
                newAtm.location.dims[coord] = f
            ret.atoms.append(newAtm)
        return ret
    def getCopiedPeptideFast(me,pose):
        ret = Peptide([],pose)
        bits = me.poseBit(pose)
        for atom in range(len(me.template.atoms)):
            newAtm = me.template.atoms[atom].copy()
            atomStart = 9 * atom
            for coord in range(3):
                coordStart = atomStart + 3 * coord
                val = int.from_bytes(bits[coordStart:coordStart+3], byteorder='big', signed=True)
                #f = me.intToFloat(val)
                f = val / 1000.0
                newAtm.location.dims[coord] = f
            ret.atoms.append(newAtm)
        return ret
    def getPosePDB(me,pose,addTer=True):
        return me.getPeptide(pose).makePDBString(addTer)
    def getCoordMatrixPeptideObject(me,pose):
        ret = Peptide([],pose)
        ret.atoms = me.getPoseAry(pose)
        return ret
    def getPoseStringArray(me,pose):
        pep = me.getPeptide(pose)
        s = []
        for i in range(len(me.template.atoms)):
            s.append(pep.atoms[i].toPDBLine()+"\n")
            #s.append(me.replaceAtom(i,mat[i]).toPDBLine()+"\n")
        return s
    def linkToFixvar(me,fxvr):
        me.fixvar = fxvr
        me.linked = True
        me.fixvar.poseData = me
    def getLength(me):
        return me.size
    def applySavedCull(me,filename):
        #print("Applying ",filename)
        kept = readGoodPoseFile(filename)
        me.whitelist(kept)
    def saveCullResults(me,filename):
        goods = me.goodPoseList()[0]
        writeGoodPoseFile(filename,goods)
    def numAtoms(me):
        return me.atomNum
    def whitelist(me,keep,override = False):
        #print("Before ",me.numGood())
        if (me.linked):
            me.fixvar.whitelist(keep)
        cur = 0
        for i in range(len(me.poseGood)):
            if (cur < len(keep)):
                seek = keep[cur] - 1
                if (i == seek):
                    if (override):
                        me.poseGood[seek] = True
                    else:
                        me.poseGood[seek] = me.poseGood[seek]
                    cur+=1
                else:
                    me.poseGood[i] = False
            else:
                me.poseGood[i] = False
        #print("After ",me.numGood())
    def firstIndex(me):
        for i in range(len(me.poseGood)):
            if (me.poseGood[i]):
                return (i+1)
    def allGood(me):
        if (me.linked):
            return me.goodPoseList()[1]
        for g in me.poseGood:
            if not g:
                return False
        return True
    def good(me,poseNum,override = False):
        if (me.linked):
            me.fixvar.good(poseNum)
        if (override):
            seek = poseNum - 1
            me.poseGood[seek] = True
    def bad(me,poseNum):
        if (me.linked):
            me.fixvar.bad(poseNum)
        seek = poseNum - 1
        me.poseGood[seek] = False
        if (me.poseGood[seek] != me.fixvar.poses[seek].retain):
            print("Mismatch! "+str(seek))
        else:
            pass
    def numGood(me):
        count = 0
        for g in me.poseGood:
            if(g):
                count+=1
        if (me.linked):
            otherCount = me.fixvar.numGood()
            if (count != otherCount):
                outf = open("MistatchDump.txt",'w')
                for i in range(0,len(me.fixvar.poses)):
                    if (me.poseGood[i] != me.fixvar.poses[i].retain):
                        outf.write("X "+str(i)+"\n")
                        print("X "+str(i))
                    else:
                        outf.write("V "+str(i)+" "+str(me.poseGood[i])+"\n")
                outf.close()
                raise Exception("Fixvar and PeptidePermutations count of good poses does not match\nFixvar:"+str(otherCount)+"\nPeptidePermuations:"+str(count)+"\n")
        return count
    def applyEnergyToFixvar(currentPoseData):
        import ETable
        table = ETable.EnergyTable.readTableFile("Etable.tsv")
        for i in range(table.length()):
            pose = table.getParameter("Pose",i)
            energy = table.getParameter("Energy",i)
            if (energy is None):
                currentPoseData.bad(i+1)
            else:
                currentPoseData.fixvar.setEnergy(pose,energy)
        currentPoseData.fixvar.energyInfo = True
        #back = currentPoseData.fixvar.confirmAllEnergySet()
        #if back:
            #PRINT("Fixvar has energy set")
        #else:
            #raise Exception("Fixvar did not have all energy values set")    
def fileToStringArray(file):
    fileRead = open(file,'r')
    back = []
    for line in fileRead:
        back.append(line)
    fileRead.close()
    return back   
def StringArrayToString(array):
    s = ""
    for line in array:
        s+=line
        if not line[-1:] == "\n":
            #print("|"+line[-1:]+"|")
            s+="\n"          
    return s
def newLoc(angle):
    dx = math.cos(angle)
    dy = math.sin(angle)
def padString(s,numChars,front = True):
    r = s
    while (len(r) < numChars):
        if (front):
            r = " "+r
        else:
            r = r+" "
    return r
def pause(length):
    start = time.time()
    while (time.time() < start + length):
        pass
def foundStreak(inAry,start):
    cur = start
    while (len(inAry) - 1 > cur) and (inAry[cur+1] == inAry[cur] + 1):
        cur = cur + 1
    return (cur - start)
def readGoodPoseFile(filename):
    outarray = []
    inf = open(filename,'r')
    pastHeader = False
    for line in inf:
        if (pastHeader):
            instruction = line[:1]
            if (instruction == "#"):
                outarray.append(int(line[1:]))
            else:
                parameters = line[1:].split(",")
                parameters = [int(x) for x in parameters]
                if (instruction == "+"):
                    for n in range(parameters[0],parameters[0]+parameters[1]+1):
                        outarray.append(n)
                if (instruction == "-"):
                    for n in range(parameters[0],parameters[1]+1):
                        outarray.append(n)
        else:
            pastHeader = True
    return outarray
def dumpArray(file,array):
    outf = open(file,'w')
    for n in array:
        outf.write(str(n)+"\n")
    outf.close()
def writeGoodPoseFile(filename,inAry):
    outf = open(filename,'w')
    cur = 0
    CUTOFF = 10
    outf.write(str(len(inAry))+"\n")
    while (cur < len(inAry)):
        streak = foundStreak(inAry,cur)
        if (streak == 0):
            outf.write("#"+str(inAry[cur])+"\n")
            cur+=1
        elif (streak < CUTOFF):
            outf.write("+"+str(inAry[cur])+","+str(streak)+"\n")
            cur+=streak+1
        elif (streak >= CUTOFF):
            outf.write("-"+str(inAry[cur])+","+str(inAry[cur]+streak)+"\n")
            cur+=streak+1
        else:
            raise Exception("Serious Error. Streak value not valid: "+str(streak))


