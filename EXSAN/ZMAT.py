import PDBTools
import Fixvar
import Vector
import KDTree
import multiprocessing
import HBond
from BaseWaters import isAcceptor,getAcceptorStem,getPlaneStem
commonFolder = None
def setCommonFolder(cons):
    global commonFolder
    commonFolder = cons["commonFolder"]
recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
class ZMAT_Atom:
    def __init__(me,line):
        me.fixvarNum = int(line[:6])
        me.atomType = line[7:10].replace(" ","")
        me.distAtom = int(line[10:16])
        me.angAtom = int(line[16:22])
        me.torAtom = int(line[22:28])
        me.dist = float(line[29:36])
        me.ang = float(line[37:44])
        me.tor = float(line[45:52])
        me.unk1 = int(line[53:54])
        me.unk2 = int(line[55:56])
        me.resType = line[57:60].replace(" ","")
        me.resNum = int(line[61:67])
        me.backboneOrder = int(line[67:73])
    @staticmethod
    def blank():
        return ZMAT_Atom("     0 X       0     0     0    0.000   0.00    0.00 0 0 BAD      0     0")
    def toZMATLine(me):
        propertyTuple = (me.fixvarNum, me.atomType,me.distAtom, me.angAtom, me.torAtom, me.dist, me.ang,me.tor,me.unk1,me.unk2,me.resType,me.resNum,me.backboneOrder)
        return "%6i %-3s%6i%6i%6i%8.2f%8.2f%8.2f %1i %1i %3s %6i%6i"%propertyTuple
    def identifier(me):
        return me.resType+"-"+str(me.resNum)+" "+me.atomType
    def generateAtom(atmData,i = 0,chain = "L"):
        newAtom = PDBTools.Atom.blankAtom()
        newAtom.atomNumber = i + 1
        newAtom.atomType = atmData.atomType
        newAtom.residueType = atmData.resType
        newAtom.chain = chain
        newAtom.residueNumber = atmData.resNum
        return newAtom
def createZMAT(filepath):
    import os
    if os.path.isfile(filepath):
        return
    filename = os.path.basename(filepath)
    import MakeZMAT
    import shutil
    firstSplit = filename.split("_")
    sequence = firstSplit[0]
    secondSplit = firstSplit[1].split(".")
    number = int(extractNumber(secondSplit[0]))
    PRINT("Making ZMAT for %s starting at residue %i "%(sequence,number))
    MakeZMAT.SeqToZMAT(sequence,number,filename)
    shutil.copy(filename,"%s/%s"%(commonFolder,filename))
class ZMAT:
    def __init__(me,filename):
        splitFilename = filename.split('/')
        me.filename = splitFilename[-1]
        createZMAT(filename)
        inf = open(filename,'r')
        inf.readline()
        me.zAtoms = []
        for line in inf:
            a = ZMAT_Atom(line)
            me.zAtoms.append(a)
            #print(a.identifier())
            
        inf.close()
    def getAtomData(me,res,aType,asAtomObject = False):
        for a in me.zAtoms:
            if (a.resNum == res) and (a.atomType == aType):
                if (asAtomObject):
                    return a.generateAtom()
                return a
    def setReference(me,dist,ang,tor):
        me.dist = dist
        me.ang = ang
        me.tor = tor
    def readReference(me,CONSTANTS):
        try:
            target = PDBTools.readAllInPDB("out.pdb")
        except:
            target = PDBTools.readAllInPDB("%s/out.pdb"%CONSTANTS["commonFolder"])
        distNum = CONSTANTS["referenceAtomDistance"] - 1
        dist = target[distNum]
       
        angNum = CONSTANTS["referenceAtomAngle"] - 1
        ang = target[angNum]
        torNum = CONSTANTS["referenceAtomTorsion"] - 1
        tor = target[torNum]
        me.setReference(dist,ang,tor)
    def convertToCartesian(me,fxvr,poseNum,chain = "L"):
        fPose = fxvr.poses[poseNum].val
        return me.fixvarPoseToCartesian(fxvr,fPose,chain)
    def fixvarPoseToCartesian(me,fxvr,fPose,chain = "L"):
        cartAtoms = [None,me.dist,me.ang,me.tor]
        atomParams = []
        for i in range(len(me.zAtoms)):
            a = me.zAtoms[i]
            entry = [a.tor,a.ang,a.dist]
            atomParams.append(entry)
        #print(fxvr.atom)
        #print(fxvr.kind)
        for i in range(len(fxvr.atom)):
            atom = fxvr.atom[i] - 4
            kind = fxvr.kind[i]
            value = fPose[i]
            #print(atom,kind,value)
            if (kind == 2):
                atomParams[atom][kind] = value / 10.0
            else:
                atomParams[atom][kind] = value

                       
        for i in range(len(atomParams)):
            #print((i+4),atomParams[i])
            atmData = me.zAtoms[i]
            num = i + 4
            newAtom = atmData.generateAtom(i,chain)

            aDist = cartAtoms[atmData.distAtom]
            #print(aDist.toPDBLine())
            aAng = cartAtoms[atmData.angAtom]
            #print(aAng.toPDBLine())
            aTor = cartAtoms[atmData.torAtom]
            #print(aTor.toPDBLine())
            newAtom.setAtomLocationByZMAT(aDist,aAng,aTor,atomParams[i][2],atomParams[i][1],atomParams[i][0])
            #print(newAtom.toPDBLine()+"\n")
            cartAtoms.append(newAtom)
        return cartAtoms[4:]
    def modifyCartesian(me,fxvr,cartAtoms,atomNumbers,torsions,residue):
        atomNumbersRef = []
        
        atomParams = []
        for i in range(len(me.zAtoms)):
            if (me.zAtoms[i].resNum == residue):
                a = me.zAtoms[i]
                entry = [a.tor,a.ang,a.dist,a.fixvarNum,i]
                atomParams.append(entry)
                #print(entry)

        for j in range(len(atomNumbers)):                       
            for i in range(len(atomParams)):
                if (atomParams[i][3] == atomNumbers[j]):
                    atomParams[i][0] = torsions[j]

        #print(atomParams)
        returnAtoms = []
        for i in range(len(atomParams)):
            num = atomParams[i][4]
            #print((i+4),atomParams[i])
            newAtom = cartAtoms[atomParams[i][3]-4]
            atmData = me.zAtoms[num]
            aDist = cartAtoms[atmData.distAtom-4]
            #print(aDist.toPDBLine())
            aAng = cartAtoms[atmData.angAtom-4]
            #print(aAng.toPDBLine())
            aTor = cartAtoms[atmData.torAtom-4]
            #print(aTor.toPDBLine())
            newAtom.setAtomLocationByZMAT(aDist,aAng,aTor,atomParams[i][2],atomParams[i][1],atomParams[i][0])
            #print(newAtom.toPDBLine()+"\n")
            cartAtoms[atomParams[i][3]-4]=newAtom
            returnAtoms.append(newAtom)
        return returnAtoms
            
    def arrangeAtoms(me,cartAtoms):
        cartAtoms = sorted(cartAtoms,key=lambda a: sortVal(a.atomType))
        cartAtoms = sorted(cartAtoms,key = lambda a: a.residueNumber)
        for i in range(len(cartAtoms)):
            cartAtoms[i].atomNumber = i+1
        return cartAtoms
    def makePDBString(me,fxvr,poseNum,chain = "L",addTer=True):
        atms = me.convertToCartesian(fxvr,poseNum,chain)
        atms = me.arrangeAtoms(atms)
        s = ""
        for a in atms:
            s+=a.toPDBLine()+"\n"
        if (addTer):
            s+="TER\n"
        return s
    def createAllPDB(me,fxvr,filename="All.pdb",chain = "L"):
        outf = open(filename,'w')
        for p in range(fxvr.getLength()):
            if (p % 5000 == 0):
                print(p)
            outf.write(me.makePDBString(fxvr,p,chain))
        outf.close()
    def createOrderMap(me):
        tzAtoms = sorted(me.zAtoms,key=lambda a: sortVal(a.atomType))
        tzAtoms = sorted(tzAtoms,key = lambda a: a.resNum)
        orderMap = [None] * len(tzAtoms)
        for i in range(len(tzAtoms)):
            order = tzAtoms[i].fixvarNum-4
            #print(order,i)
            orderMap[order] = i
            #print(tzAtoms[i].identifier())
        return orderMap
    def createPeptidePermutations(me,fxvr):
        orderMap = me.createOrderMap()
        exAtms = me.convertToCartesian(fxvr,0,"L")
        exAtms = me.arrangeAtoms(exAtms)
        example = PDBTools.Peptide(exAtms,0)
        example.makeDictionary()
        pp = PDBTools.PeptidePermutations(example,fxvr.getLength())
        for pose in range(fxvr.getLength()):
            if (pose % 5000 == 0):
                print(pose, "  3D-Coord")
        #for pose in range(1):
            #print("Pose ",pose)
            ary = me.convertToIntArray(fxvr,pose,orderMap)
            #print(ary)
            pp.feedPoseIntArray(ary,pose+1)
        return pp
class ThreeDCoord_Thread(multiprocessing.Process):
    def __init__(me,zobj,fxvr,CONSTANTS,tree,inQueue,outQueue):
        multiprocessing.Process.__init__(me)
        me.zObj = zobj
        me.orderMap = me.zObj.createOrderMap()
        me.fxAtom = fxvr.atom
        me.fxKind = fxvr.kind
        me.CONSTANTS = CONSTANTS
        me.tree = tree
        me.inQueue = inQueue
        me.outQueue = outQueue

        me.NTerminalHydrogen = None
        for i in range(len(me.zObj.zAtoms)):
            a = me.zObj.zAtoms[i]
            if (a.atomType == "H") and (a.resNum == 1):
                me.NTerminalHydrogen = i
            if (a.atomType == "N") and (a.resNum == 1):
                me.NTerminal = i
        me.dummyN = me.zObj.zAtoms[me.NTerminal].generateAtom()
        if (me.NTerminalHydrogen is not None):
            me.dummyH = me.zObj.zAtoms[me.NTerminalHydrogen].generateAtom()
        
        import PoseScorer
        targetByChain = PoseScorer.splitAtomsIntoChains(tree.myList)
        me.tRes = PoseScorer.getTargetDict(targetByChain,True)
        
        
    def run(me):
        HBond.HBond.feedConstants(me.CONSTANTS)
        while(True):
            task = me.inQueue.get()
            if task is None:
                # Poison pill means shutdown
                #print ('%s: Exiting' % self.name)
                me.inQueue.task_done()
                break
            poseNum = task[0]
            pose = task[1]
            coordMatrix = me.convertToIntArray(pose)
            isGood = me.validity(coordMatrix)
            me.inQueue.task_done()
            me.outQueue.put((poseNum,coordMatrix,isGood))            
            

    def validity(me,ary):
        class Wrapper:
            def __init__(me,locary):
                import Vector
                vec = Vector.Vector3D(locary)
                me.location = vec
        def wrapIntArray(poseIntArray):
            ret = [None] * len(poseIntArray)
            for i in range(len(poseIntArray)):
                ret[i] = Wrapper(poseIntArray[i])
                ret[i].atomType = me.zObj.zAtoms[i].atomType
            return ret
        aryWrapped = wrapIntArray(ary)
        poseTree = KDTree.KDTree.loadAtomArray(aryWrapped)
        sulphurClash = me.CONSTANTS["clashLimit"]+1.0
        for atom in aryWrapped:
            neighs = poseTree.radiusSearch(atom.location,900)
            bad = False
            for n in neighs:
                if (n != atom):
                    return False
            if (atom.atomType[:1].upper() == "S"):
                atomLoc = atom.location.scalarMultiplication(.001)
                neigh = me.tree.radiusSearch(atomLoc,sulphurClash)
                for n in neigh:
                    #print("Sulphur Hit")
                    #print(n)
                    
                    if  not (n.isHydrogen()):
                        #print("Sulphur Clash")
                        return False
                
                
        return True
    def convertToIntArray(me,fPose):
        #print(len(fPose.val))
        cartAtoms = [None,me.zObj.dist.location,me.zObj.ang.location,me.zObj.tor.location]
        atomParams = []
        intArray = [None] * len(me.zObj.zAtoms)
        for i in range(len(me.zObj.zAtoms)):
            a = me.zObj.zAtoms[i]
            entry = [a.tor,a.ang,a.dist]
            atomParams.append(entry)

        for i in range(len(me.fxAtom)):
            atom = me.fxAtom[i] - 4
            kind = me.fxKind[i]
            #print(len(fPose.val),i)
            value = fPose.val[i]
            #print(atom,kind,value)
            if (kind == 2):
                atomParams[atom][kind] = value / 10.0
            else:
                atomParams[atom][kind] = value

        #print()
        #print(me.dist.toPDBLine())
                       
        for i in range(len(atomParams)):
            #print((i+4),atomParams[i])
            #print("\t"+str(orderMap[i]))
            atmData = me.zObj.zAtoms[i]
            aDist = cartAtoms[atmData.distAtom]
            #print("dist ",aDist,atmData.distAtom)
            aAng = cartAtoms[atmData.angAtom]
            #print("ang ",aAng,atmData.angAtom)
            aTor = cartAtoms[atmData.torAtom]
            #print("tor ",aTor,atmData.torAtom)
            newLoc = None
            if (me.NTerminalHydrogen == i) and (me.CONSTANTS.get("addHTerminalN",False)):
                newLoc = me.seekNTerminalHydrogen(aDist,aAng)
            if (newLoc is None):
                newLoc = Vector.Vector3D.zmatToCartesian(aDist,aAng,aTor,atomParams[i][2],atomParams[i][1],atomParams[i][0])
            
            cartAtoms.append(newLoc)
            coordLine = [None] * 3
            for j in range(3):
                val = round(newLoc.dims[j],3)
                sval = PDBTools.formatNumber(val).replace(".","")
                coordLine[j] = int(sval)
            intArray[me.orderMap[i]] = coordLine
            #print(coordLine)
            #print(newLoc.dims)
        return intArray
    def seekNTerminalHydrogen(me,N,CA):
        #print("Pose")
        contacts = me.tree.radiusSearch(N,me.CONSTANTS["hBondDistCutoff"], condition = lambda x : isAcceptor(x))
        if (len(contacts) == 0):
            return None
        possible = []
        import time
        for con in contacts:
            hLoc = Vector.Vector3D.zmatToCartesian(N,CA,con.location,1.01,120.0,0)
            me.dummyH.location = hLoc
            me.dummyN.location = N

            seek = getAcceptorStem(con.residueType,con.atomType)
            stem = me.tRes[con.chain][con.residueNumber][seek]
            planeSeek = getPlaneStem(con.residueType,con.atomType)
            if planeSeek is None:
                planeStem = None
            else:
                planeStem = me.tRes[con.chain][con.residueNumber][planeSeek] 
            hb = HBond.evaluatePotentialHydrogenBond(me.dummyN,con,planeStem,stem,con,me.dummyH,me.dummyN)
            if (hb is not None):
                possible.append((hLoc,hb.energy()))
                #print("\t",con,hb.energy())
            else:
                pass
                #print("\t",con,"\tReject",me.dummyN.distance(con))
        if (len(possible) == 0):
            return None
        best = min(possible,key = lambda x : x[1])
        #print(best[0],best[1])
        #time.sleep(2)
        return best[0]
            

        return None
        
def residueNumber(filename,atomNumber):
    inf = open(filename,'r')
    inf.readline()
    for line in inf:
        num = int(line[:6])
        if (num == atomNumber):
            return int(line[61:67])
    return None


sortValDict = {"C":0,"CA":1,"N":2,"H":3,"O":4,"CB":5,"OXT":7}
def sortVal(aID):
    if aID in sortValDict:
        return sortValDict[aID]
    return 6
    
def test():    
    dist = PDBTools.Atom("ATOM    149  CZ  ARG A  18      33.608  68.944  24.091")
    ang = PDBTools.Atom("ATOM    148  NE  ARG A  18      33.597  70.169  24.611")
    tor = PDBTools.Atom("ATOM    147  CD  ARG A  18      34.648  71.160  24.410")
    fxvr = Fixvar.Fixvar()
    stem = "D:/Trump/BE9/Common/"
    zObj = ZMAT(stem+"GSV_aa3.zmat")
    zObj.setReference(dist,ang,tor)
    print(zObj.makePDBString(fxvr,1))
    #zObj.createAllPDB(fxvr,"AllDab.pdb")
    print("done")
def testDense():    
    dist = PDBTools.Atom("ATOM    149  CZ  ARG A  18      33.608  68.944  24.091")
    ang = PDBTools.Atom("ATOM    148  NE  ARG A  18      33.597  70.169  24.611")
    tor = PDBTools.Atom("ATOM    147  CD  ARG A  18      34.648  71.160  24.410")
    fxvr = Fixvar.Fixvar()
    stem = "D:/Trump/BE9/Common/"
    zObj = ZMAT(stem+"GSV_aa3.zmat")
    zObj.setReference(dist,ang,tor)
    zObj.createPeptidePermutations(fxvr)
    #print(zObj.makePDBString(fxvr,1))
    #zObj.createAllPDB(fxvr,"AllDab.pdb")
    print("done")
def zmatObj(step):
    return ZMAT(zmatLoc(step))
def zmatName(step):
    return "%s%s.zmat"%(step["zsequence"],step["zmatSuffix"])
def zmatLoc(step):
    #print("step  ",step)
    return "%s/%s"%(commonFolder,zmatName(step))

def extractNumber(full):
    digit = ""
    for c in reversed(full):
        if c in "0123456789":
            digit = c + digit
        else:
            return digit
    return digit
'''
if __name__ == "__main__":
    print(numberify("supercalifraglistic1"))
'''
#test()
#testDense()
