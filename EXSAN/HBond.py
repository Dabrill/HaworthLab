import PDBTools
import BaseWaters
import KDTree
import PoseScorer
import ZMAT
import Counter
import math
import os
import multiprocessing
recorder = None
leastEnergy = -10
useCTerminalElectrostatic = False
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
class HBond:
    MaxValue = "ERROR" 
    highRes = None
    allowableAngle = None
    IdealDistance = 1.9
    Scaling = 2
    BackboneTypes = ["H","HN","N","O","OXT","C"]
    def __init__(me,tpAtm,ligAtm,geo):
        me.tpAtm = tpAtm
        me.ligAtm = ligAtm
        me.dist = geo[0]
        me.vCos = geo[1]
        me.hCos = geo[2]
        #me.geo = -1 * me.vCos * me.hCos
        me.geo = (me.vCos * me.hCos * me.vCos * me.hCos)
        me.energyValue = None
                
    def hBondPair(me):
        return me.hBondPairTypes(me.ligAtm.atomType,me.tpAtm.atomType)
    @staticmethod
    def feedConstants(CONSTANTS):
        HBond.MaxValue = CONSTANTS["directHBond"]
        HBond.sideChainScaling = CONSTANTS.get("hBondSideFrac",0.5)
    @staticmethod
    def hBondPairTypes(one,two):
        if (one == "O" and two == "N"):
            return True
        if (one == "N" and two == "O"):
            return True
        if (one == "OXT" and two == "N"):
            return True
        if (one == "N" and two == "OXT"):
            return True
        return False
    def backboneToBackbone(me):
        return (me.tpAtm.atomType in HBond.BackboneTypes) and (me.ligAtm.atomType in HBond.BackboneTypes)
    def energy(me):
        if me.energyValue is not None:
            return me.energyValue
        Ehb_base = float(HBond.MaxValue)
        e_by_dist = Ehb_base*math.exp(-(HBond.Scaling*(me.dist-HBond.IdealDistance)**2))
        e_hb = e_by_dist * me.geo 
        if not me.backboneToBackbone():
            e_hb *= HBond.sideChainScaling
        #print(me.tpAtm,"\tGeo %3.1f  Energy: %4.0f Max %f"%(me.geo,e_hb,e_by_dist))
        me.energyValue = e_hb
        return e_hb
    def valid(me):
        return me.hBondPair() and me.properAngle()
    def properAngle(me):
        if (me.cos > 0):
            return False
        #if (me.ang < 91):
        #return False
        return True
    def __str__(entry):
        connectionType = "ToSide"
        if (entry.backboneToBackbone()):
            connectionType = "ToBack"
        vAng = math.degrees(math.acos(entry.vCos))
        hAng = math.degrees(math.acos(entry.hCos))
        tup = (entry.tpAtm,entry.ligAtm,connectionType,entry.dist,vAng,hAng,entry.geo,entry.energy())
        s = "%s\t%s\t%s\t%.3f\t%.1f\t%.1f\t%.2f\t%.2f"%tup
        return s     
def evaluatePotentialHydrogenBond(atom,contact,CA,C,O,H,Stem):
    '''
    print(atom,contact)
    print(CA,C,O,"\t",H,Stem)
    '''
    hDist = H.distance(O)
    hCos = Stem.location.cosThreePoints(H.location,O.location)
    oAng = C.angle(O,H)
    V1 = H.location.difference(Stem.location)                
    V2 = O.location.difference(C.location)
    if (CA is None):
        pnSin = hCos
    else:
        V3 = C.location.difference(CA.location)
        planeNormal = V2.crossProduct(V3)
        pnSin = math.sin(math.radians(planeNormal.angle(V1)))
        #planeScore = math.abs(planeNormal.angle(V1) / 90.0)
    vectorCos = V1.dotProduct(V2) / (V1.magnitude() * V2.magnitude())
    hb = HBond(contact,atom,(hDist,pnSin,hCos))
    hb.good = False
    hb.validGeo = True
    '''
    print("\n")
    print(atom,contact)
    print(vectorCos,oAng,hCos)
    '''
    if (vectorCos < 0) and (oAng > 100) and (hCos < 0):
        hb.validGeo = True
        energy = hb.energy()
        if (energy < leastEnergy):
            hb.good = True
    '''
    print(atom,contact,hb.good)
    print("\t",vectorCos,oAng,hCos)
    print("\t",H,Stem,O,C)
    '''
    return hb
def findAcceptorBonds(donor,H,acceptors,tRes,distCutoff):
    ret = []
    aPossible = acceptors.radiusSearch(donor.location,distCutoff,True)
    for contact,dist in aPossible:
        if (dist < distCutoff):
            #print("\t%s"%contact)
            O = contact
            stemName = BaseWaters.getAcceptorStem(contact.residueType,contact.atomType)
            planeName = BaseWaters.getPlaneStem(contact.residueType,contact.atomType)
            C = tRes[contact.chain][contact.residueNumber][stemName]
            if (planeName is None):
                CA = None
            else:
                CA = tRes[contact.chain][contact.residueNumber][planeName]
            hb = evaluatePotentialHydrogenBond(donor,contact,CA,C,O,H,donor)
            if (hb.good):
                ret.append(hb)
    return ret
def findDonorBonds(acceptor,stem,plane,donors,tRes,distCutoff):
    ret = []
    dPossible = donors.radiusSearch(acceptor.location,distCutoff,True)
    for contact,dist in dPossible:
        if (dist < distCutoff):
            H = contact
            donorName = BaseWaters.getDonorStem(contact.residueType,contact.atomType)
            donor = tRes[contact.chain][contact.residueNumber][donorName]
            hb = evaluatePotentialHydrogenBond(acceptor,donor,plane,stem,acceptor,H,donor)
            if (hb.good):
                ret.append(hb)
    return ret
def CTerminalElectrostatic(residue,donors,tRes,distCutoff):
    C = residue["C"]
    CA = residue["CA"]
    O = residue["O"]
    OXT = residue["OXT"]
    ret = []
    dPossible = donors.radiusSearch(C.location,distCutoff,True)
    for contact,dist in dPossible:
        if (dist < distCutoff):
            H = contact
            donorName = BaseWaters.getDonorStem(contact.residueType,contact.atomType)
            donor = tRes[contact.chain][contact.residueNumber][donorName]
            hb = evaluatePotentialHydrogenBond(C,contact,O,CA,C,H,donor)
            '''
            print(donor)
            print("\t",hb.validGeo,hb.geo)
            '''
            if hb.validGeo:
                distO = O.distance(H)
                distOXT = OXT.distance(H)
                
                hb.dist = 0.5 * (distO + distOXT)
                hb.energyValue = None
                energy = hb.energy()
                '''
                print("\t",distO,distOXT)
                print("\t",hb.dist)
                print("\tEnergy: ",energy)
                '''
                if (energy < -300):
                    ret.append(hb)
            
    return ret
def listHBonds(CONSTANTS,pose,donors,acceptors,tRes,sortReturn = False):
    distCutoff = CONSTANTS["hBondDistCutoff"]
    angleCutoff = CONSTANTS["hBondAngle"]
    hBondEnergy = CONSTANTS["directHBond"]
    energyCutoff = CONSTANTS.get("hBondEnergyCutoff",-50)
    addTerminalH = CONSTANTS.get("addHTerminalN",False)
    ret = []
    for i in range(len(pose.residues)):
        res = pose.residues[i]
        if ("N" in res) and ((i > 0) or addTerminalH):
            N = res["N"]
            H = res.get("H",None)
            if (H is None):
                H = res.get("HN",None)
            if (H is None):
                H = pose.implyHydrogen(N.residueNumber)
            ret+=findAcceptorBonds(N,H,acceptors,tRes,distCutoff)
        if ("O" in res):
            O = res["O"]
            C = res["C"]
            CA = res["CA"]
            ret+=findDonorBonds(O,C,CA,donors,tRes,distCutoff)
            if ("OXT" in res):
                OXT = res["OXT"]
                ret+=findDonorBonds(OXT,C,CA,donors,tRes,distCutoff)
                if useCTerminalElectrostatic:
                    ret+=CTerminalElectrostatic(res,donors,tRes,distCutoff)   
    if sortReturn:
        ret = sorted(ret,key = lambda x: x.ligAtm.atomType) 
        ret = sorted(ret,key = lambda x: x.tpAtm.atomType)
        ret = sorted(ret,key = lambda x: x.ligAtm.residueNumber)
        ret = sorted(ret,key = lambda x: x.tpAtm.residueNumber)
    return ret

def textifyHBondTable(table):
    if (table == None):
        return ""
    hbondText = ""
    for entry in table:
        hbondText+=str(entry)+"\n"
    return hbondText
class PerfectHydrogenBondSeeker:
    def __init__(me,CONSTANTS,step,fxvr,targetTree,seqStep = None):
        me.step = step
        me.seqStep = seqStep
        me.fxvr = fxvr
        me.targetTree = targetTree
        me.CONSTANTS = CONSTANTS
        HBond.feedConstants(CONSTANTS)
        PRINT("Seeking H-Bonds")
        PRINT(str(fxvr.numGood())+"  poses input")
        if seqStep is None:
            me.seqStep = step


        if (len(me.step["sequence"]) == 1):
            me.replace = True
        else:
            me.replace = False
        me.reverse = me.seqStep["reverse"] == 1
        sqZmatName = ZMAT.zmatName(me.seqStep)
        me.sqZmatLoc = "%s/%s"%(me.CONSTANTS["commonFolder"],sqZmatName)
        me.zmat = ZMAT.ZMAT(me.sqZmatLoc)
        me.zmat.readReference(me.CONSTANTS)

        me.phi = None
        me.psi = None
        me.atmNums = []   
        if me.replace:
            me.fixvarTors = me.indexOfPermutatedInIntiation(fxvr,me.zmat)
            me.angleDic = {}
        else:
            me.fixvarTors = None
        for i in range(len(step["variableTorsions"])):
            tors = step["variableTorsions"][i]
            rotAtmNum = tors["atom"]
            rotAtm = me.zmat.zAtoms[rotAtmNum-4]
            if me.replace:
                me.angleDic[(rotAtm.atomType,tors["type"])] = i
            else:
                if (not me.reverse):
                    me.rotRes = rotAtm.resNum
                    
            if (tors["type"] == 0):
                
                me.atmNums.append(rotAtmNum)
                me.resToEval = rotAtm.resNum
                if (rotAtm.atomType == "C"):
                    me.phi = tors
                if (rotAtm.atomType == "N"):
                    me.psi = tors

        
        if me.replace:
            if (me.fixvarTors[0] is None):
                PRINT("No Oxygen Torsion Set")
            if (me.fixvarTors[1] is None):
                PRINT("No Nitrogen Torsion Set")
            me.newFxvr = me.fxvr.makeChildFixvar([],me.zmat)
        else:
            me.newFxvr = me.fxvr.makeChildFixvar(me.atmNums,me.zmat)
      
        me.zmatDic = {}
        for i in range(len(me.zmat.zAtoms)):
            a = me.zmat.zAtoms[i]
            if (a.resNum not in me.zmatDic):
                me.zmatDic[a.resNum] = {}
            me.zmatDic[a.resNum][a.atomType] = i


        targetByChain = PoseScorer.splitAtomsIntoChains(targetTree.myList)
        me.tRes = PoseScorer.getTargetDict(targetByChain,True) 
        acceptorList,donorList = BaseWaters.getHBDandHBA(CONSTANTS)
        me.donorTree = KDTree.KDTree.loadAtomArray(donorList)
        me.acceptorTree = KDTree.KDTree.loadAtomArray(acceptorList)

        

        me.posePep = None
        me.dumpFile = open("HBDump.txt",'w')
    @staticmethod
    def indexOfPermutatedInIntiation(fxvr,zmat):
        ret = [None,None]
        for i in range(len(fxvr.atom)):
            kind = fxvr.kind[i]
            if (kind == 0):  
                atmNum = fxvr.atom[i]
                atm = zmat.zAtoms[atmNum-4]
                if (atm.atomType == "O"):
                    ret[0] = i
                elif(atm.atomType == "N"):
                    ret[1] = i
        return ret
    @staticmethod
    def makeUniquePoseList(fxvr,replace,dontGroup = None):
        redundancyCheckDict = {}
        posesToAnalyze = []
        goodPoseList = fxvr.goodPoseList()[0]
        numGoodPoses = len(goodPoseList)
        for i in goodPoseList:
            if replace:
                if (dontGroup is None):
                    raise Exception("Must specficy which to replace")
                myBB = ""
                for j in range(len(fxvr.kind)):
                    if j not in dontGroup:
                        myBB +=str(fxvr.poses[i-1].val[j])+"."
            else:
                myBB = fxvr.bbStr(i)
            if (myBB in redundancyCheckDict):
                answerWrapper = redundancyCheckDict[myBB]
                answerWrapper.represents.append(i)
            else:
                answerWrapper = HBSeekAnswer(i)
                redundancyCheckDict[myBB] = answerWrapper
                posesToAnalyze.append(answerWrapper)
        return posesToAnalyze
    @staticmethod
    def replacePermutationBounds(tor):
        start = tor["start"]
        stop = tor["stop"]
        tor["start"] = tor.get("hbStart",start)
        tor["stop"] = tor.get("hbStop",stop)
    def getAtom(me,residue,atype):
        if me.posePep is None:
            return None
        index = me.zmatDic[residue].get(atype,None)
        if (index is None):
            return None
        return me.posePep[index]
    def getResidue(me,residue):
        if me.posePep is None:
            return None
        return me.zmatDic[residue]

        
    def analyzeReplace(me):
        torsions = []
        Ocontacts = me.listPosibilities(me.O,me.C,me.C,me.CA,1.22,120.0)
        Otors = []
        for ocon in Ocontacts:
            torCalced = me.calcTorsionFromPossibility(ocon)
            if (torCalced is not None):
                me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(ocon[0],ocon[3],torCalced[0],ocon[2]))
                Otors.append(torCalced)
        NContacts = me.listPosibilities(me.N,me.CA,me.CA,me.C,1.46,120.0)
        Ntors = []
        for ncon in NContacts:
            torCalced = me.calcTorsionFromPossibility(None,ncon)
            if (torCalced is not None):
                me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(ncon[0],ncon[3],torCalced[1],ncon[2]))
                Ntors.append(torCalced)

                
        NTorPL = me.poseLine[me.angleDic[("N",0)]]
        oindx = me.angleDic.get(("O",0),None)
        if oindx is None:
            OTorPL = 0
        else:
            OTorPL = me.poseLine[oindx]
        added = set([(OTorPL,NTorPL)]) 

        for o in Otors:
            prop = (o[0],NTorPL)
            if prop not in added:
                torsions.append(prop)
                added.add(prop)
                #print("O ",prop)
        for n in Ntors:
            prop = (OTorPL,n[1])
            if prop not in added:
                torsions.append(prop)
                added.add(prop)
                #print("N ",(OTorPL,n[1]))
                
        if (len(Ntors) > 0) and (len(Otors) > 0):
            for o in Otors:
                for n in Ntors:
                    prop = (o[0],n[1])
                    if prop not in added:
                        torsions.append(prop)
                        added.add(prop)
                        #print("B ",prop)
        return torsions
    def analyzeReverse(me):
        def getBcontacts():
            return me.listPosibilities(me.N1,me.H1,me.CA1,me.C1,1.46,120.0)
        torsions = []
        Acontacts = []
        if (me.H2 is not None):
            Acontacts = me.listPosibilities(me.H2,me.N2,me.N2,me.CA2,1.01,120.0) + me.listPosibilities(me.O1,me.C1,me.N2,me.CA2,2.23,92.0)
        for i in range(len(Acontacts)):
            aCon = Acontacts[i]
            aTorCalced = me.calcTorsionFromPossibility(aCon)
            if (aTorCalced is not None):
                me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(aCon[0],aCon[3],aTorCalced[0],aCon[2]))
                BContacts = getBcontacts()
                for bCon in BContacts:
                    bTorCalced = me.calcTorsionFromPossibility(aCon,bCon)
                    if (bTorCalced is not None):
                        me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(bCon[0],bCon[3],bTorCalced[1],bCon[2]))
                        torsions.append(bTorCalced)
                if (len(BContacts) == 0):
                    for psiAngle in range(me.psi["start"],me.psi["stop"]+1,me.psi["step"]):
                        torPermute = (aTorCalced[0],psiAngle)
                        me.regen(torPermute)
                        if not me.poseHasClash():
                            torsions.append(torPermute)


        if (len(torsions) == 0):
            for phiAngle in range(me.phi["start"],me.phi["stop"]+1,me.phi["step"]):
                torPermute = (phiAngle,None)
                me.regen(torPermute)
                permuteO = (me.O1,me.O1.location)                        
                if not me.poseHasClash(("N","H","CA")):
                    BContacts = getBcontacts()
                    for bcon in BContacts:
                        torCalced = me.calcTorsionFromPossibility(permuteO,bcon)
                        if (torCalced is not None):
                            me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(bcon[0],bcon[3],torCalced[1],bcon[2]))
                            torsions.append(torCalced)
        return torsions
    def analyzeForward(me):
        def getOcontacts():
            return me.listPosibilities(me.O2,me.C2,me.CA2,me.N2,2.37,95.0)

        serial = 0
        torsions = []
        NContacts = []
        if (me.H2 is not None):
            NContacts = me.listPosibilities(me.H2,me.N2,me.C1,me.CA1,2.05,95.0) + me.listPosibilities(me.O1,me.C1,me.C1,me.N2,1.22,120.0)
        #print(NContacts)
        for ncon in NContacts:
            NtorCalced = me.calcTorsionFromPossibility(ncon,None)
            if (NtorCalced is not None):                
                me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(ncon[0],ncon[3],NtorCalced[0],ncon[2]))
                OContacts = getOcontacts()
                for ocon in OContacts:
                    OtorCalced = me.calcTorsionFromPossibility(ncon,ocon)
                    if (OtorCalced is not None):
                        me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(ocon[0],ocon[3],OtorCalced[1],ocon[2]))
                        torsions.append(OtorCalced)
                if (len(OContacts) == 0):
                    for phiAngle in range(me.phi["start"],me.phi["stop"]+1,me.phi["step"]):
                        torPermute = (NtorCalced[0],phiAngle)
                        me.regen(torPermute)
                        if not me.poseHasClash():
                            torsions.append(torPermute)
                
        if (len(torsions) == 0):
            for psiAngle in range(me.psi["start"],me.psi["stop"]+1,me.psi["step"]):
                torPermute = (psiAngle,None)
                #print(torPermute)
                me.regen(torPermute)
                if (me.H2 is None):
                    permuteN = (me.N2,me.N2.location)
                else:
                    permuteN = (me.H2,me.H2.location)
                if not me.poseHasClash(("C","O")):
                    OContacts = getOcontacts()
                    for ocon in OContacts:
                        OtorCalced = me.calcTorsionFromPossibility(permuteN,ocon)
                        if (OtorCalced is not None):
                            me.dumpFile.write("\t%s<->%s at %2.1f for %4.0f\n"%(ocon[0],ocon[3],OtorCalced[1],ocon[2]))
                            torsions.append(OtorCalced)
        return torsions
    def analyzeAllPoses(me):
        import Fixvar
        posesToAnalyze = PerfectHydrogenBondSeeker.makeUniquePoseList(me.fxvr,me.replace,me.fixvarTors)
        uniqueSize = len(posesToAnalyze)
        PRINT("There are %d unique hbond backbones to analyze"%(uniqueSize))
        '''
        cores = min(CONSTANTS["cores"],uniqueSize)
        for i in range(cores):
            posesToAnalyze.append(None)
        '''

        countDisp = Counter.Counter(uniqueSize,"%d/%d   hb")
        serial = 0
        for i in range(len(posesToAnalyze)):
            countDisp.disp(i)
            #input("Pose %i"%i)
            poseNum = posesToAnalyze[i].represents[0]
            me.dumpFile.write("%s\n"%posesToAnalyze[i].represents)
            me.poseLine = me.fxvr.poses[poseNum-1].val
            if me.replace:
                me.posePep = me.zmat.fixvarPoseToCartesian(me.newFxvr,me.poseLine)
                me.makeFirstResidueAtomVariables()
            else:
                me.posePep = me.zmat.fixvarPoseToCartesian(me.newFxvr,me.poseLine+[0,0])
                if me.reverse:
                    me.makeReverseAtomVariables()
                else:
                    me.makeForwardAtomVariables()
            me.pdbOrder = me.zmat.arrangeAtoms(me.posePep)
            me.torEvaled = set()          
            if me.reverse:
                torsions = me.analyzeReverse()
            elif me.replace:
                torsions = me.analyzeReplace()
            else:
                torsions = me.analyzeForward()

            '''
            for t in torsions:
                outf = open("Possible%i.pdb"%serial,'w')
                for a in me.pdbOrder:
                    outf.write(a.toPDBLine()+"\n")
                outf.close()
                print(me.CA2.location)
                print(t)
                me.regen(t)
                outf = open("PossibleRegen%i.pdb"%serial,'w')
                for a in me.pdbOrder:
                    outf.write(a.toPDBLine()+"\n")
                outf.close()
                serial+=1
            '''

     
            if (me.replace):
               for t in torsions:
                   #print("Accept ",t)
                   newPose = me.poseLine.copy()
                   for j in range(2):
                       col = me.fixvarTors[j]
                       if (col is not None):
                           newPose[col] = t[j]
                   #print(newPose)
                   me.newFxvr.poses.append(Fixvar.FixvarPoses(newPose,0))
            else:
                for n in posesToAnalyze[i].represents:
                    for t in torsions:
                        takeVal = me.fxvr.poses[n-1].val
                        newPose = takeVal.copy()
                        for k in range(2):
                            newPose.append(t[k])   
                        me.newFxvr.poses.append(Fixvar.FixvarPoses(newPose,0))
        me.dumpFile.close()
    def debugDump(me):
        me.newFxvr.createApprovedFixvar("fixvarHBondSeek.out")
    def conclude(me):
        import Log
        import Fixvar
        me.newFxvr.createApprovedFixvar("fixvarHBondSeek.out")
        numCreated = me.newFxvr.numGood()
        PRINT("\t"+str(numCreated)+" poses created")
        PRINT(Log.timeString())

        fixvarfile = Fixvar.Fixvar()
        fixvarfile.setBB(ZMAT.ZMAT(me.sqZmatLoc))
        fixvarfile.createApprovedFixvar("fixvarTMD.out")

        PRINT("Combining fixvar files")
        fixvarfile.incorporate(me.newFxvr)
        PRINT("\t"+str(fixvarfile.numGood())+" poses total")
        fixvarfile.createApprovedFixvar("fixvar.out")
        return fixvarfile
        
    def listPosibilities(me,atom,atmStem,aDist,aAng,dist,ang):
        if (atom.atomType == "N"):
            contacts = me.acceptorTree.radiusSearch(aDist.location,me.CONSTANTS["hBondDistCutoff"]+dist,True)
            donor = False
        if (atom.atomType == "H"):
            contacts = me.acceptorTree.radiusSearch(aDist.location,me.CONSTANTS["hBondDistCutoff"]+dist,True)
            donor = False
        if (atom.atomType == "O"):
            contacts = me.donorTree.radiusSearch(aDist.location,me.CONSTANTS["hBondDistCutoff"]+dist,True)
            donor = True
            planeStem = me.getAtom(atom.residueNumber,"CA")
        possible = []
        for con,conDist in contacts:
            if (conDist < (0.75 + dist)):
                continue
            atom.setAtomLocationByZMAT(aDist,aAng,con,dist,ang,0)
            if donor:
                seek = BaseWaters.getDonorStem(con.residueType,con.atomType)
            else:
                seek = BaseWaters.getAcceptorStem(con.residueType,con.atomType)
                planeSeek = BaseWaters.getPlaneStem(con.residueType,con.atomType)
                if (planeSeek is None):
                    planeStem = None
                else:
                    planeStem = me.tRes[con.chain][con.residueNumber][planeSeek]
            conStem = me.tRes[con.chain][con.residueNumber][seek]
            if (atom.atomType == "N"):
                if me.reverse:
                    H = me.H1
                elif me.replace:
                    H = me.H
                else:
                    H = me.getAtom(atom.residueNumber,"H")
                H.setAtomLocationByZMAT(atom,atmStem,con,1.01,120.0,0)
                hb = evaluatePotentialHydrogenBond(atom,con,planeStem,conStem,con,H,atom)
            else:
                hb = evaluatePotentialHydrogenBond(atmStem,con,planeStem,conStem,con,atom,atmStem)
            #print(con,hb.energy())
            #print("\t%2.1f  %3.1f  %3.1f  %3.1f"%(hb.dist,hb.vCos,hb.oCos,hb.hCos))
            if (hb.good):
                possible.append((atom,atom.location,hb.energy(),con))
        return possible
    def listStemPossibilities(me,atom,atomChild,aDist,aAng,dist,ang):
        def processContact():
            if (atom.atomType == "N"):
                atomChild.setAtomLocationByZMAT(atom,aDist,con,childDist,childAng,0)
                hb = evaluatePotentialHydrogenBond(atomChild,con,conStem,con,atomChild,atom)
                if (hb.good):
                    possible.append((atom,atom.location,hb.energy(),con))            
        if (atom.atomType == "N"):
            childDist = 1.01
            childAng = 120.0
            searchDist = me.CONSTANTS["hBondDistCutoff"]+dist+childDist
            contacts = me.acceptorTree.radiusSearch(aDist.location,searchDist,True)
            donor = False
        possible = []
        for con,conDist in contacts:
            #print(con)
            atom.setAtomLocationByZMAT(aDist,aAng,con,dist,ang,0)
            if donor:
                seek = BaseWaters.getDonorStem(con.residueType,con.atomType)
            else:
                seek = BaseWaters.getAcceptorStem(con.residueType,con.atomType)
            conStem = me.tRes[con.chain][con.residueNumber][seek]
            if (atom.distance(con) >= HBond.IdealDistance):
                processContact()
            else:
                angHAC = aAng.angle(aDist,con,False)
                angRot = ang - angHAC
                rotRange = 120  + (angRot - abs(angRot))
                print("%3.0f-%3.0f = %3.0f"%(angHAC,ang,angRot))
                print(conDist)
                zeroDist = atom.distance(con)
                deg = 0

                pX = math.cos(math.radians(angRot))
                pY = math.sin(math.radians(angRot))
                invariantX = pX * math.sin(math.radians(angHAC))
                variantX = pX - invariantX
                invariantY = pY * math.cos(math.radians(angHAC))
                variantY = pY - invariantY
                
                #dX = conDist - pX math.sin(math.radians(angHAC))
                #dY = *dist
                for deg in range(0,360,30):
                    dZ = math.sin(math.radians(deg)) * dist
                    dX = conDist - (invariantX + math.cos(math.radians(deg)) * variantX) * dist
                    dY = (invariantY + math.cos(math.radians(deg)) * variantY) * dist
                    ptDist = math.sqrt(dX**2+dY**2+dZ**2)
                    
                    atom.setAtomLocationByZMAT(aDist,aAng,con,dist,ang,deg)
                    sinASC = math.sqrt(1 - atom.location.cosThreePoints(aDist.location,con.location)**2)
                    sinSAC = math.sqrt(1 - aDist.location.cosThreePoints(atom.location,con.location)**2)
                    distCalced = conDist * sinASC / sinSAC                  
                    print(deg,atom.distance(con))
                    #print(deg,atom.angle(aDist,con),(atom.angle(aDist,con)-abs(angRot))/rotRange,"%1.2f"%(distCalced-zeroDist),conDist)#,aDist.angle(atom,con),atom.angle(con,aDist))

                    #print(distCalced,atom.distance(con)))

            '''
            if (atom.atomType == "N"):
                me.H1.setAtomLocationByZMAT(me.N1,me.CA1,con,1.01,120.0,0)
                hb = evaluatePotentialHydrogenBond(atom,con,conStem,con,me.H1,atom)
            else:
                hb = evaluatePotentialHydrogenBond(atmStem,con,conStem,con,atom,atmStem)
            #print(con,hb.energy())
            #print("\t%2.1f  %3.1f  %3.1f  %3.1f"%(hb.dist,hb.vCos,hb.oCos,hb.hCos))
            if (hb.good):
                possible.append((atom,atom.location,hb.energy(),con))
            '''
        return possible
    def calcTorsionFromPossibility(me,firstPrediction,secondPrediction = None):
        if (firstPrediction is not None):
            atomOne = firstPrediction[0]
            atomOne.location = firstPrediction[1]
        else:
            torOne = None
        if (secondPrediction is not None):
            atomTwo = secondPrediction[0]
            atomTwo.location = secondPrediction[1]
        else:
            atomTwo = None
        torTwo = None
        if me.reverse:
            if (atomOne.atomType == "H"):
                me.C1.setAtomLocationByZMAT(me.N2,me.CA2,atomOne,1.35,120.0,-180)
                me.CA1.setAtomLocationByZMAT(me.C1,me.N2,me.CA2,1.51,120.0,-180)
                me.O1.setAtomLocationByZMAT(me.C1,me.N2,me.CA1,1.22,120.0,-180)            
            if (atomOne.atomType == "O"):  
                me.C1.setAtomLocationByZMAT(atomOne,me.N2,me.CA2,1.22,32.0,180)
                if (me.H2 is not None):
                    me.H2.setAtomLocationByZMAT(me.N2,me.CA2,me.C1,1.01,120.0,-180)
                me.CA1.setAtomLocationByZMAT(me.C1,me.N2,me.CA2,1.51,120.0,-180)
            #me.N1.setAtomLocationByZMAT(me.CA1,me.C1,me.N2,1.46,109.42,-180)
            #me.H1.setAtomLocationByZMAT(me.N1,me.CA1,me.C1,1.01,120.0,0)
            torOne = int(me.C2.torsion(me.CA2,me.N2,me.C1))
            if (atomTwo is not None):
                torTwo = int(me.N2.torsion(me.C1,me.CA1,atomTwo))
        elif me.replace:
            cDistPL = me.poseLine[me.angleDic[("C",2)]] / 10.0
            cAngPL = me.poseLine[me.angleDic[("C",1)]]
            cTorPL = me.poseLine[me.angleDic[("C",0)]]
            caAngPL = me.poseLine[me.angleDic[("CA",1)]]
            caTorPL = me.poseLine[me.angleDic[("CA",0)]]
            me.C.setAtomLocationByZMAT(me.tpD,me.tpA,me.tpT,cDistPL,cAngPL,cTorPL)
            me.CA.setAtomLocationByZMAT(me.C,me.tpD,me.tpA,1.51,caAngPL,caTorPL)
            if (me.OXT is not None):
                me.OXT.setAtomLocationByZMAT(me.C,me.CA,me.O,1.22,120.0,180.0)
            torOne = int(me.O.torsion(me.C,me.CA,me.tpD))
            if (atomTwo is None):
                me.N.setAtomLocationByZMAT(me.CA,me.C,me.tpD,1.51,109.42,me.poseLine[me.angleDic[("N",0)]])
            else:
                torTwo = int(me.tpD.torsion(me.C,me.CA,atomTwo))
            if (me.H is not None):
                me.H.setAtomLocationByZMAT(me.N,me.CA,me.C,1.01,120.0,-180)
        else:
            if (atomOne.atomType == "H"):
                me.N2.setAtomLocationByZMAT(me.C1,me.CA1,atomOne,1.35,120.0,0)
                me.O1.setAtomLocationByZMAT(me.C1,me.CA1,me.N2,1.22,120.0,-180)                 
            if (atomOne.atomType == "O"):
                me.N2.setAtomLocationByZMAT(me.C1,me.CA1,atomOne,1.35,120.0,-180)
                if (me.H2 is not None):
                    me.H2.setAtomLocationByZMAT(me.N2,me.C1,me.O1,1.01,120.0,-180)
            me.CA2.setAtomLocationByZMAT(me.N2,me.C1,me.CA1,1.46,120.0,180.0)
            torOne = int(me.N2.torsion(me.C1,me.CA1,me.N1))
            if (atomTwo is not None): 
                me.C2.setAtomLocationByZMAT(me.CA2,me.N2,atomTwo,1.51,109.4,0)
                torTwo = int(me.C2.torsion(me.CA2,me.N2,me.C1))
        return me.torValidation((torOne,torTwo))
    def torValidation(me,tor):
        SeekGrain = 5
        Cutoff = -45
        if tor in me.torEvaled:
            return None
        me.torEvaled.add(tor)
        offset = SeekGrain
        proposedTor = None

        if me.reverse:
            doIgnore = ("N",)
        elif not me.replace:
            doIgnore = ("C","O")
        while me.poseHasClash():
            if tor[1] is None:
                proposedTor = (tor[0]+offset,None)
            else:
                deg = 0
                proposedTor = (tor[0],tor[1]+offset)
            me.regen(proposedTor)
            if (offset > 0):
                offset*=-1
            else:
                offset-=SeekGrain
            if (offset <= Cutoff):
                return None
        if proposedTor is not None:
            if proposedTor not in me.torEvaled:
                if tor[1] is None:
                    me.torEvaled.add(proposedTor[0])
                else:
                    me.torEvaled.add(proposedTor)
                return me.torsInBounds(proposedTor)
        return me.torsInBounds(tor)
    def torsInBounds(me,tor):       
        torOne = None
        if (tor[0] is not None):
            torOne = tor[0] % 360
        torTwo = None
        if (tor[1] is not None):
            torTwo = tor[1] % 360
        
            
        if (me.replace):
            if (torTwo is not None) and ((torTwo < me.psi["start"]) or (torTwo > me.psi["stop"])):
                return None
        elif (me.reverse):
            if (torOne is not None) and ((torOne < me.phi["start"]) or (torOne > me.phi["stop"])):
                return None
            if (torTwo is not None) and ((torTwo < me.psi["start"]) or (torTwo > me.psi["stop"])):
                return None
        else:
            if (torOne is not None) and ((torOne < me.psi["start"]) or (torOne > me.psi["stop"])):
                return None
            if (torTwo is not None) and ((torTwo < me.phi["start"]) or (torTwo > me.phi["stop"])):
                return None
        return (torOne,torTwo)
            
                
            
    def makeForwardAtomVariables(me):
        r = me.rotRes
        b = me.rotRes - 1
        me.O2 = me.getAtom(r,"O")
        me.CA2 = me.getAtom(r,"CA")
        me.N2 = me.getAtom(r,"N")
        me.C2 = me.getAtom(r,"C")
        me.H2 = me.getAtom(r,"H")
        me.CA1 = me.getAtom(b,"CA")
        me.N1 = me.getAtom(b,"N")
        me.C1 = me.getAtom(b,"C")
        me.O1 = me.getAtom(b,"O")
        me.H1 = me.getAtom(b,"H")

    def makeReverseAtomVariables(me):
        me.CA2 = me.getAtom(2,"CA")
        me.N2 = me.getAtom(2,"N")
        me.C2 = me.getAtom(2,"C")
        me.H2 = me.getAtom(2,"H")
        me.CA1 = me.getAtom(1,"CA")
        me.N1 = me.getAtom(1,"N")
        me.C1 = me.getAtom(1,"C")
        me.O1 = me.getAtom(1,"O")
        me.H1 = me.getAtom(1,"H")
    def makeFirstResidueAtomVariables(me):
        me.CA = me.getAtom(1,"CA")
        me.N = me.getAtom(1,"N")
        me.C = me.getAtom(1,"C")
        me.O = me.getAtom(1,"O")
        me.H = me.getAtom(1,"H")
        me.OXT = None
        last = me.posePep[len(me.posePep)-1]
        if (last.atomType == "OXT"):
            me.OXT = last
        me.tpD = me.zmat.dist
        me.tpA = me.zmat.ang
        me.tpT = me.zmat.tor
        
    def regen(me,tor):
        if me.reverse:
            me.C1.setAtomLocationByZMAT(me.N2,me.CA2,me.C2,1.35,120.0,tor[0])
            if (me.H2 is not None):
                me.H2.setAtomLocationByZMAT(me.N2,me.CA2,me.C1,1.01,120.0,-180)
            me.CA1.setAtomLocationByZMAT(me.C1,me.N2,me.CA2,1.51,120.0,-180)
            me.O1.setAtomLocationByZMAT(me.C1,me.N2,me.CA1,1.22,120.0,-180)
            if (tor[1] is not None):
                me.N1.setAtomLocationByZMAT(me.CA1,me.C1,me.N2,1.46,109.42,tor[1])
                me.H1.setAtomLocationByZMAT(me.N1,me.CA1,me.C1,1.01,120.0,-180)
        elif me.replace:
            cDistPL = me.poseLine[me.angleDic[("C",2)]] / 10.0
            cAngPL = me.poseLine[me.angleDic[("C",1)]]
            cTorPL = me.poseLine[me.angleDic[("C",0)]]
            caAngPL = me.poseLine[me.angleDic[("CA",1)]]
            caTorPL = me.poseLine[me.angleDic[("CA",0)]]
            nTorPL = me.poseLine[me.angleDic[("N",0)]]
            me.C.setAtomLocationByZMAT(me.tpD,me.tpA,me.tpT,cDistPL,cAngPL,cTorPL)
            me.CA.setAtomLocationByZMAT(me.C,me.tpD,me.tpA,1.46,caAngPL,caTorPL)
            if (tor[1] is None):
                me.N.setAtomLocationByZMAT(me.CA,me.C,me.tpD,1.51,109.42,nTorPL)
            else:
                me.N.setAtomLocationByZMAT(me.CA,me.C,me.tpD,1.51,109.42,tor[1])
            me.O.setAtomLocationByZMAT(me.C,me.CA,me.tpD,1.51,120.0,tor[0])
            if (me.OXT is not None):
                me.OXT.setAtomLocationByZMAT(me.C,me.CA,me.O,1.22,120.0,180.0)
            if (me.H is not None):
                me.H.setAtomLocationByZMAT(me.N,me.CA,me.C,1.01,120.0,-180)
        else:
            me.N2.setAtomLocationByZMAT(me.C1,me.CA1,me.N1,1.35,120.0,tor[0])
            me.CA2.setAtomLocationByZMAT(me.N2,me.C1,me.CA1,1.46,120.0,-180)
            if (me.H2 is not None):
                me.H2.setAtomLocationByZMAT(me.N2,me.C1,me.CA2,1.01,120.0,-180)
            me.O1.setAtomLocationByZMAT(me.C1,me.CA1,me.N2,1.22,120.0,-180)
            if (tor[1] is not None):         
                me.C2.setAtomLocationByZMAT(me.CA2,me.N2,me.C1,1.51,120.0,tor[1])
                me.O2.setAtomLocationByZMAT(me.C2,me.CA2,me.N2,1.22,120.0,-180)

            
            
         
            
    def poseHasClash(me,ignore = None):
        clashLimit = me.CONSTANTS["clashLimit"]
        for key in me.getResidue(me.resToEval):
            if (key == "H"):
                continue
            if (ignore is not None) and(key in ignore):
                continue
            atm = me.getAtom(me.resToEval,key)
            neigh, dist = me.targetTree.nearestNeighbor(atm.location,True)
            if (dist < clashLimit):
                return True
        return False
        

    
class HBSeekAnswer:
    serial = 0
    def __init__(me,poseNum):
        me.represents = [poseNum]
        me.answer = None
        HBSeekAnswer.serial+=1
        me.serial = HBSeekAnswer.serial
