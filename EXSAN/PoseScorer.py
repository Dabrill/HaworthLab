import multiprocessing
import HBond
import KDTree
import PDBTools
import Network
import BaseWaters
import Permutations
import os

from enum import Enum
class Kind(Enum):
    HYDROPHOBIC = 1
    WRONG_POCKET = 2
    BULK_EXPOSED = 3
    SALT_BRIDGE = 4
    PLANAR_AMIDE = 5
    LYSINE_DP = 6
class ColumnValue(Enum):
    KIND = 0
    T_RES = 1
    L_RES = 2
    TEXT = 3
    DIST = 4
    ENERGY = 5
class OxyReplaceStatus(Enum):
    No = 0
    Side = 1
    Back = 4
def getTargetTree(cons):
    commonFolder = cons["commonFolder"]
    targetPDBFile = "%s/%s"%(commonFolder,cons["targetProtein"])

    targetByAtomsIncludingRepeats = PDBTools.readPDBSelectively(targetPDBFile,cons["targetProteinChain"])
    targetByAtoms = []
    atomsSeen = set()
    for a in targetByAtomsIncludingRepeats:
        atmId = str(a)
        if (atmId not in atomsSeen):
            atomsSeen.add(atmId)
            targetByAtoms.append(a)


    targetWithH = PDBTools.readPDBSelectively(BaseWaters.getBaseWaterFilename(cons),cons["targetProteinChain"])
    
    targetTree = KDTree.KDTree.loadAtomArray(targetByAtoms)
    targetTreeWithH = KDTree.KDTree.loadAtomArray(targetWithH)
    return targetTree,targetTreeWithH
class PoseScorer(multiprocessing.Process):
    def __init__(self, examplePose, task_queue, result_queue, hBondQueue, alignQueue,sideChainQueue, constants, zmat, baseWatersWithH):
        multiprocessing.Process.__init__(self)
        self.examplePose = examplePose
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.hBondInfo = hBondQueue
        self.alignInfo = alignQueue
        self.sideChainQueue = sideChainQueue
        self.CONSTANTS = constants
        self.Zmat = zmat
        self.baseWatersWithH = baseWatersWithH
        self.baseWaterKeeperFine = BaseWaters.BaseWaterKeeper(baseWatersWithH,2.2)
        self.baseWaterKeeperGross = BaseWaters.BaseWaterKeeper(baseWatersWithH,6.0)
        self.bwH = []
        self.baseWaters = []
        for w in self.baseWatersWithH:
            for i in range(len(w)):
                if (i % 3 == 0):
                    self.baseWaters.append(w[i])
                else:
                    self.bwH.append(w[i])
    def run(self):
        #Directory Variables
        topFolder = os.path.abspath('')

        #Target Variables
        targetTree,targetTreeWithH = getTargetTree(self.CONSTANTS)       

        targetHBondableList = []
        targetInDistance = []
        distCenter = self.Zmat.dist
        #distAllowed = self.CONSTANTS["cullRadiusMeasure"]
        HB_ALLOWED_DIST = self.CONSTANTS["hBondDistCutoff"]
        for a in targetTree.myList:
            #if (a.distance(distCenter) < distAllowed) and a.atomType in ["O","N","OXT"]:
            if a.atomType in ["O","N","OXT"]:
                targetHBondableList.append(a) 
        bondableTree = KDTree.KDTree.loadAtomArray(targetHBondableList)
        tRes = getTargetDict(splitAtomsIntoChains(targetTreeWithH.myList),True)
        '''
        for c in tRes:
            print(c)
            for r in tRes[c]:
                print("\t",r)
                '''
    
    
        #Scoring Parameters / table variables
        #BULK_VALUE = self.CONSTANTS["bulkWaterEnergy"]
        #RETURN_TO_BULK_FACTOR = self.CONSTANTS["returnToBulkFactor"]
        HYDROPHOBIC_VALUE = self.CONSTANTS["hydrophobicEnergy"]
        BULK_FACE_VALUE = self.CONSTANTS["bulkFaceEnergy"]
        SALT_VALUE = self.CONSTANTS["saltBridgeEnergy"]
        HBond.HBond.feedConstants(self.CONSTANTS)
        
        #Base Waters

        tarOxygen, tarHydrogen = BaseWaters.getHBDandHBA(self.CONSTANTS)
        possibleDonor = tarHydrogen + self.bwH
        possibleAcceptor = tarOxygen + self.baseWaters
        hTree = KDTree.KDTree.loadAtomArray(possibleDonor)
        oTree = KDTree.KDTree.loadAtomArray(possibleAcceptor)
        hTreeTPOnly = KDTree.KDTree.loadAtomArray(tarHydrogen)
        oTreeTPOnly = KDTree.KDTree.loadAtomArray(tarOxygen)
        bwTree = KDTree.KDTree.loadAtomArray(self.baseWaters)

        ignoreWater = self.CONSTANTS.get("ignoreWater",False)

        tpCenter = BaseWaters.calcAveragePoint(targetTree.myList)
        
        print(self.name+" booted")
        pose= self.examplePose
        pose.makeDictionary()
        while(True):
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print ('%s: Exiting' % self.name)
                self.task_queue.task_done()
                break
            poseNum = next_task[0]
            poseBytes = next_task[1]
            Permutations.PeptidePermutations.applyBytesToPeptide(self.examplePose,poseBytes)
            #pose = next_task[1]
            poseWatersAll = next_task[2]
            poseWaters = [None] * (len(poseWatersAll) // 3)
            for i in range(0,len(poseWatersAll) // 3):
                poseWaters[i] = poseWatersAll[i*3]
            wTree = KDTree.KDTree.loadAtomArray(poseWaters)
            poseTree = KDTree.KDTree.loadAtomArray(pose.atoms)
            results = []

            #if (poseNum % 50 == 0):
            #print("Crunching ",poseNum)

            #Hydrogen Bonds
            #foundHBonds = []
            foundHBonds = HBond.listHBonds(self.CONSTANTS,pose,hTreeTPOnly,oTreeTPOnly,tRes,True)
            #print(poseNum,len(foundHBonds))
            #import time
            #time.sleep(0.5)
            hbFileText = "["+str(poseNum)+"]\n"+textifyHBondTable(foundHBonds)
            self.hBondInfo.put(hbFileText)
            hbEnergyTotal = 0
            for hBond in foundHBonds:
                hbEnergyTotal+=hBond.energy()
            results.append((len(foundHBonds),"#HBond"))
            results.append((hbEnergyTotal,"EHBond"))


            alignStr = "["+str(poseNum)+"]\n"
            if not ignoreWater:
                #Water Alignment
                expected,neighWaters = self.expectedAndNeighborhoodWaters(poseTree,poseWaters)
                align = alignWater(bwTree,pose,poseWaters)
                #align = []

                alignStr+="Present in ligand structure:\n"
                alignDistTotal = 0
                alignCount = 0.0
                for j in range(0,len(align)):
                    alignStr+=str(j+1)+"\t"+str(align[j])
                    if (align[j] > 0):
                        left = poseWaters[j]
                        right = self.baseWaters[align[j]-1]
                        dist = left.distance(right)
                        alignDistTotal+=dist
                        alignCount+=1.0
                        alignStr+="\t"+str(round(dist,2))
                    alignStr+="\n"
                if (alignCount > 0):
                    alignAverage = alignDistTotal / alignCount
                else:
                    alignAverage = 0.0
                alignStr+="Number of alignments: "+str(int(alignCount))+"\n"
                alignStr+="Sum of alignment distances: "+str(alignDistTotal)+"\n"
                alignStr+="Average alignment distance: "+str(alignAverage)+"\n"
                alignStr+="Only Present in base:\n"
                #expected = []
                missing = missingWaters(expected,align)
                oxygenReplacement(missing,poseTree,self.CONSTANTS.get("ignoreORep",False))

                missingEnergy = 0
                alignStr+=" [Ring,Salt,Phobic,Other,Bulk,ORep,watHB,aaHB,scoreHB] -> Score\n"
                for m in missing:
                    if (self.CONSTANTS["ORepZeroMode"]):
                        oRepScore = (m.replacesOxygen is OxyReplaceStatus.Back)  * self.CONSTANTS["waterOxygenReplacement"] + m.hbScore * (m.replacesOxygen == OxyReplaceStatus.No)
                    else:
                        oRepScore = (m.replacesOxygen.value / OxyReplaceStatus.Back.value) * self.CONSTANTS["waterOxygenReplacement"] + m.hbScore * (m.replacesOxygen == OxyReplaceStatus.No)

                    nearestResidue = poseTree.nearestNeighbor(m.location).resStr()
                    #print((m.replacesOxygen is OxyReplaceStatus.Side)  * self.CONSTANTS["waterOxygenReplacement"], (m.replacesOxygen is OxyReplaceStatus.Side))
                    wScoreTotal = m.score+oRepScore
                    wScoreTotal = min(wScoreTotal,-5)
                    wString = "\tWat%d @ %s [%d,%d,%d,%d,%s,%s,%d,%d,%d] -> %d\n"%(m.residueNumber,nearestResidue,m.phobicRing,m.saltBridge,m.hydrophobic,m.otherContact,m.bulkWater,m.replacesOxygen.name,m.hbWat,m.hbAA,m.hbScore,wScoreTotal)
                    alignStr+=wString
                    missingEnergy += wScoreTotal#m.score + oRepScore
            else:
                missing = []
                missingEnergy = 0
            
            results.append((len(missing),"#MissWat"))
            results.append((missingEnergy,"EMissWat"))
            self.alignInfo.put(alignStr)
            
            #SideChain Stuff
            sideChainStr = "["+str(poseNum)+"]\n"
            sideChainData = []
            sideChainData += findHydrophobicContact(self.CONSTANTS,pose,targetTree)
            sideChainData += findHydrophilicBulkExposed(pose,tpCenter,targetTree,BULK_FACE_VALUE)
            sideChainData += findSaltBridges(self.CONSTANTS,pose,targetTree,tRes,wTree,poseWatersAll)
            sideChainData += findPlanarAmide(self.CONSTANTS,pose,targetTree,tRes)
            sideChainData += findLysineDP(self.CONSTANTS,pose,targetTree)
            sideChainStr = self.stringifySidechainData(sideChainStr,sideChainData)
            self.sideChainQueue.put(sideChainStr)
            hydrophobicNumber = 0
            badPocket = 0
            bulkNumber = 0
            saltNumber = 0
            pamideNumber = 0
            lysineDPNumber = 0
            sideEnergyTotal = 0.0
            for datum in sideChainData:
                rowKind = datum[ColumnValue.KIND.value]
                rowEnergy = datum[ColumnValue.ENERGY.value]
                if (rowKind == Kind.HYDROPHOBIC):
                    hydrophobicNumber+=1
                if (rowKind == Kind.WRONG_POCKET):
                    badPocket+=1
                if (rowKind == Kind.BULK_EXPOSED):
                    bulkNumber+=1
                if (rowKind == Kind.SALT_BRIDGE):
                    saltNumber+=1
                if (rowKind == Kind.PLANAR_AMIDE):
                    pamideNumber+=1
                if (rowKind == Kind.LYSINE_DP):
                    lysineDPNumber+=1
                sideEnergyTotal+=rowEnergy
                #print(datum)
            results.append((hydrophobicNumber,"#HyPhobic"))
            results.append((bulkNumber,"#BulkFace"))
            results.append((badPocket,"#BadPocket"))
            results.append((saltNumber,"#Salt"))
            results.append((lysineDPNumber,"#LysDP"))
            results.append((pamideNumber,"#PAmide"))
            results.append((sideEnergyTotal,"ESide"))
            
            totalEnergy = missingEnergy + hbEnergyTotal + sideEnergyTotal
            results.append((totalEnergy,"Energy"))
            
            self.task_queue.task_done()
            self.result_queue.put((poseNum,results))
            #import time
            #time.sleep(5)
        return
    def stringifySidechainData(me,string,data):
        for r in data:
            row = ""
            for i in range(1,len(r)):
                col = r[i]
                row+=str(col)+"\t"
            string+=row[:-1]+"\n"
        return string
    def expectedAndNeighborhoodWaters(me,poseTree,poseWaters,LIG_DIST_CUTOFF = 2.2,NEIGH_CUTOFF = 6.0):
        relevantWaters = me.baseWaterKeeperFine.getWatersForPose(poseTree.myList,LIG_DIST_CUTOFF)
        neigh = me.baseWaterKeeperFine.getWatersForPose(poseWaters,2.3)
        #print(len(relevantWaters),len(neigh))
        return relevantWaters,neigh
        #return [],[]
def splitAtomsIntoChains(atomList):
    aDict = {}
    for a in atomList:
        if (a.chain not in aDict):
            aDict[a.chain] = [a]
        else:
            aDict[a.chain].append(a)
    for key in aDict:
        aDict[key] = PDBTools.Peptide(aDict[key])
    return aDict
def findHydrophobicContact(CONSTANTS,pose,target):
    def scoreContacts(con):
        
        strong = []
        weak = []
        phile = []
        tooClose = []
        for crn in con:
            c = con[crn]
            if (c[2] < DistanceLimit):
                tooClose.append(c)
            elif (c[0].residueType in PDBTools.Atom.hydrophobicTable):
                if (c[1] < StrongRequired):
                    #print("Weak ",c[0],c[1],c[2])
                    #weak+=1
                    weak.append(c)
                else:
                    #print("Strong ",c[0],c[1],c[2])
                    #strong+=1
                    strong.append(c)
            #elif (c[0].residueType in PDBTools.Atom.hydrophilicTable):
            elif (c[0].residueType in ["ARG","LYS","ASP","GLU","GLN","ASN"]):
                #print("Phile ",c[0],c[1],c[2])
                #phile+=1
                phile.append(c)
        for s in tooClose:
            element = (Kind.HYDROPHOBIC,str(s[0]),str(s[3]),"Hydrophobic_Too_Close",s[2],0)
            ret.append(element)
        for s in strong:
            element = (Kind.HYDROPHOBIC,str(s[0]),str(s[3]),"Hydrophobic_Strong",s[2],2*HydrodrophicValue)
            ret.append(element)
        soFar = len(strong) + len(tooClose)
        
        for w in weak:
            mIndex = min(sizeMTable,soFar)
            mult = MultiplierTable[mIndex]
            element = (Kind.HYDROPHOBIC,str(w[0]),str(w[3]),"Hydrophobic_Weak",w[2],mult*HydrodrophicValue)
            ret.append(element)
            soFar+=1

        for p in phile:
            element = (Kind.WRONG_POCKET,str(p[0]),str(p[3]),"WrongPocket",p[2],-1*HydrodrophicValue)
            ret.append(element)
        total = weak + strong
        
        #print(a,strong,weak,phile)
    DistanceLimit =  CONSTANTS.get("HydrophobicContactDistanceLimit",2.3)
    Radius =  CONSTANTS.get("HydrophobicContactRadius",4.5)
    StrongRequired =  CONSTANTS.get("HydrophobicContactStrongNumberRequired",4)
    HydrodrophicValue = CONSTANTS.get("hydrophobicEnergy",-100)
    MultiplierTable = [1,1,1,1]
    #MultiplierTable = [1,0.75,0.5,0.25]
    sizeMTable = len(MultiplierTable)-1
    #(Contact Atom, Number of atoms seen, lowest distance, closest atom)
    curRes = None
    ret = []
    contacts = {}
    for a in pose.atoms:
        if (curRes is None) or (a.residueNumber > curRes):
            curRes = a.residueNumber
            scoreContacts(contacts)
            contacts = {}
        if (a.isHydrophobicAtom()):
            within = target.radiusSearch(a.location,Radius,returnWithDist=True)
            for conA,dist in within:
                if conA.isHydrophobicAtom() or (conA.isHydrophilicAtom() and conA.residueType!="TYR"):
                    otherResNum = conA.residueNumber
                    if otherResNum in contacts:
                        contacts[otherResNum][1]+=1
                        if (dist < contacts[otherResNum][2]):
                            contacts[otherResNum][2] = dist
                            contacts[otherResNum][3] = a
                    else:
                        contacts[otherResNum]=[conA,1,dist,a]
                        
    scoreContacts(contacts)
    return ret
def findHydrophilicBulkExposedOld(pose,center,target,value):
    #return []
    RADIUS = 4.5
    contacts = []
    lastUsed = -1
    for atm in pose.atoms:
        if (atm.residueNumber > lastUsed):
            lastUsed = atm.residueNumber
            key = atm.getHydrophilicAtom()
        if (key == atm.atomType):
            within = target.radiusSearch(atm.location,RADIUS,returnWithDist=True)
            valid = True
            nearbyPhiles = []
            #print(atm)
            for neigh,dist in within:
                #print("\t",neigh,dist)
                if (neigh.isHydrophilicRes()):
                    if (neigh.isBackbone()):
                        valid = False
                        break
                    if (neigh.residueNumber not in nearbyPhiles):
                        nearbyPhiles.append(neigh.residueNumber)
                else:
                    valid = False
                    break
            if (len(nearbyPhiles) > 2):
                valid = False
            if (valid):
                element = (Kind.BULK_EXPOSED,"",atm.identifier(),"Hydrophilic Bulk Exposed",str(len(nearbyPhiles)),value)
                contacts.append(element)
    return contacts  
def findHydrophilicBulkExposed(pose,center,tree,value):
    def ignoreConditions(atm):
        if atm.residueType in ["O","C","CA","N","H","HN"]:
            return False
        if atm.residueType in PDBTools.Atom.hydrophobicTable:
            return False
        return True
    #return []
    #RADIUS = 4.5
    AngleCutoff = 105
    NeighborRadius = 7.0
    contacts = []
    lastUsed = -1
    for atm in pose.atoms:
        if (atm.residueNumber > lastUsed):
            lastUsed = atm.residueNumber
            key = atm.getHydrophilicAtom()
        if (key == atm.atomType):
            neighbors = tree.radiusSearch(atm.location,NeighborRadius)
            good = True
            for neigh in neighbors:
                ang = center.angleThreePoints(atm.location,neigh.location)
                #print(ang)
                if (ang > AngleCutoff) and (not ignoreConditions(neigh)):
                    good = False
                    break          
            #if BaseWaters.isBulkExposed(target,atm,center=center):
            if good:
                element = (Kind.BULK_EXPOSED,"",atm.identifier(),"Hydrophilic Bulk Exposed","",value)
                contacts.append(element)
    return contacts  
def findSaltBridges(CONSTANTS,pose,target,tRes,wTree,poseWaters):
    def getAtomKey(atm,key):
        if (atm.chain == "L"):
            if pose.getSpecificAtom(atm.residueNumber,key) is None:
                raise Exception("Key %s not found in residue %i"%(key,atm.residueNumber))
            return pose.getSpecificAtom(atm.residueNumber,key)
        else:
            return tRes[atm.chain][atm.residueNumber][key]
    def getSaltBridgeComplements(atm):
        if (atm.isCTerminalCarbon()):
            ret = ["O","OXT"]
        else:
            ret = []
        ret+=BaseWaters.SB_CHARGED_ATOM.get(atm.residueType,[])
        return ret
    def getCompliments(atm):
        sbComplementNames = getSaltBridgeComplements(atm)
        sbComplements = []
        for s in sbComplementNames:
            sbComplements.append(getAtomKey(atm,s))
        return sbComplements
    def getBridgeWaters(cAtm):
        atmid = cAtm.identifier()
        entry = bridgeWaters.get(atmid,None)
        if entry is None:
            entry = wTree.radiusSearch(cAtm.location,bwRADIUS,returnWithDist=True)
            bridgeWaters[atmid] = entry
        return entry
    #def bondExistsBetweenChargedAtoms(ligChargedList,tpChargedList):
    def seekNonDirectSaltBridge(lAtom,tpAtom):
        lc = getCompliments(atm)
        tpc = getCompliments(otherAtom)
        for ligCharged in lc:
            for tpCharged in tpc:
                chargeDist = ligCharged.distance(tpCharged)
                #print(ligCharged,tpCharged,chargeDist)
                if atm.isCTerminalCarbon() and (chargeDist < ctRADIUS):
                    cTAng = atm.angle(ligCharged,tpCharged)
                    if (cTAng > angCUTOFF):
                        element = (Kind.SALT_BRIDGE,otherAtom.resStr(),atm.resStr(),"Salt Bridge: C-Terminal","%1.2f,%2.1f"%(chargeDist,cTAng),value)
                        contacts.append(element)
                        return
                else:
                    waters = getBridgeWaters(ligCharged)
                    #print("\tThere are %i waters"%len(waters)) 
                    for w,wDist in waters:
                        tpDist = w.distance(tpCharged)
                        #print("\t",w,wDist,tpDist)
                        if (tpDist < bwRADIUS):
                            element = (Kind.SALT_BRIDGE,"%s,%s"%(otherAtom.resStr(),w.resStr()),atm.resStr(),"Salt Bridge: SingleBridge",chargeDist,value)
                            contacts.append(element)
                            #print("\t\tWater mediated salt bridge")
                            return
        
    value = CONSTANTS.get("saltBridgeEnergy",-400)
    #return []
    cTerminal = pose.getCTerminalCarbon()
    wRADIUS = 7.0 * 1.15
    dRADIUS = 4.0 * 1.1
    ctRADIUS = 6.0
    bwRADIUS = 3.8
    angCUTOFF = 135
    contacts = []
    #bridgeMatchings = []
    bridgeWaters = {}
    #doubleBridgeWaters = {}
    for atm in pose.atoms:
        lType = atm.acidBase()
        #print(atm,lType)
        if (lType < "N"):
            if atm.isHydrophilicAtom():
                if (lType == "A"):
                    keyType = "B"
                else:
                    keyType = "A"
                within = target.radiusSearch(atm.location,wRADIUS,returnWithDist=True)
                for otherAtom,dist in within:
                    if otherAtom.isHydrophilicAtom():
                        if (otherAtom.acidBase() == keyType):
                            if (dist <= dRADIUS):
                                element = (Kind.SALT_BRIDGE,otherAtom.resStr(),atm.resStr(),"Salt Bridge: Direct",dist,value)
                                #print("\tDirect salt bridge")
                                contacts.append(element)
                            else:
                                seekNonDirectSaltBridge(atm,otherAtom)
                        elif ((otherAtom.acidBase() == lType) and (dist <= dRADIUS)):
                            element = (Kind.SALT_BRIDGE,otherAtom.resStr(),atm.resStr(),"Salt Bridge: Repulsive",dist,-value)
                            #print("\tDirect salt bridge")
                            contacts.append(element)                            
                            
    return contacts
def findPlanarAmide(CONSTANTS,pose,target,tRes):
    def mapKeys(source):
        atom = source["CA"]
        resData = amideData[atom.residueType]
        ret = {}
        for key in resData:
            ret[key] = source[resData[key]]
        return ret
    value = CONSTANTS['planarAmideEnergy']
    contactLength = 3.2
    searchLength = 4.4
    contacts = []
    amideData = {"ASN" : {"D" : "ND2", "A": "OD1", "C": "CG", "S" : "CB"}, "GLN": {"D" : "NE2", "A": "OE1", "C": "CD", "S" : "CG"}}
    for res in pose.residues:
        CA = res["CA"]
        if CA.residueType in amideData:
            ligData = mapKeys(res)
            stem = ligData["C"]
            seenResidues = set()
            #print(stem)
            within = target.radiusSearch(stem.location,searchLength,True, condition = lambda x : x.residueType in amideData)
            for neigh,cDist in within:
                seenKey = (neigh.chain,neigh.residueNumber)
                if (seenKey == (CA.chain,CA.residueNumber)):
                    continue
                if seenKey not in seenResidues:
                    #print("\t",neigh,seenKey)
                    seenResidues.add(seenKey)
                    try:
                        tpData = mapKeys(tRes[neigh.chain][neigh.residueNumber])
                    except:
                        #print("\tIncomplete")
                        continue

                    tpV = tpData["S"].location.difference(tpData["C"].location).unitVector()
                    lC_tS_V = ligData["C"].location.difference(tpData["S"].location).unitVector()
                    #print("\t\tVec: ",ligV.dotProduct(tpV))
                    if (lC_tS_V .dotProduct(tpV) < -0.7):
                        ligV = ligData["S"].location.difference(ligData["C"].location).unitVector()
                        if (ligV.dotProduct(tpV) < -0.8):
                            ligNormal = ligData["S"].location.planeNormal(ligData["A"].location,ligData["D"].location)
                            tpNormal = tpData["S"].location.planeNormal(tpData["A"].location,tpData["D"].location)
                            #print("\t\tPlane: ",abs(ligNormal.dotProduct(tpNormal)))
                            if (abs(ligNormal.dotProduct(tpNormal)) >= 0.5):
                                distA = ligData["D"].distance(tpData["A"])
                                if (distA <= contactLength):
                                    distB = tpData["D"].distance(ligData["A"])
                                    if (distB <= contactLength):
                                        #print(distA,distB)
                                        contacts.append((Kind.PLANAR_AMIDE,tpData["C"].resStr(),ligData["C"].resStr(),"Planar Amide",cDist,value))
    return contacts
    
def findLysineDP(CONSTANTS,pose,target):
    value = CONSTANTS['lysineDPEnergy']
    contactLength = 3.2
    contacts = []
    for res in pose.residues:
        CA = res["CA"]
        if (CA.residueType == "LYS"):
            NZ = res["NZ"]
            #print(NZ)
            within = target.radiusSearch(NZ.location,contactLength,True, condition = lambda x : x.atomType[0] == "O")
            if (len(within) >= 2):
                if (within[0][0].residueNumber == within[1][0].residueNumber):
                    del within[1]
            if (len(within) >= 2):
                dist = max(within, key =  lambda x : x[1])[1]
                contacts.append((Kind.LYSINE_DP,"%s & %s"%(within[0][0].resStr(),within[1][0].resStr()),NZ.resStr(),"Lysine DP",dist,value))
                #print("\tHas DP")
    return contacts
def oxygenReplacement(waters,poseTree,ignore):
    for w in waters:
        #print(w)
        w.replacesOxygen = OxyReplaceStatus.No
        if (not ignore):
            if w.residueType == "WAT":
                neigh = poseTree.radiusSearch(w.location,1.0)
                for n in neigh:
                    # print("\t",n,"\t",n.atomType[0] == "O")
                    if (n.atomType[0] == "O"):
                        if (n.atomType in ["O","OXT"]):
                            w.replacesOxygen = OxyReplaceStatus.Back
                            break
                        w.replacesOxygen = OxyReplaceStatus.Side
            
def waterBox(baseWaters,pose,boxSize):
    if (len(pose) < 1):
        return []
    #Finding highest and lowest value in each dimension
    firstAtomLoc = pose[0].location.dims
    minV = [firstAtomLoc[0],firstAtomLoc[1],firstAtomLoc[2]]
    maxV = [firstAtomLoc[0],firstAtomLoc[1],firstAtomLoc[2]]
    for atm in pose:
        ary = atm.location.dims
        for d in range(0,3):
            if (ary[d] < minV[d]):
                minV[d] = ary[d]
            if (ary[d] > maxV[d]):
                maxV[d] = ary[d]
    for d in range(0,3):
        minV[d]-=boxSize
        maxV[d]+=boxSize
    return watersInBox(baseWaters,minV,maxV)
def watersInBox(baseWaters,minV,maxV):
    #Creating a 'box' of waters that are close enough to the peptide to measure distance to nearest atom
    boxWaters = []
    for w in baseWaters:
        good = True
        ary = w.location.dims
        for d in range(0,3):
            between = (minV[d] <ary[d]) and (ary[d] < maxV[d])
            #print("\t",minV[d],ary[d],maxV[d])
            good = good and between
        if (good):
            boxWaters.append(w)
    return boxWaters
def expectedWaters(baseWaters,poseTree,LIG_DIST_CUTOFF = 2.2):
    #Calculating the distance between each water and it's nearest neighbor on the peptide
    #Then retaining that water if close enough for consideration
    ret = []
    #print("\nPose ",pose.indexNum)
    for wat in baseWaters:
        pNeigh,pDist = poseTree.nearestNeighborWithDist(wat.location)
        #print("Wat"+str(wat.residueNumber),pDist)
        if (pDist <= LIG_DIST_CUTOFF):
            ret.append(wat)
    return ret
        

def missingWaters(expected,align):
    if (len(align) == 0):
        return expected
    expectedSet = set([e.residueNumber for e in expected if e.residueNumber != 0])
    alignedSet = set([a for a in align if a != 0])
    missingNumbers = expectedSet.difference(alignedSet)
    return list(filter(lambda e: e.residueNumber in missingNumbers, expected))
def alignWater(baseWaterTree,pose,poseWaters, cutoff = 2.2, sumCutoff = 3.5):
    def addEdge(graph,left,right,dist = None):
        Xstr = "X"+str(left.residueNumber)
        graph.addNode(Xstr)            
        if dist is None:
            dist = left.distance(right)
        Ystr = "Y"+str(right.residueNumber)
        graph.addNode(Ystr)
        graph.addEdge(Xstr,Ystr,dist)
    baseWaters = baseWaterTree.myList
    align = [0] * len(poseWaters)
    poseSeen = [[] for  i in range(len(poseWaters))]
    baseSeen = [0] * len(baseWaters)
    for pI in range(len(poseWaters)):
        pw = poseWaters[pI]
        neighbors = baseWaterTree.radiusSearch(pw.location,cutoff,True)
        for neigh, neighDist in neighbors:
            bI = neigh.residueNumber
            poseSeen[pI].append((bI,neighDist))
            baseSeen[bI-1]+=1 
    graph = Network.Network()
    for i in range(len(poseSeen)):
        X = poseWaters[i]
        if (len(poseSeen[i]) == 1):
            mateIndex = poseSeen[i][0][0]
            mateSeen = baseSeen[mateIndex - 1]   
            if (mateSeen == 1):
                align[i] = mateIndex
                continue
        for j in range(len(poseSeen[i])):
            BWIndex, d = poseSeen[i][j]
            Y = baseWaters[BWIndex-1]
            addEdge(graph,X,Y,d)
    allSubgraphs = graph.splitIntoSubgraphs()
    for i in range(len(allSubgraphs)):
        subgraph = allSubgraphs[i]
        subgraph.cutoff = sumCutoff
        subgraph.addNode("S")
        subgraph.addNode("T")
        for subNodeID in subgraph.nodes:
            if (subNodeID[0] == "X"):
                subgraph.addEdge("S",subNodeID,0)
            elif (subNodeID[0] == "Y"):
                subgraph.addEdge(subNodeID,"T",0)
        result = subgraph.leastCostMatching(verbose=False)
        for match in result:
            xn = int(match[0][1:])
            yn = int(match[1][1:])
            align[xn-1] = yn
    return align
def getTargetDict(targetByChain,addH = False):
    grand = {}
    for chain in targetByChain:
        tRes = {}
        target = targetByChain[chain]
        firstRes = target.atoms[0].residueNumber
        for a in target.atoms:
            if a.residueNumber not in tRes:
                tRes[a.residueNumber] = {}
            tRes[a.residueNumber][a.atomType] = a
        if (addH):
            for resNum in tRes:
                if "H" not in tRes[resNum]:
                    if (resNum > firstRes):
                        try:
                            N = tRes[resNum]["N"]
                        except:
                            if (("C" not in tRes[resNum]) and ("CA" not in tRes[resNum])):
                                print(resNum," may be a cofactor: ",tRes[resNum].keys())
                                continue
                            else:
                                print(resNum)
                                print(tRes[resNum])
                        CA = tRes[resNum]["CA"]
                        try:
                            C = tRes[resNum-1]["C"]
                            tRes[resNum]["H"] = PDBTools.Atom.getH(N,CA,C)
                        except:
                            C = tRes[resNum]["C"]
                            tRes[resNum]["H"] = PDBTools.Atom.getH(N,CA,C,True)
        grand[chain] = tRes
                    
    return grand
def textifyHBondTable(table):
    if (table == None):
        return ""
    hbondText = ""
    for entry in table:
        hbondText+=str(entry)+"\n"
    return hbondText
def between(a,b,c):
    if (a < c):
        if (a > b):
            return False
        if (b > c):
            return False
        return True
    elif (a > c):
        if (a < b):
            return False
        if (b < c):
            return False
        return True
    elif (a == c):
        if (a == b):
            return True
    raise Exception("Strange error in between")
