import Vector
import Fixvar
import PDBTools
import os
#import shutil
import math 
recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
def dateToStamp(dateString):
    import datetime
    import time
    return time.mktime(datetime.datetime.strptime(dateString, "%d/%m/%Y").timetuple())
#updateCutoffDate = '06/05/2017'
updateCutoffDate = '06/10/2017'
updateCutoffStamp = dateToStamp(updateCutoffDate)

ACCEPTOR_DATA = {
    #"ARG":["NH1","NH2","NE"],
    "THR":["OG1"],
    "SER":["OG"],
    #"HIS":["NE2","ND1"],
    #"LYS":["NZ"],
    "ASP":["OD1","OD2"],
    "GLU":["OE1","OE2"],
    "ASN":["OD1","ND2"],
    "GLN":["OE1","NE2"],
    "TYR":["OH"]
    #"TRP":["NE1"]
    }
ACCEPTOR_STEM = {
    #"ARG":["NH1","NH2","NE"],
    "THR":{"OG1":"CB"},
    "SER":{"OG":"CB"},
    #"HIS":["NE2","ND1"],
    #"LYS":["NZ"],
    "ASP":{"OD1":"CG","OD2":"CG"},
    "GLU":{"OE1":"CD","OE2":"CD"},
    "ASN":{"OD1":"CG","ND2":"CG"},
    "GLN":{"OE1":"CD","NE2":"CD"},
    "TYR":{"OH":"CZ"}
    #"TRP":["NE1"]
    }
PLANE_STEM = {
    "THR":{"OG1":None},
    "SER":{"OG":None},
    "ASP":{"OD1":"CB","OD2":"CB"},
    "GLU":{"OE1":"CG","OE2":"CG"},
    "ASN":{"OD1":"CB","ND2":"CB"},
    "GLN":{"OE1":"CG","NE2":"CG"},
    "TYR":{"OH":"CE1"}
    }
    
DONOR_DATA = {
    "ARG":["1HH1","1HH2","2HH1","2HH2","HNE"],
    "THR":["HOG"],
    "SER":["HOG"],
    "HIS":["HND"],
    "LYS":["HZ1","HZ2","HZ3"],
    "ASN":["HND1","HND2"],
    "GLN":["HE21","HE22","1HE2","2HE2"],
    "TYR":["HOH"],
    "TRP":["HNE"]
    }
            
DONOR_STEM = {
    "ARG":{"1HH1":"NH1","2HH1":"NH1","1HH2":"NH2","2HH2":"NH2","HNE":"NE"},
    "THR":{"HOG":"OG1"},
    "SER":{"HOG":"OG"},
    "HIS":{"HND":"ND1"},
    "LYS":{"HZ1":"NZ","HZ2":"NZ","HZ3":"NZ"},
    "ASN":{"HND1":"ND2","HND2":"ND2"},
    "GLN":{"HE21":"NE2","HE22":"NE2","1HE2":"NE2","2HE2":"NE2"},
    "TYR":{"HOH":"OH"},
    "TRP":{"HNE":"NE1"}
    }
SB_CHARGED_ATOM = {
    "ASP":["OD1","OD2"],
    "GLU":["OE1","OE2"],
    "LYS":["NZ"],
    "ARG":["NH1","NH2","NE"]
    }
    
def getBaseWaterFilename(CONSTANTS):
    return "%s/%sbase.pdb"%(CONSTANTS["commonFolder"],CONSTANTS["tpName"])
def shouldGenerateBaseWaterFile(CONSTANTS):
    commonFolder = CONSTANTS["commonFolder"]
    fname = getBaseWaterFilename(CONSTANTS)
    if os.path.isfile(fname):
        creation = os.stat(fname).st_mtime
        shouldRegen = updateCutoffStamp > creation
        #input("Created:%i\nCutoff%i\nUpdate?%s"%(creation,updateCutoffStamp,str(shouldRegen)))
        if (shouldRegen):
            cutoffString = updateCutoffDate.replace("/","_")
            #input(cutoffString)
            os.rename(fname,"%s/%sbase_%s.pdb"%(CONSTANTS["commonFolder"],CONSTANTS["tpName"],cutoffString))
            PRINT("\nOutdate Base Water File\n")
        return shouldRegen
    return True
def getBaseWaters(CONSTANTS):
    fname = getBaseWaterFilename(CONSTANTS)
    if os.path.isfile(fname):
        return PDBTools.readAtomsInChain(fname,"W","O")
    else:
        return None
    return None
def getBaseWatersAllAtoms(CONSTANTS):
    commonFolder = CONSTANTS["commonFolder"]
    ignoreGenWaters = CONSTANTS.get("ignoreGenWaters",False)
    fname = commonFolder+"/"+CONSTANTS["tpName"]+"base.pdb"
    ret = []
    if os.path.isfile(fname):
        atoms = PDBTools.readAtomsInChain(fname,"W")
        if ignoreGenWaters:
            atoms = list(filter(lambda a:a.residueType != "GEN",atoms))
        wat = []
        for i in range(0,len(atoms),3):
            wat.append(tuple(atoms[i:i+3]))
        return tuple(wat)
    else:
        return None
def isAcceptor(a):
    if (a.atomType in ["O","OXT"]):
        return True
    if (a.residueType in ACCEPTOR_DATA):
        if (a.atomType in ACCEPTOR_DATA[a.residueType]):
            return True
    return False
def isDonor(a):
    if (a.atomType[:2] == "HN"):
        return True
    if (a.residueType in DONOR_DATA):
        if (a.atomType in DONOR_DATA[a.residueType]):
            return True
    return False

def getHBDandHBA(CONSTANTS):
    commonFolder = CONSTANTS["commonFolder"]
    fname = getBaseWaterFilename(CONSTANTS)
    chainList = []
    for c in CONSTANTS["targetProteinChain"]:
        chainList.append(c)
    #print("Chains: ",chainList)

    atoms = []
    for c in chainList:
        atoms = atoms + PDBTools.readAtomsInChain(fname,c)
        
    acceptor = []
    hydrogen = []
    for a in atoms:
        if isAcceptor(a):
            acceptor.append(a)
        if isDonor(a):
            hydrogen.append(a)
    return acceptor, hydrogen
def getCenterAlphaCarbon(ot,resNum,radius):
    base = None
    for CA in ot.CAList:
        if (CA.residueNumber == resNum):
            base = CA
            break
    if (base == None):
        raise Exception("ERROR. Target Residue Not Found")
    return base
def reduceBaseWaters(wats,ot,resNum,radius):
    close = []
    PRINT(str(len(wats))+" waters presented")
    base = getCenterAlphaCarbon(ot,resNum,radius)
    for w in wats:
        if (w.distance(base) < radius):
            close.append(w)
    PRINT(str(len(close))+" waters in radius")
    return close

def residueMap(target):
    if (len(target.atoms) < 1):
        return {}
    curRes = target.atoms[0].residueNumber
    res = [curRes]
    curChain = target.atoms[0].chain
    ret = {}
    ret[curChain] = res
    for a in target.atoms:
        if (a.chain != curChain):
            curRes = a.residueNumber
            res = [curRes]
            curChain = a.chain
            ret[curChain] = res
            
        if (curRes < a.residueNumber):
            res.append(a.residueNumber)
            curRes = a.residueNumber
            
    return ret
def makeBaseWaters(CONSTANTS,tree,target):
    import ZMAT
    import subprocess
    import WaterLayer
    rmap = residueMap(target)
    commonFolder = CONSTANTS["commonFolder"]
    #Find a point far from any atoms in the target
    GRAIN = 10
    HIGH = 99
    LOW = 12
    LIMIT = 20  
    center = Vector.Vector3D(LOW,LOW,LOW)
    neigh,dist = tree.nearestNeighborWithDist(center)
    while (dist < LIMIT):
        center.dims[0]+=GRAIN
        if(center.dims[0] > HIGH):
            center.dims[0] = LOW
            center.dims[1]+=GRAIN
            if(center.dims[1] > HIGH):
                center.dims[1] = LOW
                center.dims[2]+=GRAIN
        neigh,dist = tree.nearestNeighborWithDist(center)
    PRINT(str(center)+"\t\t"+str(dist))

    #Load ZMAT with shortest zmat file in common folder
    zmatChoices = []
    for file in os.listdir(commonFolder):
        if file.endswith(".zmat"):
            zmatChoices.append(file)
    zmatChoices = sorted(zmatChoices, key=lambda x: len(x))
    smallest = zmatChoices[0]

 
    locationZmat = commonFolder+"/"+smallest
    #print(locationZmat)
    zmat = ZMAT.ZMAT(locationZmat)


    #Make reference atoms for ZMAT object, based on far point calculated above
    distAtm = PDBTools.Atom.blankAtom()
    distAtm.location = center
    angLoc = center.translate([1,0,0])
    angAtm = PDBTools.Atom.blankAtom()
    angAtm.location = angLoc
    torLoc = center.translate([1,1,0])
    torAtm = PDBTools.Atom.blankAtom()
    torAtm.location = torLoc
    zmat.setReference(distAtm,angAtm,torAtm)


    #Make fake fixvar with no variance, then use to make a dummy peptide far from target protein
    blankFxvr = Fixvar.Fixvar(None)
    blankFxvr.kind.append(0)
    blankFxvr.atom.append(4)
    nullPose = Fixvar.FixvarPoses([0],1)
    blankFxvr.poses.append(nullPose)
    pep = zmat.makePDBString(blankFxvr,0,addTer=False)
    #print("Calculated peptide atoms\n",pep)

    #Processing Target for Watgen
    #deHTarget = target.createNoHCopy()
    PRINT("Creating Base Water Structure using dummy peptide centered at "+str(center))
#CONSTANTS["targetProteinChain"]

    p = subprocess.Popen(["java","-jar",CONSTANTS["universalFolder"]+"/"+CONSTANTS["WGProbe"],"A","6","2"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #PRINT(str(currentPoseData.firstIndex()))
    #PRINT("First "+str(currentPoseData.firstIndex()))
    #PRINT(currentPoseData.getPoseStringArray(currentPoseData.firstIndex()))
    ball = str(target)+pep+"__DONE__\r\n"
    '''
    ballDump = open("BallDump.out",'w')
    ballDump.write(ball)
    ballDump.close()
    '''
    #PRINT(ball)
    #print(ball)

    backTuple = p.communicate(ball.encode())
    #print(backTuple)
    back = backTuple[0].decode('utf-8')
    backError = backTuple[1].decode('utf-8')
    if (len(backError) > 1):
        PRINT("Java Gave Error:\n"+backError)
    back = back.replace("\r","")
    #PRINT(">"+str(back))
    '''
    ballDump = open("BackDump.out",'w')
    ballDump.write(back)
    ballDump.close()
    '''
    backAry = back.split("\n")
    '''
    for key in rmap:
        print(key)
        cChain = rmap[key] 
        for i in range(len(cChain)):
            res = cChain[i]
            print("\t",i,"->",res)
    '''
    
    watGenFileName = "%s/%sWatgenBase.pdb"%(commonFolder,CONSTANTS["tpName"])
    outf = open(watGenFileName,'w')
    watgenValid = False
    receivedAtoms = False
    for line in backAry:
        if (line[:4] == "ATOM"):
            a = PDBTools.Atom(line)
            if (a.valid):
                receivedAtoms = True
            if (a.chain in CONSTANTS["targetProteinChain"]):
                #print(a.toPDBLine(),"  ->   ",a.chain,",",a.residueNumber)
                a.residueNumber = rmap[a.chain][a.residueNumber-1]
            if (a.residueType == "WAT"):
                watgenValid = True
            outf.write(a.toPDBLine()+"\n")          
    outf.close()
    if (not receivedAtoms):
        raise Exception("Watgen did not return any atoms")
    if (not watgenValid):
        raise Exception("Watgen did not generate any base waters")
    PRINT("Watgen Base Waters Created. Adding Deeper Layers")
    finalFileName = "%s/%sBase.pdb"%(commonFolder,CONSTANTS["tpName"])
    layerRoutine = WaterLayer.WaterLayer(watGenFileName,CONSTANTS["targetProteinChain"])
    layerRoutine.generateWaterLayer()
    layerRoutine.outputFile(finalFileName)
def getPlaneStem(residue,atom):
    if (atom == "O"):
        return "CA"
    return PLANE_STEM[residue][atom]
def getAcceptorStem(residue,atom):
    if (atom == "O"):
        return "C"
    return ACCEPTOR_STEM[residue][atom]
def getDonorStem(residue,atom):
    if (atom[:2] == "HN"):
        if (len(atom) == 2) or (atom[2] in ["1","2","3"]):
            return "N"
    if (atom == "H"):
        return "N"
    return DONOR_STEM[residue][atom]
def setupBaseWaters(CONSTANTS):
    import PoseScorer
    import KDTree
    commonFolder = CONSTANTS["commonFolder"]
    targetPDBFile = "%s/%s"%(commonFolder,CONSTANTS["targetProtein"])
    baseWatersWithH = getBaseWatersAllAtoms(CONSTANTS)
    if baseWatersWithH is None:
        raise Exception("No base water file")
    bwH = []
    baseWaters = []
    for w in baseWatersWithH:
        for i in range(len(w)):
            if (i % 3 == 0):
                baseWaters.append(w[i])
            else:
                bwH.append(w[i])

    waterTree = KDTree.KDTree.loadAtomArray(baseWaters)
    
    #target = PDBTools.readPDBSelectively(targetPDBFile,self.CONSTANTS["targetProteinChain"])
    targetByAtoms = list(filter(lambda x : not x.isHydrogen(),PDBTools.readPDBSelectively(targetPDBFile,CONSTANTS["targetProteinChain"])))
    targetTreeNoH = KDTree.KDTree.loadAtomArray(targetByAtoms)
    targetWithH = PDBTools.readPDBSelectively(getBaseWaterFilename(CONSTANTS),CONSTANTS["targetProteinChain"])
    targetByChain = PoseScorer.splitAtomsIntoChains(targetWithH)    
    targetTreeWithH = KDTree.KDTree.loadAtomArray(targetWithH)
    tarOxygen, tarHydrogen = getHBDandHBA(CONSTANTS)
    possibleDonor = tarHydrogen + bwH
    possibleAcceptor = tarOxygen + baseWaters
    hTree = KDTree.KDTree.loadAtomArray(possibleDonor)
    oTree = KDTree.KDTree.loadAtomArray(possibleAcceptor)
    tRes = PoseScorer.getTargetDict(targetByChain,True)
    bwDumpFile = open(CONSTANTS["commonFolder"]+"/baseWaterDump.txt",'w')
    
    classifyDepth(waterTree,targetTreeNoH)
    countWaterHB(CONSTANTS,baseWatersWithH,targetTreeWithH,oTree,hTree,tRes,dumpfile=bwDumpFile)
    return baseWatersWithH
def countWaterHB(CONSTANTS,waters,targetTree,oTree,hTree,tRes,dumpfile = None):
    ANG_CUTOFF = 140
    TARGET_DIST_CUTOFF = 3.5
    RADIUS = 2.5
    SB_ANG_CUTOFF = 120
    def exclude(aFrom,contacts,this): 
        mid = [c for c in contacts if (c.residueNumber != this)]
        ang = [d for d in mid if (getAngle(aFrom,d) >= ANG_CUTOFF)]
        if(dumpfile is not None):
            for hb in ang:
                dumpfile.write("\tHydrogen bond: %s\n"%(hb))
        return ang
    def classifyHB(allHB):
        wat = 0
        aa = 0
        for hb in allHB:
            if hb.residueType in ["WAT","HOH"]:
                wat+=1
            else:
                aa+=1
        return wat,aa
    def getAngle(aBase,aOther):      
        if (aBase.atomType[:1] == "H"):
            baseNum = aBase.residueNumber
            angO = waters[baseNum-1][0]
            ang = angO.angle(aBase,aOther)
            return ang
        elif (aBase.atomType[:1] == "O"):
            otherNum = aOther.residueNumber
            if (aOther.chain == "W"):
                otherW = waters[otherNum-1]
                angO = otherW[0]
                ang = angO.angle(aOther,aBase)
                return ang
            else:
                try:
                    otherRes = tRes[aOther.chain][aOther.residueNumber]
                    seek = getDonorStem(aOther.residueType,aOther.atomType)
                    stem = otherRes[seek]
                    ang = stem.angle(aOther,aBase)
                except:
                    print(aOther)
                    print("\t",aBase)
                    print("\tOther atom is from target")
                
 
                    print("[",end='')
                    for o in otherRes:
                        print(o+",",end='')
                    print("]")
                    print("We want ",seek)
                    print("We got ",stem)   

                    print("\t",stem)
                
                    print(ang)

                return ang

    for w in waters:
        O = w[0]
        H1 = w[1]
        H2 = w[2]
        wNum = O.residueNumber
        O.phobicRing = 0
        O.hydrophobic = 0
        O.otherContact = 0
        O.saltBridge = 0

        '''
        #Grabs all TP atoms in 15A radius, then calculates the average point to use as an anchor
        tpVicinity = tree.radiusSearch(w.location,15.0)
        center = calcAveragePoint(tpVicinity)
        '''
        O.bulkWater = isBulkExposed(targetTree,O,NeighborRadius=5.0)
        if (O.bulkWater):
            dumpfile.write("\tBulkExposed")
        
        #input("Press enter")
        if(dumpfile is not None):
            dumpfile.write("\n%s\n"%(O))
            #print("\n\n%s"%(O))

        tNeigh = targetTree.radiusSearch(O.location,TARGET_DIST_CUTOFF,returnWithDist=True)
        neighDict = {}
        for i in range(len(tNeigh)):
            a = tNeigh[i][0]
            resNum = a.residueNumber
            if (resNum in neighDict):
                neighDict[resNum].append(i)
            else:
                neighDict[resNum] = [i]

        for resNum in neighDict:
            residue = neighDict[resNum]
            indexFirstAtom = residue[0]
            firstAtom = tNeigh[indexFirstAtom][0]
            resType = firstAtom.residueType
            SaltBridgeResidues = ("ASP","GLU","LYS","ARG")
            HydrophobicRingResidues = ("PHE","TYR","TRP")
            saltBridge = False
            if (resType in SaltBridgeResidues):
                #print("Check salt bridge to "+firstAtom.resStr())
                chargedAtoms = []
                sbKind = firstAtom.acidBase()
                for aIndx in residue:
                    neigh = tNeigh[aIndx][0]
                    if (sbKind == "A"):
                        if (neigh.atomType in ACCEPTOR_DATA[resType]):
                            chargedAtoms.append(neigh)
                            #print("\t%s"%neigh)
                    if (sbKind == "B"):
                        if (neigh.atomType in DONOR_DATA[resType]):
                            chargedAtoms.append(neigh)                       
                            #print("\t%s"%neigh)
                if (len(chargedAtoms) > 0):
                    for ca in chargedAtoms:
                        if (sbKind == "A"):
                            for wH in [H1,H2]:
                                dH = wH.distance(ca)
                                if (dH < RADIUS):
                                    ang = ca.angle(wH,O)
                                    #print("\tH in distance %f with angle %f"%(dH,ang))
                                    if (ang >= SB_ANG_CUTOFF):
                                        saltBridge = True
                                        #print("\tSalt Bridge: %s\n"%(ca.resStr()))
                                        break
                        if (sbKind == "B"):
                            dO = O.distance(ca)
                            if (dO < RADIUS):
                                seek = getDonorStem(resType,ca.atomType)
                                stem = tRes[ca.chain][resNum][seek]
                                ang = stem.angle(ca,O)
                                #print("\tO in distance %f with angle %f"%(dO,ang))
                                if (ang >= SB_ANG_CUTOFF):
                                    saltBridge = True
                                    #print("\tSalt Bridge: %s\n"%(ca.resStr()))
                                    break
                if (saltBridge):
                    if(dumpfile is not None):
                        dumpfile.write("\tSalt Bridge: %s\n"%(ca.resStr()))
                    O.saltBridge+=1
            ringContacts = 0
            phobicContacts = 0
            for aIndx in residue:
                neigh = tNeigh[aIndx][0]
                #print(neigh)
                #neighDist = tNeigh[aIndx][1]
                if (neigh.isHydrophobicAtom()):
                    phobicContacts+=1
                    if (resType in HydrophobicRingResidues):
                        ringContacts+=1
            '''
            if (resType in HydrophobicRingResidues):
                print("\t%d Ring contacts found from %s"%(ringContacts,firstAtom.resStr()))
            if (firstAtom.isHydrophobicRes()):
                print("\t%d Phobic contacts found from %s"%(phobicContacts,firstAtom.resStr()))
            '''
            if (ringContacts >=4):
                if(dumpfile is not None):
                    dumpfile.write("\tRing contact: %s\n"%(neigh.resStr()))
                O.phobicRing+=1
            elif (phobicContacts > 0):
                if(dumpfile is not None):
                    dumpfile.write("\tHydrophobic: %s\n"%(neigh.resStr()))
                O.hydrophobic+=1
            else:
                if(dumpfile is not None):
                    dumpfile.write("\tOther contact: %s\n"%(neigh.resStr()))
                O.otherContact+=1
                
        acc1 = oTree.radiusSearch(H1.location,RADIUS)
        acc2 = oTree.radiusSearch(H2.location,RADIUS)
        don = hTree.radiusSearch(O.location,RADIUS)
        acc2 = exclude(H2,acc2,wNum)
        acc1 = exclude(H1,acc1,wNum)
        don = exclude(O,don,wNum)
        O.hbWat, O.hbAA = classifyHB(don+acc1+acc2)
        '''
        O.don = len(don)
        O.acc = len(acc1) + len(acc2)
        O.countHB = O.don + O.acc
        '''
        #depthScore = O.bulkDepth * (O.bulkDepth + 1) * 5
        depthScore = 0
        O.hbScore = CONSTANTS["waterHBondWater"] * O.hbWat +  CONSTANTS["waterHBondAA"] * O.hbAA + CONSTANTS["waterSaltBridge"] * O.saltBridge
        O.score = CONSTANTS["waterPhobic"] * O.hydrophobic + CONSTANTS["waterOtherContact"] * O.otherContact + CONSTANTS["waterRingContact"] * O.phobicRing - depthScore + CONSTANTS["waterBulkExposure"] * O.bulkWater
        if(dumpfile is not None):
            if (O.wasStranded):
                dumpfile.write("\tStranded Water\n")
            dumpfile.write("\tDepth: %i Value: %i\n"%(O.bulkDepth,-depthScore))
            dumpfile.write("\tScore: %d\n"%(O.hbScore + O.score))
            
def isBulkExposed(tree,atom,center = None,AngleCutoff=105,NeighborRadius=7.0,CenterPointRadius=15.0):
    '''
    CenterPointRadius = 15.0
    NeighborRadius = 7.0
    AngleCutoff = 105
    '''
    if center is None:
        #Grabs all TP atoms in 15A radius, then calculates the average point to use as an anchor
        tpVicinity = tree.radiusSearch(atom.location,CenterPointRadius)
        center = calcAveragePoint(tpVicinity)
    #Grabs all TP atoms in 6A radius to see if any are "away" from center relative to the water, defined by more than 115 degrees
    neighbors = tree.radiusSearch(atom.location,NeighborRadius)
    for neigh in neighbors:
        ang = center.angleThreePoints(atom.location,neigh.location)
        #print(ang)
        if (ang > AngleCutoff):
            return False
    return True
def classifyDepth(waterTree,tree):
    import Network
    waters = waterTree.myList
    for w in waters:
        w.wasStranded = False
        w.bulkDepth = 0
    return

    #Graph Object, treats each water as a node
    net = Network.Network()
    #Because the waters are already fed in as a KDTree, this simply grabs the list of all waters from the water KDTree.
    #The "tree" variable is the Target Protein but doesn't need to be accesed in order
    waters = waterTree.myList
    #Adds an imaginary "start" node that all bulk waters connect to
    net.addNode("S")
    #Adds a node for each water
    for water in waters:
        net.addNode(water.resStr())
    for i in range(len(waters)):       
        water = waters[i]
        watID = water.resStr()
        isBulk = isBulkExposed(tree,water)
        #Bulk waters connect to start node
        if (isBulk):
            net.addEdge("S",watID,1)
        #Water nodes connect to all neighboring waters, that is, within 2.8A
        neighborWaters = waterTree.radiusSearch(water.location,2.8,condition = lambda x: x!=water)
        #If no other water is within 2.8A, water connects to nearest other water
        if (len(neighborWaters) == 0):
            neighborWaters = [waterTree.nearestNeighbor(water.location,condition = lambda x: x!=water)]
        #Adds edges between the neighboring waters
        for nw in neighborWaters:
            net.addMutualEdge(watID,nw.resStr(),1)
    #Calculates each water node's distance from Start node
    net.mapDepth("S")
    #Applies the depth data to the water object, instead of the Network node object. Sets any water that did not find any connection to bulk as "-1"
    maxDepth = 0
    stranded = []
    for water in waters:
        calcDepth = net.nodes[water.resStr()].depth
        if (calcDepth):
            water.bulkDepth =  calcDepth - 1
            water.wasStranded = False
            if (water.bulkDepth > maxDepth):
                maxDepth = water.bulkDepth
        else:
            stranded.append(water)
    for water in stranded:
        water.bulkDepth = maxDepth
        water.wasStranded = True
def calcAveragePoint(atoms):
    import Vector
    size = len(atoms)
    if (size == 0):
        return None
    dim = [0,0,0]
    for atm in atoms:
        for d in range(3):
            dim[d]+=atm.location.dims[d]
    for d in range(3):
        dim[d]/=size
    return Vector.Vector3D(dim)
class BaseWaterKeeper:
    def __init__(me,allWaters,dist):
        me.grain = dist
        me.allWaters = allWaters
        first = me.allWaters[0][0].location.dims
        me.low = [first[0],first[1],first[2]]
        me.high = [first[0],first[1],first[2]]
        for i in range(1,len(me.allWaters)):
            w = me.allWaters[i][0].location.dims
            for d in range(3):
                if (w[d] > me.high[d]):
                    me.high[d] = w[d]
                if (w[d] < me.low[d]):
                    me.low[d] = w[d]
        me.dRange = [None,None,None]
        for d in range(3):
            me.dRange[d] = math.ceil((me.high[d] - me.low[d]) / me.grain)
        '''
        print(me.low)
        print(me.dRange)
        print(me.high)
        '''
        me.master = {}
        for w in me.allWaters:
            O = w[0]
            loc = me.place(O)
            #print(O,loc)
            group = me.master.get(loc,None)
            if group is None:
                me.master[loc] = [O]
            else:
                group.append(O)

    def place(me,atom):
        loc = [None,None,None]
        for d in range(3):
            loc[d] = math.floor((atom.location.dims[d] - me.low[d]) / me.grain)
        return tuple(loc)
    def getWatersForPose(me,pose,dist):
        boxes = {}
        for a in pose:
            loc = me.place(a)
            for x in range(max(loc[0]-1,0),min(loc[0]+2,me.dRange[0])):
                for y in range(max(loc[1]-1,0),min(loc[1]+2,me.dRange[1])):
                    for z in range(max(loc[2]-1,0),min(loc[2]+2,me.dRange[2])):
                        pLoc = (x,y,z)
                        box = boxes.get(pLoc,None)
                        if box is None:
                            boxes[pLoc] = [a]
                        else:
                            box.append(a)
        ret = []
        for box in boxes:
            waters = me.master.get(box,[])
            for wat in waters:
                for atm in boxes[box]:
                    if (atm.distance(wat) < dist):
                        ret.append(wat)
                        break
            
        return ret
        
        
