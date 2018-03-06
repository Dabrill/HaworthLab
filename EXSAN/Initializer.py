import HBond
from ZMAT import zmatObj,zmatName
from KDTree import KDTree
from PoseScorer import splitAtomsIntoChains,getTargetDict
from BaseWaters import getHBDandHBA,getAcceptorStem,getPlaneStem
import math
import Counter
import multiprocessing
MinDist = 1.0
MaxDist = 2.8
OCbondLength = 1.23
CarbonylAtom = 4

CD = (30,37,1)
CA = (90,179,1)
CT = (0,359,1)
'''
CD = (30,31,1)
CA = (100,112,1)
CT = (70,90,1)
'''

CaA = (140,179,10)
CaT = (0,359,30)
OT = (0,179,30)
'''
CD = (30,37,1)
CA = (90,179,5)
CT = (0,359,5)
CaA = (140,179,10)
CaT = (0,359,30)
OT = (0,179,45)
'''
HB_IdealDistance = 1.9
HB_Scaling = 2
def hbDistEnergy(dist,base):
    e_by_dist = base*math.exp(-(HB_Scaling*(dist-HB_IdealDistance)**2))
    return e_by_dist
def initCTerminal(step,cons,targetAll):
    print("\n\n")
    from Fixvar import Fixvar,FixvarPoses
    getOutDotPDB(cons)
    
    cores = cons["cores"]
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    workers = [BestAngleFinder(step,cons,targetAll,tasks,results) for i in range(cores)]
    first = workers[0]
    
    cLocs = first.getCarbonylLocations()
    numberOfPoses = len(cLocs)
    countDisp = Counter.Counter(len(cLocs),"%d/%d   initial")
    count = 0
    newFxvr = Fixvar(None)
    newFxvr.atom = [4,4,4,5,5,8,6]
    newFxvr.kind = [2,1,0,1,0,0,0]
    headerZmatName = zmatName(step)
    headStr = "%8s %s"%(len(newFxvr.kind),headerZmatName)
    newFxvr.header = [headStr]
    print("There are %i cLocs"%len(cLocs))

    allUnique = set()
    BUFFER = 300
    fed = 0
    recieved = 0
    while (fed < BUFFER) and (fed < numberOfPoses):
        poseGeoFeed = cLocs[fed]
        tasks.put(poseGeoFeed)
        fed+=1
    for w in workers:
        w.prepareAtoms()
        w.start()
    #input("Extendibility Cull loaded and prepared")
    countDisp = Counter.Counter(numberOfPoses,"%d/%d CTerminal Initiation")
    count = 0
    backPoses = []
    while (recieved < numberOfPoses):
        if (tasks.qsize() < BUFFER):
            if (fed < numberOfPoses):
                poseGeoFeed = cLocs[fed]
                tasks.put(poseGeoFeed)
                fed+=1            
            elif (fed < numberOfPoses + cores):
                tasks.put(None)
                fed+=1
        if not results.empty():
            comm = results.get()
            if (comm is not None):
                backPoses+=comm
                count+=len(comm)
            recieved+=1
            countDisp.disp(recieved)
    print("Poses created: %i, now culling"%count)
    backPoses = sorted(backPoses,key = lambda x: x[1])
    finalPoses = []
    cutoff = step.get("energyCutoff",0.8)
    i = 0
    maxEnergy = backPoses[0][1]
    for i in range(len(backPoses)):
        pose = backPoses[i]
        frac = pose[1] / maxEnergy
        if (frac < cutoff):
            break
        finalPoses.append(pose[0])
    
    print("Final count: ",i)    
    for fp in finalPoses:
        newFxvr.poses.append(FixvarPoses(fp,0))
    newFxvr.sortByAngles()
    newFxvr.createApprovedFixvar("fixvar.out")
def getOutDotPDB(cons):
    import os
    import shutil
    import TMD
    top = os.getcwd()
    currentFolder = "%s/TMD"%top
    try:
        os.stat(currentFolder)
    except:
        os.mkdir(currentFolder)
    os.chdir(currentFolder)

    step = {"type":"TMD",
         "sequence":"G",
         "zsequence":"G#",
         "folder":"1-G",
         "zmatSuffix":"_aa1",
         "nclsrs":0,
         "numFixvar":0,
         "nDistanceCheck":0,
         "reverse":0,
         "variableTorsions":[],
         "criteria":[]}
    
    TMD.runTMD(cons,step,None)
    shutil.copy("out.pdb","%s/out.pdb"%top)
    os.chdir(top)
    shutil.rmtree(currentFolder)
class BestAngleFinder(multiprocessing.Process):
    def __init__(me,step,cons,targetAll,task_queue,result_queue):
        multiprocessing.Process.__init__(me)
        me.cons = cons
        me.step = step
        me.targetAll = targetAll
        me.targetAll.makeDictionary()
        me.targetTree = KDTree.loadAtomArray(list(filter(lambda x: not x.isHydrogen(),targetAll.atoms)))
        me.centerResNum = step["centerResidue"]
        me.referenceResNums = step["referenceResidues"]
        
        if (me.centerResNum not in me.referenceResNums):
            me.referenceResNums.append(me.centerResNum)
        for rNum in me.referenceResNums:
            #print(targetAll.getSpecificAtom(rNum,"N").toPDBLine())
            me.targetAll.implyHydrogen(rNum)
        targetByChain = splitAtomsIntoChains(me.targetAll.atoms)
        me.tRes = getTargetDict(targetByChain,True)
        acceptorList,donorList = getHBDandHBA(cons)
        me.donorTree = KDTree.loadAtomArray(donorList)
        me.acceptorTree = KDTree.loadAtomArray(acceptorList)

        me.tpChain = step.get("chain",None)
        if (me.tpChain is None):
            if (len(cons["targetProteinChain"]) == 1):
                me.tpChain = cons["targetProteinChain"][0]
            else:
                raise Exception("No chain given for initiator routine and none can be assumed")
        me.tpDist = targetAll.getSpecificAtom(me.centerResNum,"N")
        me.tpAng = targetAll.getSpecificAtom(me.centerResNum,"CA")
        me.tpTor = targetAll.getSpecificAtom(me.centerResNum,"C")

        me.zmatObj = zmatObj(step)


        me.task_queue = task_queue
        me.result_queue = result_queue
    def getCarbonylLocations(me):
        C = me.zmatObj.getAtomData(1,"C",True)
        good = []
        cMin = int(10*(OCbondLength + MinDist))
        cMax = int(10*(OCbondLength + MaxDist))
        ###XXX
        cMin = CD[0]
        cMax = CD[1]
        for intDist in range(cMin,cMax+1):
            dist = intDist / 10.0
            for ang in range(CA[0],CA[1]+1,CA[2]):
                for tor in range(CT[0],CT[1]+1,CT[2]):
                    geo = (intDist,ang,tor)
                    C.setAtomLocationByZMAT(me.tpDist,me.tpAng,me.tpTor,dist,ang,tor)
                    neigh, closestDist = me.targetTree.nearestNeighbor(C.location,True)
                    if (closestDist > me.cons["clashLimit"]):
                        good.append(geo)
        return good
    def prepareAtoms(me):
        me.C = me.zmatObj.getAtomData(1,"C",True)
        
        me.CAdata = me.zmatObj.getAtomData(1,"CA")
        me.CA = me.CAdata.generateAtom(2)
        
        me.Odata = me.zmatObj.getAtomData(1,"O")
        me.O = me.Odata.generateAtom(4)
        
        me.OXTdata = me.zmatObj.getAtomData(1,"OXT")
        me.OXT = me.OXTdata.generateAtom(5)

        me.Ndata = me.zmatObj.getAtomData(1,"N")
        me.N = me.Ndata.generateAtom(5)
        
        me.Hdata = me.zmatObj.getAtomData(1,"H")
        me.H = me.Hdata.generateAtom(5)

        me.tpHlist = []
        for rNum in me.referenceResNums:
            H = me.targetAll.getSpecificAtom(rNum,"H")
            N = me.targetAll.getSpecificAtom(rNum,"N")
            me.tpHlist.append((H,N))
    def run(self):
        HBond.HBond.feedConstants(self.cons)
        while(True):
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                #print ('%s: Exiting' % self.name)
                self.task_queue.task_done()
                break
            geo = next_task
            result = self.bestAngles(geo)
            self.result_queue.put(result)
        return
    def bestAngles(me,geo):
        import time
        sTime = time.time()
        bestScore = 0
        bestGeo = (None,None)
        me.C.setAtomLocationByZMAT(me.tpDist,me.tpAng,me.tpTor,geo[0] / 10.0,geo[1],geo[2])
        #print(geo)
        interest = []
        total = 0
        for ang in range(CaA[0],CaA[1]+1,CaA[2]):  
            for tor in range(CaT[0],CaT[1]+1,CaT[2]):
                relocateAtom(me.CA,me.CAdata,me.C,me.tpDist,me.tpAng,tor,ang)
                '''
                oRange = me.OSpread()
                print(oRange)
                for oTor in range(oRange[0],oRange[1]+1,OT[2]):
                '''
                for oTor in range(OT[0],OT[1]+1,OT[2]):
                    total+=1
                    relocateAtom(me.O,me.Odata,me.C,me.CA,me.tpDist,oTor)
                    relocateAtom(me.OXT,me.OXTdata,me.C,me.CA,me.O)
                    NOcontacts,HB_Score = me.numberNOContacts()
                    #print(NOcontacts,oTor)
                    if (NOcontacts >= 3):
                        interest.append((ang,tor,oTor,HB_Score,NOcontacts))
                    #me.dumpPose("Overlay.pdb",[me.C,me.O,me.OXT,me.CA])
                #print(noClash,attempt,total)
        if (len(interest) == 0):
            return None
        interest = sorted(interest,key = lambda x: x[3])
        bestScore = interest[0][3]
        interest = list(filter(lambda x : x[3] < bestScore * .75,interest))
        interest = uniqueCA(interest)
        for i in range(len(interest)):
            pose = interest[i]
            optPose = me.optimizePose(pose)
            interest[i] = optPose
            #input("\n\n")
        interest = sorted(interest,key = lambda x: x[3])      
        interest = uniqueCA(interest)
        #print("Time ",time.time()-sTime)

        if not allUnique(interest):
            print("Pose has uniqueness issues after optimization")
            for repose in interest:
                print(repose)
        #print("Best Geo ",bestGeo)


        for i in reversed(range(3)):
            interest = sorted(interest,key = lambda x: x[i])

        '''
        groups = []
        currentGroup = []
        for pose in interest:
            if (len(currentGroup) == 0):
                currentGroup.append(pose)
            else:
                for i in reversed(range(3)):
                    match = False
                    for member in currentGroup:
                        match = match or (abs(member[i] - pose[i]) <= 5)
                    if not match:
                        break
                if match:
                    currentGroup.append(pose)
                else:
                    groups.append(currentGroup)
                    currentGroup = [pose]
        if (len(currentGroup) > 0):
            groups.append(currentGroup)
        '''
        groups = []
        for pose in interest:
            match = False
            for curGroup in groups:
                for i in reversed(range(3)):
                    match = False
                    for member in curGroup:
                        match = match or (abs(member[i] - pose[i]) <= 5)
                    if not match:
                        break
                if match:
                    curGroup.append(pose)
                    break
            if not match:
                groups.append([pose])
        '''
        for i in range(len(groups)):
            print("Group ",i)
            A = groups[i]
            for B in A:
                print("\t",B)
        '''
        final = []
        for g in groups:
            final.append(min(g,key = lambda x:x[3]))
                        
        #print(len(final)," vs ",len(interest))       
        #print("Final")
        withN = []
        for i in range(len(final)):   
            bestGeo = final[i]
            #print("\t",bestGeo)
            withN+=me.seekN(geo+bestGeo)
        return withN
    def optimizePose(me,pose):
        caAng = pose[0]
        caTor = pose[1]
        oTor = pose[2]
        val = pose[3]
        change = True
        relocateAtom(me.CA,me.CAdata,me.C,me.tpDist,me.tpAng,pose[1],pose[0])
        oTor,val = me.optimize(pose,2)
        curPose = (caAng,caTor,oTor,val)
        #print("First O ",oTor,val)
        roundCount = 1
        
        while(change):
            change = False
            #print("\nRound %i"%roundCount)
            #print(curPose)
            #print("\t\t",me.CA.location)
            #print("\t\t",me.O.location)
            roundCount+=1
            priorVal = curPose[3]
            caAng,val = me.optimize(curPose,0)
            #print("\t\t",me.CA.location)
            #print("\t\t",me.O.location)
            curPose = (caAng,caTor,oTor,val)
            #print("CA ",caAng,val)
            oTor,val = me.optimize(curPose,2)
            curPose = (caAng,caTor,oTor,val)
            #print("CA O ",oTor,val)
            if (val - priorVal < -1):
                change = True

            priorVal = curPose[3] 
            caTor,val = me.optimize(curPose,1)
            curPose = (caAng,caTor,oTor,val)
            #print("caT ",caTor,val)
            oTor,val = me.optimize(curPose,2)
            curPose = (caAng,caTor,oTor,val)
            #print("caT O ",oTor,val)
            if (val - priorVal < -1):
                change = True
            #print("\t\t",me.O.location)
        return curPose
        
        
    def OSpread(me):
        low = None
        high = None
        for i in range(len(me.tpHlist)):
            H,N = me.tpHlist[i]
            relocateAtom(me.O,me.Odata,me.C,me.CA,H,0)
            tor = int(me.O.torsion(me.C,me.CA,me.tpDist)) % 180
            print("\t",tor)
            '''
            if (tor == 0):
                relocateAtom(me.OXT,me.OXTdata,me.C,me.CA,me.O)
                me.dumpPose("OSeekZero.pdb",[me.C,me.O,me.OXT,me.CA])
            '''
            if (low is None):
                low = tor
                high = tor
            if (tor < low):
                low = tor
            if (tor > high):
                high = tor
        print("LH ",low,high)
        input("\n")
        return (max(low - 20,0),min(high + 20,179))
    def seekN(me,geo):
        priorAngles = 6
        me.C.setAtomLocationByZMAT(me.tpDist,me.tpAng,me.tpTor,geo[0] / 10.0,geo[1],geo[2])
        relocateAtom(me.CA,me.CAdata,me.C,me.tpDist,me.tpAng,geo[4],geo[3])
        relocateAtom(me.O,me.Odata,me.C,me.CA,me.tpDist,geo[5])
        relocateAtom(me.OXT,me.OXTdata,me.C,me.CA,me.O)
        dist = me.Ndata.dist
        contacts = me.acceptorTree.radiusSearch(me.CA.location,me.cons["hBondDistCutoff"]+dist,True)
        #print(me.H.location)
        #input(me.H)
        tors = []
        for con,conDist in contacts:
            if (conDist < (0.75 + dist)):
                continue
            relocateAtom(me.N,me.Ndata,me.CA,me.C,con,0)
            seek = getAcceptorStem(con.residueType,con.atomType)
            planeSeek = getPlaneStem(con.residueType,con.atomType)
            if (planeSeek is None):
                planeStem = None
            else:
                planeStem = me.tRes[con.chain][con.residueNumber][planeSeek]
            conStem = me.tRes[con.chain][con.residueNumber][seek]
            me.H.setAtomLocationByZMAT(me.N,me.CA,con,1.01,120.0,0)
            if not me.clash([me.C,me.CA,me.O,me.OXT,me.N]):
                hb = HBond.evaluatePotentialHydrogenBond(me.N,con,planeStem,conStem,con,me.H,me.N)
                if (hb.good):
                    nTorsion = int(me.N.torsion(me.CA,me.C,me.tpDist))
                    if (nTorsion not in tors):
                        tors.append(nTorsion)
        results = []
        if (len(tors) > 0):
            seekResults = list(map(lambda x : ((geo[:priorAngles]+(x,)),geo[priorAngles]),tors))
            results+=seekResults
        
        for ntor in range(0,359,me.step["NGrain"]):
            newGeo = (geo[:priorAngles]+(ntor,))
            entry = (newGeo,geo[priorAngles])
            if (newGeo not in results):
                results.append(entry)
        return results
    def optimize(me,pose,column):
        def scoreTor(checkTor):
            if (column == 2):
                #print("\t\tCheck ",checkTor)
                
                relocateAtom(me.O,me.Odata,me.C,me.CA,me.tpDist,checkTor)
                relocateAtom(me.OXT,me.OXTdata,me.C,me.CA,me.O)
                #print("\t\t",me.O.location)
            else:
                if(column == 1):
                    relocateAtom(me.CA,me.CAdata,me.C,me.tpDist,me.tpAng,checkTor,pose[0])
                else:
                    #print("\t",checkTor,pose[0])
                    #print("\t",me.CA.location)
                    relocateAtom(me.CA,me.CAdata,me.C,me.tpDist,me.tpAng,pose[1],checkTor)
                    #print("\t",me.CA.location)
                relocateAtom(me.O,me.Odata,me.C,me.CA,me.tpDist,pose[2])
                relocateAtom(me.OXT,me.OXTdata,me.C,me.CA,me.O)
            return me.numberNOContacts()[1]
        def sanitize(tor):
            if (column == 1):
                return tor % 360
            if (tor < 1):
                return 1
            if (tor > 179):
                return 179
            return tor
        startTor = pose[column]
        centerVal = scoreTor(startTor)
        leftTor = sanitize(startTor-1)
        leftVal = scoreTor(leftTor)
        rightTor = sanitize(startTor+1)
        rightVal = scoreTor(rightTor)
        #print("\n\n")
        '''
        if (column == 2):
            print("\t\t",pose)
            print("\t\t",leftTor,startTor,rightTor)
            print("\t\t%5.1f, %5.1f, %5.1f"%(leftVal,centerVal,rightVal))
        '''
        increment = 0
        if (centerVal <= leftVal) and (centerVal <= rightVal):
            return startTor,scoreTor(startTor)
        if (rightVal < centerVal):
            if (rightVal < leftVal):
                increment = 1
                curTor = rightTor
                curVal = rightVal
            if (leftVal < rightVal):
                increment = -1
                curTor = leftTor
                curVal = leftVal
        elif (leftVal < centerVal):
            increment = -1
            curTor = leftTor
            curVal = leftVal
        if (increment == 0):
            if (leftTor == rightTor):
                if (rightTor <= startTor):
                    increment = -1
                    curTor = leftTor
                    curVal = leftVal
                elif (leftTor >= startTor):
                    increment = 1
                    curTor = rightTor
                    curVal = rightVal
        if (increment == 0):
            raise Exception("Logical case not considered %5.1f,%5.1f,%5.1f"%(leftVal,centerVal,rightVal))
        
        while(True):
            newTor = sanitize(curTor + increment)
            if (newTor == curTor):
                return curTor,scoreTor(curTor)
            newVal = scoreTor(newTor)
            #print("\t",newTor,newVal)
            if (newVal >= curVal):
                return curTor,scoreTor(curTor)
            curTor = newTor
            curVal = newVal
        
            
    def clash(me,atomList = None):
        if (atomList is None):
            atomList = [me.C,me.CA,me.O,me.OXT,me.N]
        for atm in atomList:
            neigh, closestDist = me.targetTree.nearestNeighbor(atm.location,True)
            if (closestDist <= me.cons["clashLimit"]):
                return True
        return False
    def numberNOContacts(me,debug=False):
        con = 0
        score = 0.0
        for atm in [me.O,me.OXT]:
            for i in range(len(me.tpHlist)):
                H,N = me.tpHlist[i]
                dist = H.distance(atm)
                if (dist >= MinDist) and (dist <= MaxDist):
                    V1 = H.location.difference(N.location)
                    VHO = H.location.difference(atm.location)
                    cosNHO = V1.cosTwoVectors(VHO)
                    if (cosNHO < 0):
                        cAng = me.C.angle(atm,H)
                        if (cAng > 100):             
                            V2 = atm.location.difference(me.C.location)
                            vectorCos = V1.dotProduct(V2) / (V1.magnitude() * V2.magnitude())
                            if (vectorCos < 0):
                                V3 = me.C.location.difference(me.CA.location)
                                planeNormal = V2.crossProduct(V3)
                                pnSin = math.sin(math.radians(planeNormal.angle(V1)))
                                geoScale = (pnSin * cosNHO)**2
                                hbEnergy = hbDistEnergy(dist,me.cons["directHBond"]) * geoScale
                                score+= hbEnergy
                                if (debug):
                                    print("\t%s->%s %2.2f, %5.1f, %5.1f"%(atm,N,geoScale,hbEnergy,score))
                                    print("\t%5.1f, %2.2f, %2.2f"%(hbDistEnergy(dist,me.cons["directHBond"]),cosNHO,pnSin))
                                con+=1
        if HBond.useCTerminalElectrostatic:
            for i in range(len(me.tpHlist)):
                H,N = me.tpHlist[i]
                dist1 = H.distance(me.O)
                dist2 = H.distance(me.OXT)
                dist = 0.5 * (dist1 + dist2)
                if (dist >= MinDist) and (dist <= MaxDist):
                    V1 = H.location.difference(N.location)
                    VHO = H.location.difference(me.C.location)
                    cosNHO = V1.cosTwoVectors(VHO)
                    if (cosNHO < 0):
                        cAng = me.CA.angle(atm,H)
                        if (cAng > 100):
                            V2 = atm.location.difference(me.C.location)
                            vectorCos = V1.dotProduct(V2) / (V1.magnitude() * V2.magnitude())
                            if (vectorCos < 0):
                                V3 = me.C.location.difference(me.CA.location)
                                planeNormal = V2.crossProduct(V3)
                                pnSin = math.sin(math.radians(planeNormal.angle(V1)))
                                geoScale = (pnSin * cosNHO)**2
                                hbEnergy = hbDistEnergy(dist,me.cons["directHBond"]) * geoScale
                                score+= hbEnergy
                                if (debug):
                                    print("\t%s->%s %2.2f, %5.1f, %5.1f"%(atm,N,geoScale,hbEnergy,score))
                                    print("\t%5.1f, %2.2f, %2.2f"%(hbDistEnergy(dist,me.cons["directHBond"]),cosNHO,pnSin))
                                con+=1
        return con,score          
    def dumpPose(me,filename,atomList):
        if (atomList is None):
            atomList = [me.C,me.CA,me.O,me.OXT,me.N,me.H]
        outf = open(filename,'a')
        for atm in atomList:
            outf.write("%s\n"%atm.toPDBLine())
        outf.write("TER\n")
        outf.close()
def relocateAtom(atm,data,aDist,aAng,aTor,tor = None, ang = None,dist = None):
    if (dist is None):
        dist = data.dist
    if (ang is None):
        ang = data.ang
    if (tor is None):
        tor = data.tor
    atm.setAtomLocationByZMAT(aDist,aAng,aTor,dist,ang,tor)            
def uniqueCA(ary):
    unique = set()
    ret = []
    for pose in ary:
        key = pose[:2]
        if key not in unique:
            ret.append(pose)
            unique.add(key)
    return ret
def allUnique(ary):
    unique = set(ary)
    if (len(unique) != len(ary)):
        return False
    return True
if __name__ == "__main__":
    for x in range(0,50):
        d = x / 10.0
        print(d,"%4.1f"%energy(d,-1200))
    '''
    for n in range(CD[0],CD[1],CD[2]):
        print(n)
    '''
