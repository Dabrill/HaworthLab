import PDBTools
import Vector
import math
import KDTree
#Dab Exsan Mode
maxWaterDistance = 4.0
waterChains = ["W"]
#Jason Tetrahedral Mode
#maxWaterDistance = 10.0
#waterChains = ["W","I"]
class WaterLayer:
    def __init__(me,filename,chain):
        inf = open(filename)
        me.tp = []
        watgenOxy = []
        watgenWaterAtoms = []
        me.header = inf.read()
        inf.seek(0)
        for line in inf:
            atm = PDBTools.Atom(line)
            if (atm.chain in chain) and (not atm.isHydrogen()):
                me.tp.append(atm)
            if (atm.chain in waterChains):
                watgenWaterAtoms.append(atm)
                if (atm.atomType == "O"):
                    watgenOxy.append(atm)        
        me.proteinTree = KDTree.KDTree.loadAtomArray(me.tp)        
        me.tpTree = KDTree.KDTree.loadAtomArray(me.tp+watgenOxy)
        me.watgenWaters = []
        for i in range(0,len(watgenWaterAtoms),3):
            me.watgenWaters.append(Water.grabFromArray(watgenWaterAtoms,i))
        me.generatedWaters = []
        me.writeWaters = []
        me.waterCount = watgenOxy[len(watgenOxy)-1].residueNumber
        me.runAdded = 0
    def createCandidatesForWater(me,inWater):
        def makeWater(torsion):
            O = makeNewOxygen(torsion)
            H1 = newBlankHydrogen(1)
            H1.setAtomLocationByZMAT(O,inWater.O,inWater.H1,1.00,0.01,180)
            H2 = newBlankHydrogen(2)
            H2.setAtomLocationByZMAT(O,H1,inWater.O,1.00,109.4,70)
            ret = Water(O,H1,H2)
            ret.type = "CANDIDATE"
            ret.donors = [inWater]
            return ret
        def newBlankHydrogen(number):
            ret = PDBTools.Atom.blankAtom()
            ret.residueType = "GEN"
            ret.residueNumber = -1
            ret.atomType = "H%i"%number
            ret.atomNumber = 0
            ret.chain = "W"
            return ret
        def makeNewOxygen(torsion):
            ret = PDBTools.Atom.blankAtom()
            ret.residueType = "GEN"
            ret.residueNumber = -1
            ret.atomType = "O"
            ret.atomNumber = 0
            ret.chain = "W"
            ret.setAtomLocationByZMAT(inWater.O,inWater.H1,inWater.H2,2.76,109.4,torsion)
            return ret
        return makeWater(90),makeWater(-90)
    def validateCandidate(me,candidate):
        neigh,dist = me.tpTree.nearestNeighborWithDist(candidate.O.location)
        neighs = me.tpTree.radiusSearch(candidate.O.location,1.5)
        #if (len(neighs) == 0):
        if (dist >= 1.5):
            protNeigh, protDist = me.proteinTree.nearestNeighborWithDist(candidate.O.location)
            if (protDist <= maxWaterDistance):
                watNeigh = me.closestGeneratedWater(candidate,2.0)
                if (watNeigh is None):
                    me.runAdded+=1
                    return True
                else:
                    #print("\n\n%s\nis close to a candidate water\n%s"%(neigh.O.toPDBLine(),candidate.O.toPDBLine()))
                    watNeigh.averageWater(candidate)
            
        return False
    def neighbors(me,candidate,dad):
        neigh,dist = me.tpTree.nearestNeighborWithDist(candidate.O.location)
        protNeigh, protDist = me.proteinTree.nearestNeighborWithDist(candidate.O.location)
        watNeigh = me.closestGeneratedWater(candidate,1.5)
        print("\n%s:"%candidate.O.residueNumber)
        print("Spawned by %s"%dad.O)
        print("%s at %f"%(neigh,dist))
        print("%s at %f"%(protNeigh,protDist))
        if (watNeigh is not None):
            print("%s at %f"%(watNeigh,watNeigh.distance(candidate.O)))
    def closestGeneratedWater(me,location,cutoff=None,includeDistance=False):
        closest = None
        closestDistance = None
        for w in me.generatedWaters:
            dist = w.distance(location)
            if ((cutoff is None) or (dist < cutoff)) and ((closest is None) or (dist < closestDistance)):
                closest = w
                closestDistance = dist
        if (includeDistance):
            return closest,closestDistance
        else:
            return closest
    def prepareNewRound(me):
        me.watgenWaters = me.watgenWaters + me.generatedWaters
        me.generatedWaters = []
        me.runAdded = 0
        watgenOxy = []
        for w in me.watgenWaters:
            watgenOxy.append(w.O)
        me.tpTree = KDTree.KDTree.loadAtomArray(me.tp+watgenOxy)
        #me.dumpList(me.tp+watgenOxy,"Round%iDump.txt"%me.roundCount)
    def dumpList(me,ary,toFile):
        f = open(toFile,'w')
        for l in ary:
            f.write("%s\n"%l.toPDBLine())
        f.close()
    def generateWaterLayer(me):
        me.roundCount = 1
        lock = True
        evalWaters = me.watgenWaters
        while (lock):
            for water in evalWaters:
                candidates = me.createCandidatesForWater(water)
                for c in candidates:
                    if (me.validateCandidate(c)):
                        me.waterCount+=1
                        c.accept(me.waterCount)
                        me.generatedWaters.append(c)
            
            for water in me.generatedWaters:
                water.optimize()

            me.roundCount+=1
            lock = (me.runAdded > 0)
            evalWaters = me.generatedWaters
            me.writeWaters = me.writeWaters + me.generatedWaters
            if (lock):
                me.prepareNewRound()
        print("Added %i waters in %i rounds"%(len(me.writeWaters),me.roundCount))
        
    def outputFile(me,filename):
        outf = open(filename,"w")
        outf.write(me.header)
        for i in range(len(me.writeWaters)):
            w = me.writeWaters[i]
            #w.renumber(i+1)
            outf.write("%s\n"%w)
        outf.close()
class Water:
    def __init__(me,O,H1,H2):
        me.O = O
        me.H1 = H1
        me.H2 = H2
        me.type = "DEFAULT"
        me.averages = []
    @staticmethod
    def grabFromArray(array,start):
        newbie = Water(array[start],array[start+1],array[start+2])
        newbie.type = "WATGEN"
        return newbie
    def distance(me,other):
        return me.O.distance(other.O)
    def accept(me,resNumber):
        me.type = "GENERATED"
        me.averages=[me.O.location.dims]
        for atm in [me.O,me.H1,me.H2]:
            atm.residueNumber = resNumber
    def __str__(me):
        return "%s\n%s\n%s"%(me.O.toPDBLine(),me.H1.toPDBLine(),me.H2.toPDBLine())
    def averageWater(me,other):
        me.averages.append(other.O.location.dims)
        coord = [0.0,0.0,0.0]
        for d in range(3):
            total = 0.0
            for avg in me.averages:
                total+=avg[d]
            coord[d] = total / len(me.averages)
        #print(me.O.residueNumber,len(me.averages))
        me.O.location.dims = coord
        me.donors.append(other.donors[0])
        '''
        for don in me.donors:
            print("\t%s"%don.O)
        '''
    def averageWaterVerbose(me,other):
        outf = open("AverageWater.pdb",'w')
        outf.write("%s\n"%me.donors[0])
        outf.write("%s\n"%other.donors[0])
        outf.write("%s\n"%me)
        outf.write("%s\n"%other)
        me.averages.append(other.O.location.dims)
        coord = [0.0,0.0,0.0]
        for d in range(3):
            total = 0.0
            for avg in me.averages:
                total+=avg[d]
            coord[d] = total / len(me.averages)
        #print(me.O.residueNumber,len(me.averages))
        me.O.location.dims = coord
        me.donors.append(other.donors[0])
        '''
        for don in me.donors:
            print("\t%s"%don.O)
        '''

        #print("\t%d"%(me.donors[0].O.angle(me.O,me.donors[1].O)))
        
        me.optimize()
        outf.write("%s\n"%me)
        outf.close()
    def optimize(me):
        def score(AH1,TH1,TH2):
            H1.setAtomLocationByZMAT(O,D1,HD1,1.00,AH1,TH1)
            H2.setAtomLocationByZMAT(O,H1,D1,1.00,109.4,TH2)
            AH2 = O.angle(H2,D2)
            ret =  math.cos(math.radians(AH2))-math.cos(math.radians(AH1))
            if (best is None) or (best > ret):
                pass
                #print(AH1,AH2,ret,best)
            return ret
        if (len(me.averages) == 1):
            return
        #print(me.O,len(me.donors))
        twoDegreeCos = -math.cos(math.radians(2.0))
        O = me.O
        H1 = me.H1
        H2 = me.H2
        bestParam = None
        best = None
        for i in range(len(me.donors)):
            for j in range(i+1,len(me.donors)):
                D1 = me.donors[i].O
                HD1 = me.donors[i].H1
                D2 = me.donors[j].O
                for iAH1 in range(2,90,8):
                    maxPossible = twoDegreeCos - math.cos(math.radians(iAH1))          
                    if (best is not None) and (best < maxPossible):
                        #print("Best found")
                        break
                        #pass
                    for iTH1 in range(0,360,30):
                        for iTH2 in range(0,360,30):
                            myScore = score(iAH1,iTH1,iTH2)
                            #print(myScore,best)
                            if (best is None) or (best > myScore):
                                best = myScore
                                bestParam = (iAH1,iTH1,iTH2)
                                #print("OK ",best,myScore)

        H1.setAtomLocationByZMAT(O,D1,HD1,1.00,bestParam[0],bestParam[1])
        H2.setAtomLocationByZMAT(O,H1,D1,1.00,109.4,bestParam[2])
        cAH2 = O.angle(H2,D2)
        #print("\t",bestParam,cAH2)
        #print(math.cos(math.radians(cAH2)),math.cos(math.radians(bestParam[0])))
if __name__ == "__main__":
    extendWaterNetwork()
