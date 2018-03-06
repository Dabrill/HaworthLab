recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
class logbook():
    def __init__(self):
        self.counter = 0
        self.redundant = 0
        self.perfect = True
        self.outf =open("Redundancies.txt",'w')
        self.unique = []
        self.redundant = []
    def end(self):
        if (self.perfect):
            PRINT("All solutions unique")
            self.outf.write("All solutions unique\n")
        else:
            PRINT("Percent redundant "+str(len(self.redundant)*100 / self.counter)+"%")
            self.redundant = sorted(self.redundant,key=lambda group: group[0].indexNum)
            for array in self.redundant:
                self.outf.write("."+str(array[0].indexNum)+".\n")
                for j in range(1,len(array)):
                    self.outf.write("\t."+str(array[j].indexNum)+".\n")
            self.unique.sort()
            self.outf.write("\n\n\nUnique solutions:\n")
            for z in self.unique:
                self.outf.write("\t."+str(z)+".\n")
        self.outf.close()
def redundancyCheck(data):
    redTable = []
    for n in range(1,data.getLength()+1):
        redTable.append(data.getCoordMatrixPeptideObject(n))
    PRINT("Starting redundancy check")
    counter = logbook()
    PRINT("Peptides for grouping: "+str(len(redTable)))
    procTable = groupPepRed(redTable,0)
    PRINT("Grouping Done")
    evaluatePepRed(procTable,counter)
    PRINT("Calculation done. Making output file.")
    if (counter.counter == len(redTable)):
        PRINT("Results have same number of poses as input")
    else:
        PRINT("Mismatch between input and output pose number")
        PRINT("Counter: "+str(counter.counter)+" Table: "+str(len(redTable)))
    counter.end()
    PRINT("Clearing memory")
    notRedundant = counter.perfect
    counter = None
    redTable = None
    procTable = None
    PRINT("Cull and redundancy routine done")
    return notRedundant
def groupPepRed(array,pos):
    if (len(array) < 1): return
    if (pos >= len(array[0].atoms)): return array
    array = sorted(array,key=lambda pep: pep.atoms[pos][0])
    last = 0
    #sa = SameAtom
    saGroup = 0
    saValue = array[0].atoms[pos][0]
    groups = [[array[0]]]

    for i in range(1,len(array)):
        if (array[i].atoms[pos][0] == saValue):
            found = False
            for j in range(saGroup,len(groups)):
                if (array[i].sameAtom(groups[j][0],pos)):
                    groups[j].append(array[i])
                    found = True
            if (not found):
                newGroup = [array[i]]
                groups.append(newGroup)
        else:
            newGroup = [array[i]]
            groups.append(newGroup)
            saValue = array[i].atoms[pos][0]
            saGroup=len(groups)-1
    if (len(groups) < len(array)):
        uberGroup = []
        for i in range(0,len(groups)):
            uberGroup.append(groupPepRed(groups[i],pos+1))
        return uberGroup
    if (len(groups) == len(array)):
        uber = []
        for g in groups:
            uber.append([g])
        return uber  
    return "error"
def evaluatePepRed(array,counter):
    try:
        a = len(array[0])
        for i in range(0,len(array)):
            evaluatePepRed(array[i],counter)
    except:
        counter.counter+=len(array)
        if (len(array) > 1):
            if (counter.perfect):
                PRINT("\n##################\nRedundancy found\n##################\n")
                counter.perfect = False
            counter.redundant.append(array)
        else:
            counter.unique.append(array[0].indexNum)
