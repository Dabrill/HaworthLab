recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
class EnergyTable:
    TABLE_FILENAME = "ETable.tsv"
    def __init__(me):
        me.table = []
        me.types = [int,float,int,float,int,float,int,int,int,int,int,int,float]
        me.colToParam = ["Pose","Energy","#MissWat","EMissWat","#HBond","EHBond","#HyPhobic","#BadPocket","#BulkFace","#Salt","#LysDP","#PAmide","ESide"]
        me.width = len(me.colToParam)
        me.paramToCol = {}
        me.lastWritten = -1
        for i in range(me.width):
            me.paramToCol[me.colToParam[i]] = i
        me.outf = None
        me.filename = None

    @staticmethod
    def makeBlankTable(size,filename=None,prepareOutput=True):
        if filename is None:
            filename = EnergyTable.TABLE_FILENAME
        new = EnergyTable()
        for i in range(size):
            line = [None] * new.width
            line[new.paramToCol["Pose"]] = i + 1
            new.table.append(line)
        if prepareOutput:
            new.setOutf(filename)
        return new
    @staticmethod
    def readTableFile(filename=None):
        if filename is None:
            filename = EnergyTable.TABLE_FILENAME
        new = EnergyTable()
        new.filename = filename
        inf = open(filename,'r')
        topline = inf.readline()
        header = topline.split('\t')
        normalMode = True
        if (len(header) > new.width):
            normalMode = False
            new.addControlCols()
            #raise NotImplementedError("Nonstandard table. Reading extended tables not implemented")
        for line in inf:
            a = line.split('\t')
            for j in range(new.width):
                if (a[j][:3] == "N/A"):
                    val = None
                else:
                    val = EnergyTable.cast(a[j],new.types[j])
                a[j] = val
            new.table.append(a)
        return new
    @staticmethod
    def cast(value,toType):
        if (toType == str):
            return str(value)
        if (toType == float):
            return float(value)
        if (toType == int):
            return int(value)
        return None
    def addControlCols(me):
        me.addCols(["E_Rank","RMS_Rank","RMS","RMS_BB_Rank","RMS_BB"],[float,float,float,float,float])
    def addCols(me,labels,typing):
        if (len(labels) != len(typing)):
            raise Exception("Must be same number of column labels and type indictators")
        else:
            size = len(labels)
        me.types = me.types + typing
        me.colToParam = me.colToParam + labels
        me.width = len(me.colToParam)
        me.paramToCol = {}
        for i in range(me.width):
            me.paramToCol[me.colToParam[i]] = i
        for i in range(me.length()):
            for j in range(size):
                me.table[i].append(None)
            
    def dump(me):
        s = me.makeHeader()
        for row in range(me.length()):
            s+=me.stringifyLine(row)
        return s              
    def makeHeader(me):
        s = ""
        for i in range(me.width):
            s+=me.colToParam[i]
            if (i < me.width - 1):
                s+="\t"
            else:
                s+="\n"
        return s
    def length(me):
        return len(me.table)
    def stringifyLine(me,row):
        s = ""
        for i in range(me.width):
            val = me.table[row][i]
            if val is not None:
                s+=str(val)
            else:
                s+="N/A"
            if (i < me.width - 1):
                s+="\t"
            else:
                s+="\n"
        return s
    def printLine(me,row):
        s = me.stringifyLine(row)
        me.outf.write(s)
        me.lastWritten = row
    def feedParameter(me,value,paramName,row):
        row = row - 1
        cell = me.table[row][me.paramToCol[paramName]]
        #if (cell.getParameter["Pose"] == row):
        #if cell is not None:
        #print("Value overwrite Row"+str(row)+", Col"+str(me.paramToCol[paramName])+"/"+paramName)
        me.table[row][me.paramToCol[paramName]] = value
    def getParameter(me,paramName,row):
        row = row - 1
        return me.table[row][me.paramToCol[paramName]]
    def doneFeeding(me):
        while (me.lastWritten < me.length() - 1):   
            me.printLine(me.lastWritten+1)
        me.outf.close()
    def setOutf(me,filename):
        me.lastWritten = -1
        if me.outf is not None:
            me.outf.close()
        me.outf = open(filename,'w')
        me.outf.write(me.makeHeader())
        me.filename = filename
    def writeToFile(me,filename = None):
        if filename is not None:
            me.setOutf(filename)
        else:
            me.outf.close()
            me.outf = open(me.filename,'w')
        for i in range(me.length()):
            me.printLine(i)
    def rankTable(me,colInLabel,colOutLabel):
        def resolveTie(table,tieStart,tieEnd):
            tieRank = (tieStart+tieEnd+2.0)/2.0
            for y in range(tieStart,tieEnd+1):
                table[y][colOut] = tieRank
        colIn = me.paramToCol[colInLabel]
        print("ColIn ",colIn)
        colOut = me.paramToCol[colOutLabel]
        me.table = sorted(me.table,key=lambda j: colVal(j,colIn))
        tieMode = False
        for i in range(len(me.table)):
            if (i > 0):
                if (me.table[i][colIn] == me.table[i-1][colIn]):
                    if not tieMode:
                        tieStart = i-1
                    tieMode = True
                else:
                    if (tieMode):
                        tieEnd = i - 1
                        resolveTie(me.table,tieStart,tieEnd)
                    tieMode = False
            me.table[i][colOut] = float(i+1)
        if (tieMode):
            resolveTie(me.table,tieStart,len(me.table)-1)
        me.table = sorted(me.table,key=lambda j: j[0])
    def makeSubtable(me,filename,poses):
        import copy
        sub = EnergyTable()
        sub.table = []
        sub.types = copy.copy(me.types)
        sub.colToParam = copy.copy(me.colToParam)
        sub.width = len(me.colToParam)
        sub.paramToCol = copy.copy(me.paramToCol)
        sub.lastWritten = -1
        sub.outf = open(filename,'w')
        sub.filename = filename
        sub.outf.write(me.makeHeader())
        for p in poses:
            sub.table.append(me.table[p-1])
        return sub
def colVal(row,col):
    v = row[col]
    if (v is None):
        return 9999999
    return v
def extendETable(CONSTANTS,poseData,anchor):
    import PDBTools
    table = EnergyTable.readTableFile("Etable.tsv")
    table.addControlCols()
    PRINT("ETable loaded")
    table.rankTable("Energy","E_Rank")
    
    model = PDBTools.readPDBFile(CONSTANTS["commonFolder"]+"/"+CONSTANTS["knownStruc"])[0]
    firstPep = poseData.getPeptide(0)
    overlaps = model.findOverlaps(firstPep)
    start = PDBTools.pickOverlap(overlaps,anchor)
    if (start < 0):
        PRINT("Problem aligning peptide to X-ray structure. Number of potential overlaps: "+str(len(overlaps)))
        if (len(overlaps) == 0):
            PRINT("Peptide: "+firstPep.getSequence())
            PRINT("X-ray: "+model.getSequence())
    
    model.makeDictionary()
    PRINT("PDB loaded")
    missingAtoms = firstPep.compareMissing(model,start)
    if(len(missingAtoms) > 0):
        PRINT("Following atoms could not be located in Xray control:")
        for ma in missingAtoms:
            PRINT("\t"+ma.identifier())
    else:
        PRINT("All atoms found")
    
    if (poseData.getLength() != table.length()):
        PRINT("Mismatch in length in extending ETable")
        PRINT("ETable: "+str(table.length()))
        PRINT("Pose Data "+str(poseData.getLength()))

    for i in range(table.length()):
        poseNum = table.getParameter("Pose",i)
        pep = poseData.getPeptide(poseNum)
        if (start >= 0):
            LMS, LMS_BB = pep.compare(model,start)
            table.feedParameter(LMS,"RMS",i)
            table.feedParameter(LMS_BB,"RMS_BB",i)
        else:
            table.feedParameter(start,"RMS",i)
            table.feedParameter(start,"RMS_BB",i)
        if (i% 25000 == 0):
            print(str(i))
    PRINT("RMS Calculated")
    table.rankTable("RMS","RMS_Rank")
    PRINT("RMS Ranked")
    table.rankTable("RMS_BB","RMS_BB_Rank")
    PRINT("RMS Backbone Ranked")
    table.writeToFile("ConTable.tsv")
    table.doneFeeding() 
