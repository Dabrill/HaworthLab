import math
class Vector3D:
    def __init__(self,first,Y=None,Z=None):
        if isinstance(first,list):
            self.dims = first
        else:
            self.dims = [first,Y,Z]
    def X(self):
        return self.dims[0]
    def Y(self):
        return self.dims[1]
    def Z(self):
        return self.dims[2]
    def copy(self):
        return Vector3D(self.dims[0],self.dims[1],self.dims[2])
    def toArray(self):
        return self.dims
    def translate(self,other):
        ret = [None] * 3
        if isinstance(other,Vector3D):
            for i in range(3):
                ret[i] = self.dims[i]+other.dims[i]
        else:
            for i in range(len(other)):
                ret[i] = self.dims[i]+ other[i]
        return Vector3D(ret)
    def unitVector(self):
        mag = self.magnitude()
        return Vector3D(self.dims[0]/mag,self.dims[1]/mag,self.dims[2]/mag)
    def dotProduct(me,other):
        return me.dims[0]*other.dims[0]+me.dims[1]*other.dims[1]+me.dims[2]*other.dims[2]
    def crossProduct(a,b):
        x = a.dims[1]*b.dims[2] - a.dims[2] *b.dims[1]
        y = a.dims[2]*b.dims[0] - a.dims[0] *b.dims[2]
        z = a.dims[0]*b.dims[1] - a.dims[1] *b.dims[0]
        return Vector3D([x,y,z])
    def pointInPlane(point,planeA,planeB,planeC):
        v1 = point.difference(planeA)
        l1 = planeA.difference(planeB).unitVector()
        l2 = planeB.difference(planeC).unitVector()
        Normal = l1.crossProduct(l2).unitVector()
        dist = v1.dotProduct(Normal)
        return point.translate(Normal.scalarMultiplication(-dist))
    '''
    def planeIntersectionThreePoints(V1,planeA,planeB,planeC):
        V1 = planeA.difference(planeB).unitVector()
        V2 = planeB.difference(planeC).unitVector()
        Normal = V1.crossProduct(V2).unitVector()
        return V1.planeIntersectionVector(Normal)
    def planeIntersectionVectorAndPoint(V1,Normal,point):
        V2 = 

        
        dist = V1.dotProduct(Normal)
        return point.translate(Normal.scalarMultiplication(-dist))
    '''
    def planeNormal(a,b,c):
        vert1 = a.difference(b).unitVector()
        vert2 = c.difference(b).unitVector()
        return vert1.crossProduct(vert2)
    def scalarMultiplication(self,scalar):
        return Vector3D(self.dims[0]*scalar,self.dims[1]*scalar,self.dims[2]*scalar)
    def angleThreePoints(a,b,c,rad=False):
        cos = a.cosThreePoints(b,c)
        if (rad):
            return math.acos(cos)
        return math.degrees(math.acos(cos))
    def cosTwoVectors(a,b):
        return a.unitVector().dotProduct(b.unitVector())
    def cosThreePoints(a,b,c):
        vert1 = a.difference(b).unitVector()
        vert2 = c.difference(b).unitVector()
        return vert1.dotProduct(vert2)
    def dihedral(self,a,b,c,roundAnswer = True):
        b1 = self.difference(a).unitVector()
        b2 = a.difference(b).unitVector()
        b3 = b.difference(c).unitVector()
        n1 = b1.crossProduct(b2)
        n2 = b2.crossProduct(b3)
        m = n1.crossProduct(b2)

        x = n1.dotProduct(n2)
        y = m.dotProduct(n2)
        answer = math.degrees(math.atan2(y,x))
        if (roundAnswer):
            return cleanAngle(answer)
        return answer
    @staticmethod
    def zmatToCartesian(aDist,aAng,aTor,distance,angle,torsion):
        if ((angle <= 0) or (angle >= 180)):
            raise ValueError("Angle must be 0<angle<180 degrees")
        angleRad = math.radians(angle)
        torsionRad = math.radians(torsion)
        l1 = aAng.difference(aTor).unitVector()
        l2 = aDist.difference(aAng).unitVector()
        vp = l1.crossProduct(l2)
        norm = math.sqrt(1-l1.dotProduct(l2)**2)
        l3 = vp.scalarMultiplication(1/norm)
        l4 = l3.crossProduct(l2)

        cosA = math.cos(angleRad)
        sinA = math.sin(angleRad)
        cosD = math.cos(torsionRad)
        sinD = math.sin(torsionRad)
        sinAcosD = sinA * cosD
        sinAsinD = sinA * sinD

        vj = Vector3D([None,None,None])
        for i in range(3):
            vj.dims[i] = distance*(-l2.dims[i]*cosA + l4.dims[i]*sinAcosD + l3.dims[i]*sinAsinD)
        return vj.translate(aDist)
    
    def angle(self,other,rad=False):
        dotP = self.dotProduct(other) / (self.magnitude() * other.magnitude())
        if (rad):
            return math.acos(dotP)
        else:
            return math.degrees(math.acos(dotP))
    def negative(self):
        return Vector3D(-self.dims[0],-self.dims[1],-self.dims[2])
    def magnitudeSQ(self):
        return self.dims[0]**2+self.dims[1]**2+self.dims[2]**2
    def magnitude(self):
        return math.sqrt(self.dims[0]**2+self.dims[1]**2+self.dims[2]**2)
    def difference(self,other):
        return Vector3D(self.dims[0]-other.dims[0],self.dims[1]-other.dims[1],self.dims[2]-other.dims[2])
    def distance(self,other):
        return self.difference(other).magnitude()
    def distanceSQ(self,other):
        return (self.dims[0]-other.dims[0])**2+(self.dims[1]-other.dims[1])**2+(self.dims[2]-other.dims[2])**2

    def toPDBFormat(self):
        return formatNumber(self.dims[0])+" "+formatNumber(self.dims[1])+" "+formatNumber(self.dims[2])
    def cleanCoords(self):
        return Vector3D(float(formatNumber(self.dims[0])),float(formatNumber(self.dims[1])),float(formatNumber(self.dims[2])))
    def __str__(me):
        return "["+str(me.dims[0])+","+str(me.dims[1])+","+str(me.dims[2])+"]"
def cleanAngle(d):
    d = d%360.0
    if (d > 359.9):
        d = 0
    return round(d,0)
def formatNumber(d):
    s = str(d)
    dot = s.find('.')
    if (dot == -1):
        s+="."
        dot = s.find('.')
    for i in range(0,3-dot):
        s=" "+s
    l = len(s)
    for i in range(0,7-l):
        s=s+"0"
    return s[:7]
