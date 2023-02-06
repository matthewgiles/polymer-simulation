#!/usr/bin/env python3

# Program to read in a dump file frame by frame
# and do some calculations.

# Example here is radius of gyration
# Here Rg of all atoms in the file are calculated

# Assumes the dump file has the following format
# id type x y z ix iy iz
# were x  etc. are atom coordinates
#      ix etc. are image flags for taking periodic boundaries into account
# need to adjust these to range from -L/2 to L/2


import numpy as np
import operator

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def unit(v):
    return v / np.linalg.norm(v)

def cross(x, y):
    return np.cross(x, y)

def dot(x, y):
    return np.dot(x, y)

class Quaternion(object):
    def __init__(self, q0, q1, q2, q3):
        self.__q0 = q0
        self.__q1 = q1
        self.__q2 = q2
        self.__q3 = q3
        self.__x = self.__y = self.__z = None
    
    def xAxis(self):
        if self.__x is None:
            x1 = self.__q0 * self.__q0 + self.__q1 * self.__q1 - self.__q2 * self.__q2 - self.__q3 * self.__q3
            x2 = 2.0 * self.__q1 * self.__q2 + 2.0 * self.__q0 * self.__q3
            x3 = 2.0 * self.__q1 * self.__q3 - 2.0 * self.__q0 * self.__q2
            self.__x = np.array([x1, x2, x3])
        return self.__x
    
    def yAxis(self):
        if self.__y is None:
            y1 = 2.0 * self.__q1 * self.__q2 - 2.0 * self.__q0 * self.__q3
            y2 = self.__q0 * self.__q0 - self.__q1 * self.__q1 + self.__q2 * self.__q2 - self.__q3 * self.__q3
            y3 = 2.0 * self.__q2 * self.__q3 + 2.0 * self.__q0 * self.__q1
            self.__y = np.array([y1, y2, y3])
        return self.__y
    
    def zAxis(self):
        if self.__z is None:
            z1 = 2.0 * self.__q1 * self.__q3 + 2.0* self.__q0 * self.__q2
            z2 = 2.0 * self.__q2 * self.__q3 - 2.0 * self.__q0 * self.__q1
            z3 = self.__q0 * self.__q0 - self.__q1 * self.__q1 - self.__q2 * self.__q2 + self.__q3 * self.__q3
            self.__z = np.array([z1, z2, z3])
        return self.__z

class Atom(object):
    """ A Class for storing atom information """

    def __init__(self, id, type, rx, ry, rz, ix = 0, iy = 0, iz = 0, q0 = 0, q1 = 0, q2 = 0, q3 = 0):
        """ Initialise the class """
        self.id = id                             # id of the atom
        self.type = type                         # type of the atom
        self.r = np.array([rx, ry, rz], dtype = np.float64)     # position of the atom
        self.image = np.array([ix, iy, iz], dtype = np.int32)         # image flags for atoms
        self.quaternion = Quaternion(q0, q1, q2, q3)
        self.unwrap_flag = False

    def sep(self, B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.r[0]-B.r[0])**2 + 
                        (self.r[1]-B.r[1])**2 + (self.r[2]-B.r[2])**2 )

    def minus(self,B):
        """ Subtract B.r vector from self.r vector """
        return self.r - B.r

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def q(self):
        return self.quaternion

    def unwrap(self, L):
        """ Unwraps the coordinates for periodic box, and overwrites x """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.r[j] = self.r[j] + self.image[j]*L[j] # unwrap
            self.unwrap_flag = True


class Processor(object):

    def __init__(self):
        self.__file = self.__N = self.__dump = None
    
    def set(self, file, N):
        self.__file = file
        self.__N = N

    def noLines(self):
        """ Get the number of lines in the file """

        with open(self.__file, 'r') as dump:
            for i, _ in enumerate(dump):
                pass

        return i + 1
    
    def readFrame(self):
        atoms = []
        L = []

        for i in range(9):
            line = self.__dump.readline().split()
            if i >= 5 and i <= 7:
                # get the box size
                L.append(np.float64(line[1]) - np.float64(line[0]))
        
        for i in range(self.__N):
            line = self.__dump.readline().split()
            a = Atom(*line[0:12])
            a.unwrap(L)
            atoms.append(a)
            
        atoms.sort(key = operator.attrgetter('id'))

        return atoms, L

    def process(self, outputFile = 'r_g_1.dat'):
        if self.__file is None:
            return
        
        noFrames = int(self.noLines() / (self.__N + 9))
        
        out = open(outputFile, 'w')  
        out.write("# frame number, radius of gyration\n")
        
        self.__dump = open(self.__file, 'r')
        
        for frame in range(noFrames):
            atoms, L = self.readFrame()
            rG = self.radiusOfGyration(atoms)
            #Tw = self.twist(atoms)
            out.write("%i %.5f\n" % (frame + 1, rG))

        self.__dump.close()
        out.close()

    def radiusOfGyration(self, atoms):
        r_mean = np.array([0.0, 0.0, 0.0])
        # Unwrap each atom and calculate the CoM by summing the position of all atoms, and
        # dividing by the number of atoms
        for a in atoms:
            r_mean += a.r
        r_mean /= self.__N
        r_mean = Atom(-1, -1, r_mean[0], r_mean[1], r_mean[2])
        
        r_2 = 0.0
        # Radius of gyration squared is the sum of the squares of the separation of each atom and the CoM
        for a in atoms:
            r_2 += a.sep(r_mean) ** 2
        return np.sqrt(r_2/self.__N)

    def twist(self, atoms):
        tangents = []
        for i in range(0, self.__N - 1):
            tangents.append(atoms[i + 1].r - atoms[i].r)
            tangents[i] = unit(tangents[i])
        
        tangents.append(atoms[0].r - atoms[self.__N - 1].r)
        unit(tangents[self.__N - 1])

        perpendiculars = []
        for i in range(0, self.__N):
            perpendiculars.append(unit(cross(tangents[i], atoms[i].q().xAxis())))
        
        binomials = [unit(cross(tangents[self.__N - 1], tangents[0]))] 
        angles = [np.arccos(dot(tangents[self.__N - 1], tangents[0]))]
        mshift = [dot(perpendiculars[self.__N - 1], rotation_matrix(binomials[0], angles[0]))]
        for i in range(1, self.__N):
            binomials.append(unit(cross(tangents[i - 1], tangents[i])))
            angles.append(np.arccos(dot(tangents[i - 1], tangents[0])))
            mshift.append(dot(perpendiculars[i - 1], rotation_matrix(binomials[i], angles[i])))

        Tw = 0
        for i in range(0, self.__N - 1):
            Tw += dot(tangents[i], cross(mshift[i], perpendiculars[i]))
        
        Tw /= (2.0 * np.pi)





p = Processor()
p.set('dump.DNA', 100)
p.process()