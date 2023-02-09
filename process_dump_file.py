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
from supercoiling_cython import calculate_Writhe

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
        self.__f = self.__v = self.__u = None

    def fAxis(self):
        if self.__f is None:
            f1 = self.__q0 * self.__q0 + self.__q1 * self.__q1 - self.__q2 * self.__q2 - self.__q3 * self.__q3
            f2 = 2.0 * self.__q1 * self.__q2 + 2.0 * self.__q0 * self.__q3
            f3 = 2.0 * self.__q1 * self.__q3 - 2.0 * self.__q0 * self.__q2
            self.__f = np.array([f1, f2, f3])
        return self.__f
    
    def yAxis(self):
        if self.__v is None:
            v1 = 2.0 * self.__q1 * self.__q2 - 2.0 * self.__q0 * self.__q3
            v2 = self.__q0 * self.__q0 - self.__q1 * self.__q1 + self.__q2 * self.__q2 - self.__q3 * self.__q3
            v3 = 2.0 * self.__q2 * self.__q3 + 2.0 * self.__q0 * self.__q1
            self.__v = np.array([v1, v2, v3])
        return self.__v
    
    def uAxis(self):
        if self.__u is None:
            u1 = 2.0 * self.__q1 * self.__q3 + 2.0* self.__q0 * self.__q2
            u2 = 2.0 * self.__q2 * self.__q3 - 2.0 * self.__q0 * self.__q1
            u3 = self.__q0 * self.__q0 - self.__q1 * self.__q1 - self.__q2 * self.__q2 + self.__q3 * self.__q3
            self.__u = np.array([u1, u2, u3])
        return self.__u

class Atom(object):
    """ A Class for storing atom information """

    def __init__(self, id, type, rx, ry, rz, ix = 0, iy = 0, iz = 0, q0 = 0, q1 = 0, q2 = 0, q3 = 0):
        """ Initialise the class """
        self.id = id                     # id of the atom
        self.type = type                   # type of the atom
        self.x = np.array([rx, ry, rz], dtype = np.float64)     # position of the atom
        self.image = np.array([ix, iy, iz], dtype = np.int32)         # image flags for atoms
        self.q = Quaternion(q0, q1, q2, q3)
        self.unwrap_flag = False

    def sep(self, B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 + 
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x vector from self.x vector """
        return self.x - B.x

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self, L):
        """ Unwraps the coordinates for periodic box, and overwrites x """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.x[j] = self.x[j] + self.image[j]*L[j] # unwrap
            self.unwrap_flag = True


class Processor(object):

    def __init__(self):
        self.__file = self.__N = self.__dump = self.__atoms = self.__L = None
    
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
        self.__atoms = []
        self.__L = []

        for i in range(9):
            line = self.__dump.readline().split()
            if i >= 5 and i <= 7:
                # get the box size
                self.__L.append(np.float64(line[1]) - np.float64(line[0]))
        
        for i in range(self.__N):
            line = self.__dump.readline().split()
            a = Atom(*np.array(line[0:12], dtype=np.float64))
            a.unwrap(self.__L)
            self.__atoms.append(a)

        self.__atoms.sort(key=operator.attrgetter('id'))

    def process(self, outputFile = 'r_g_3.dat'):
        if self.__file is None:
            return
        
        noFrames = int(self.noLines() / (self.__N + 9))
        
        out = open(outputFile, 'w')  
        out.write("# frame number, radius of gyration\n")
        
        self.__dump = open(self.__file, 'r')
        
        for frame in range(noFrames):
            self.readFrame()
            rG = self.radiusOfGyration()
            print(self.twist())
            out.write("%i %.5f\n" % (frame + 1, rG))

        self.__dump.close()
        out.close()

    def radiusOfGyration(self):
        r_mean = np.array([0.0, 0.0, 0.0])
        # Unwrap each atom and calculate the CoM by summing the position of all atoms, and
        # dividing by the number of atoms
        for a in self.__atoms:
            r_mean += a.x
        r_mean /= self.__N
        r_mean = Atom(-1, -1, r_mean[0], r_mean[1], r_mean[2])
        
        r_2 = 0.0
        # Radius of gyration squared is the sum of the squares of the separation of each atom and the CoM
        for a in self.__atoms:
            r_2 += a.sep(r_mean) ** 2
        return np.sqrt(r_2 / self.__N)

    def twist(self):
        tangents = []
        for i in range(0, self.__N - 1):
            tangents.append(unit(self.__atoms[i + 1].x - self.__atoms[i].x))
        
        tangents.append(unit(self.__atoms[0].x - self.__atoms[self.__N - 1].x))

        perpendiculars = []
        for i in range(0, self.__N):
            perpendiculars.append(unit(cross(tangents[i], self.__atoms[i].q.fAxis())))
        
        binomials = [unit(cross(tangents[self.__N - 1], tangents[0]))] 
        angles = [np.arccos(round(dot(tangents[self.__N - 1], tangents[0]), 10))]
        mshift = [np.matmul(rotation_matrix(binomials[0], angles[0]), perpendiculars[self.__N - 1])]
        for i in range(1, self.__N):
            binomials.append(unit(cross(tangents[i - 1], tangents[i])))
            angles.append(np.arccos(round(dot(tangents[i - 1], tangents[0]), 10)))
            mshift.append(np.matmul(rotation_matrix(binomials[i], angles[i]), perpendiculars[i - 1]))

        Tw = 0
        for i in range(0, self.__N - 1):
            Tw += dot(tangents[i], cross(mshift[i], perpendiculars[i]))
        
        return Tw / (2.0 * np.pi)

p = Processor()
p.set('dump3.DNA', 100)
p.process()