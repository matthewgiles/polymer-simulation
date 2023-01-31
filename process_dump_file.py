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


class Atom(object):
    """ A Class for storing atom information """

    def __init__(self, id, type, rx, ry, rz):
        """ Initialise the class """
        self.id = id                             # id of the atom
        self.type = type                         # type of the atom
        self.r = np.array([rx, ry, rz], dtype = np.float64)     # position of the atom
        self.image = np.array([0,0,0], dtype = np.int32)         # image flags for atoms
        self.unwrap_flag = False

    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.r[0]-B.r[0])**2 + 
                        (self.r[1]-B.r[1])**2 + (self.r[2]-B.r[2])**2 )

    def minus(self,B):
        """ Subtract B.r vector from self.r vector """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
            AminusB[i] = self.r[i] - B.r[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self, L):
        """ Unwraps the coordinates for periodic box, and overwrites x """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.r[j] = self.r[j] + self.image[j]*L[j] # unwrap
            unwrap_flag = True


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
            a = Atom(id = line[0], type = line[1], rx = line[2], ry = line[3], rz = line[4])
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
            rG = self.radiusOfGyration(atoms, L)
            out.write("%i %.5f\n" % (frame + 1, rG))

        self.__dump.close()
        out.close()

    def radiusOfGyration(self, atoms, L):
        return 0.0

p = Processor()
p.set('dump.DNA', 100)
p.process()