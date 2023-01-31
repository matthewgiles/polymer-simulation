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


class Atom:
    """ A Class for storing atom information """

    def __init__(self):
        """ Initialise the class """
        self.id = 0                              # id of the atom
        self.type = 0                            # type of the atom
        self.x = np.array([0.0,0.0,0.0],dtype=np.float64)     # position of the atom
        self.image = np.array([0,0,0],dtype=np.int32)         # image flags for atoms
        self.unwrap_flag = False

    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 + 
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x vector from self.x vector """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
            AminusB[i] = self.x[i] - B.x[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self):
        """ Unwraps the coordinates for periodic box, and overwrites x """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.x[j] = self.x[j] + self.image[j]*L[j] # unwrap
            unwrap_flag = True



def readframe_unwrap(infile,N):
    """ Read a single frame of atoms from a dump file 
        Rescale coordinates to be in rnage -L/2 to L/2
        DOES NOT Unwrap corrdinates for periodic box """

    atoms = []
    L = []

    # read in the 9 header lines of the dump file
    for i in range(9):
        line = infile.readline()
        if i==5 or i==6 or i==7:
            # get the box size
            line = line.split()
            L.append( np.float64(line[1]) - np.float64(line[0]) )

    # now read the atoms
    for i in range(N):
        line = infile.readline()
        line = line.split()
        newatom = Atom()
        newatom.id = int(line[0])
        newatom.type = int(line[1])
        for j in range(3):
            newatom.x[j] = np.float64(line[j+2]) # scale
        atoms.append(newatom)

    # make sure atoms are sorted by id
    atoms.sort(key=operator.attrgetter('id'))

    return atoms,L


def lines_in_file(filename):
    """ Get the number of lines in the file """

    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def radius_of_gyration(atoms,L):
    """ Calculate the radius of gytation -- Rg^2 = (1/N) sum ( r_k - r_mean )^2 
    remember to unwrap periodic boundaries """

    # add code here
    
    return 0.0




############################################################################
#### Start of the main program

dumpfilename = 'dump.DNA'   # this is hardcoded here, could instead read as command line argument
Natoms = 200   # this is hardcoded here, could instead read as command line argument

outfile_Rg = 'r_g_1.dat'

Nlines = lines_in_file(dumpfilename)  # get length of file
Nframes = int(Nlines / (Natoms+9))         # there are 9 header lines in each frame

# open the intput file
inf = open(dumpfilename, 'r')  

# open the output files and print a header
ouf_rg = open(outfile_Rg, 'w')  
ouf_rg.write("# frame number, radius of gyration\n")

# go through the file frame by frame
for frame in range(Nframes):
    # read the frame, unwrapping periodic coordinates
    atoms, L = readframe_unwrap(inf,Natoms)

    # unwarp period boundary coordinates -- needed for radius of gyration
    for i in range(len(atoms)):
        atoms[i].unwrap()

    # calculate radius of gyration
    Rg = radius_of_gyration(atoms,L)

    # output some results
    ouf_rg.write( "%i %.5f\n"%(frame+1,Rg) )


# close the files
inf.close()

ouf_rg.close()


# Finished!
