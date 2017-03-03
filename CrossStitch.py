import string
import sys
import math


NUMBER_OF_LETTERS = len(string.ascii_lowercase)


def eprint(s="", end="\n"):
    sys.stderr.write(s + end)
    #pass

def get_letter_index(letter):
    """
    Convert letter (char) in to number from 0 (a) to 25 (z).
    "a" is represented by 97 in ascii.
    """
    letter_index = ord(letter) - 97 
    return letter_index

def revert_get_letter_index(letter_index):
    return chr(letter_index + 97)
    


class Point:
    def __init__(self, x=-1, y=-1, l=""):
        self.x = x
        self.y = y
        self.l = l # l stands for Letter

    def __str__(self):
        return "x: %2s, y: %2s, l: %2s" % (str(self.x), str(self.y), self.l)

    @staticmethod
    def distance(p1, p2):
        xd = math.pow(p1.x - p2.x, 2)
        yd = math.pow(p1.y - p2.y, 2)

        return math.sqrt( xd + yd )


class PointCollection:
    def __init__(self):
        self.point_collection = [[] for i in range(NUMBER_OF_LETTERS)]

    def __getitem__(self, index):
        return self.point_collection[index]

    def eprint(self):
        for i in range(NUMBER_OF_LETTERS):
            eprint("Printing points for letter: %s" % revert_get_letter_index(i))
            for j in range(len(self.point_collection[i])):
                eprint(str(self.point_collection[i][j]))

    def fill_collection(self, pattern):

        self.n = len(pattern)

        for i in range(self.n):
            for j in range(self.n):

                if pattern[i][j] != ".":
                    letter_index = get_letter_index(pattern[i][j])
                    p = Point(i, j, pattern[i][j])       
                    self.point_collection[letter_index].append(p)


    def get_path_length_in_collection(self, letter_index):
        n_pixels = len(self.point_collection[letter_index])
        
        if n_pixels == 0:
            return -1.0

        length = 0.0
        for i in range(1, n_pixels):
            p1 = self.point_collection[letter_index][i - 1]
            p2 = self.point_collection[letter_index][i]

            d = Point.distance(p1, p2)
            
            length = length + d

        return length


    def eprint_collection_path_lengths(self):
        
        for letter_index in range(NUMBER_OF_LETTERS):
            pl = self.get_path_length_in_collection(letter_index)
            letter = revert_get_letter_index(letter_index)
            eprint(letter + " path length: " + str(pl))


    def minimize_path_lengths(self):
        pass


    def


class CrossStitch:
    def embroider(self, pattern):

        N = len(pattern)
        # point_collection = [[] for i in range(NUMBER_OF_LETTERS)]
        point_collection = PointCollection()
        point_collection.fill_collection(pattern)

        # point_collection.eprint()

        point_collection.eprint_collection_path_lengths()

        ret = []
        S = len(pattern)
        for col in string.ascii_lowercase:

            first = True
            for r in range(S):
                for c in range(S):
                    if pattern[r][c] == col:
                        if first:
                            first = False
                            ret.append(col)
                        ret.append(str(r+1) + " " + str(c))
                        ret.append(str(r) + " " + str(c+1))
                        ret.append(str(r+1) + " " + str(c+1))
                        ret.append(str(r) + " " + str(c))

        #for r in ret:
        #    if r == "a" or r == "b":
        #        eprint("\n" + r)
        #    else:
        #        eprint(r, " ")
        #eprint()
        #eprint(str(ret))

        return ret

# -------8<------- end of solution submitted to the website -------8<-------

S = int(input())
pattern = []
for i in range(S):
    pattern.append(input().strip())

cs = CrossStitch()
ret = cs.embroider(pattern)
print(len(ret))
for st in ret:
    print(st)
sys.stdout.flush()




