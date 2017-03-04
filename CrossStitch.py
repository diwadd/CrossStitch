import string
import sys
import math
import random
import time

LOCAL_RANDOM_SEED = 1
LOCAL_INFINITY = 4294967296 - 1 # 2**32 - 1
NUMBER_OF_LETTERS = len(string.ascii_lowercase)

#random.seed(LOCAL_RANDOM_SEED)

# Due to the website's acceptance of only a single file
# all functions and classes are placed in one file.

def eprint(s="", end="\n", verbose=True):
    """
    This is a custom print function since
    only logging into the standard error is allowed.
    """
    if verbose == True:
        sys.stderr.write(s + end)
    else:
        pass

def measure_time(func):
    """
    Decorator to measure the execution time of a function.
    """
    def func_wrapper(*args, **kwargs):
        start = time.time()
        ret = func(*args, **kwargs)
        stop = time.time()
        eprint("Time spent in function " + func.__name__ + ": " + str(stop - start))
        return ret
    return func_wrapper

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
    def __init__(self, x=-1, y=-1, l="", b=-1):
        self.x = x
        self.y = y
        self.l = l # l stands for Letter
        self.b = b # b stands for Badge (makes the points distinguishable)

    def __str__(self):
        return "x: %2s, y: %2s, l: %2s, b: %2s" % (str(self.x), str(self.y), self.l, str(self.b))

    @staticmethod
    def distance(p1, p2):
        xd = math.pow(p1.x - p2.x, 2)
        yd = math.pow(p1.y - p2.y, 2)

        return math.sqrt( xd + yd )

    @staticmethod
    def manhattan_distance(p1, p2):
        xd = abs(p1.x - p2.x)
        yd = abs(p1.y - p2.y)

        return xd + yd



class PointCollection:
    def __init__(self):
        self.point_collection = [[] for i in range(NUMBER_OF_LETTERS)]

    def __getitem__(self, index):
        return self.point_collection[index]

    def eprint(self):
        for letter_index in range(NUMBER_OF_LETTERS):
            eprint("Printing points for letter: %s" % revert_get_letter_index(letter_index))
            n_points = len(self.point_collection[letter_index])

            eprint("Number of points: " + str(n_points))
            for j in range(n_points):
                eprint(str(self.point_collection[letter_index][j]))
    
    def eprint_square_matrix(self, matrix):
        n = len(matrix)
        eprint("Matrix size: %d x %d" % (n, n))
        for i in range(n):
            for j in range(n):
                eprint("%4s " % (str(round(matrix[i][j],1))), " ")
            eprint()


    def get_distance_matrix(self, path):
        n_points = len(path)
        distance_matrix = [[0.0 for j in range(n_points)] for i in range(n_points)]

        for i in range(n_points):
            for j in range(n_points):
                if i <= j:
                    continue
                d = Point.distance(path[i], path[j])
                distance_matrix[i][j] = d
                distance_matrix[j][i] = d

        return distance_matrix

    @measure_time
    def fill_collection(self, pattern):

        self.n = len(pattern)

        for i in range(self.n):
            for j in range(self.n):

                if pattern[i][j] != ".":
                    letter_index = get_letter_index(pattern[i][j])
                    p = Point(i, j, pattern[i][j])       
                    self.point_collection[letter_index].append(p)

            
        # Calculate distance matrixes
        self.distance_matrixes = [None for i in range(NUMBER_OF_LETTERS)]
        for letter_index in range(NUMBER_OF_LETTERS):
            path = self.point_collection[letter_index]            
            n_points = len(path)

            if n_points == 0:
                continue

            # Badge the points (make them distinguishable).
            # During the calculation the points are swaped
            # between each other. To make use of the distance matrix
            # we need to distinguish them.
            for j in range(n_points):
                self.point_collection[letter_index][j].b = j

            self.distance_matrixes[letter_index] = self.get_distance_matrix(path) 
            # self.eprint_square_matrix(self.distance_matrixes[letter_index])


    def get_path_length_in_collection(self, letter_index):
        n_points = len(self.point_collection[letter_index])
        
        if n_points == 0:
            return -1.0

        length = 0.0
        for i in range(1, n_points):
            p1 = self.point_collection[letter_index][i - 1]
            p2 = self.point_collection[letter_index][i]

            d = Point.distance(p1, p2)
            
            length = length + d

        return length


    def eprint_collection_path_lengths(self, verbose=False):
        
        apl = 0.0 # Average Path Length
        nvp = 0 # Number of Valid Paths
        for letter_index in range(NUMBER_OF_LETTERS):
            pl = self.get_path_length_in_collection(letter_index)
            if pl == -1:
                continue

            apl = apl + pl
            nvp = nvp + 1

            letter = revert_get_letter_index(letter_index)
            
            eprint(letter + " path length: " + str(pl), "\n" ,verbose)

        eprint("Average path length: " + str(apl/nvp))

    def nearest_neighbour_single_path_optimize(self, path, letter_index):
        """
        Perform a nearest neighbour optimization of a single path.
        """

        n_points = len(path)

        if n_points == 0:
            return path
            
        minimized_path = [Point() for i in range(n_points)]
        visited_points = [False for i in range(n_points)]

        # The starting point is chosen at random.
        # It is also the first current point.
        cpi = random.randint(0, n_points - 1) # Current Point Index
        nocp = 0 # Number Of Checked Points
        current_point = path[cpi]
        minimized_path[nocp] = current_point
        visited_points[cpi] = True

        while (nocp != n_points - 1):
            
            min_d = LOCAL_INFINITY
            for i in range(n_points):
                if visited_points[i] == True:
                    continue
                
                cpb = current_point.b
                p2b = path[i].b

                d = self.distance_matrixes[letter_index][cpb][p2b]
                if d < min_d:
                    min_d = d
                    cpi = i

            current_point = path[cpi]
            visited_points[cpi] = True
            nocp = nocp + 1
            minimized_path[nocp] = current_point

        #self.point_collection[letter_index] = minimized_path

        return minimized_path

    @measure_time
    def two_opt_single_path_optimize(self, path, letter_index):
        """
        Perform a 2-opt optimization of a single path.

        Path should be a single entry (list) from self.point_collection
        i. e. when passing path = self.point_collection[letter_index]
        """

        imp = 0
        n_points = len(path)

        if n_points == 0:
            return path

        while imp < 2:

            for i in range(1, n_points - 1):
                for j in range(i + 1, n_points):
       
                    d_before_prev = 0.0
                    d_before_post = 0.0

                    d_after_prev = 0.0
                    d_after_post = 0.0

                    if j != (n_points - 1):
                        i1b = path[i - 1].b 
                        ib = path[i].b
                        d_before_prev = self.distance_matrixes[letter_index][i1b][ib]

                        jb = path[j].b 
                        j1b = path[j + 1].b
                        d_before_post = self.distance_matrixes[letter_index][jb][j1b]


                        i1b = path[i - 1].b
                        jb = path[j].b
                        d_after_prev = self.distance_matrixes[letter_index][i1b][jb]

                        ib = path[i].b
                        j1b = path[j + 1].b
                        d_after_post = self.distance_matrixes[letter_index][ib][j1b]

                    else:
                        i1b = path[i - 1].b 
                        ib = path[i].b
                        d_before_prev = self.distance_matrixes[letter_index][i1b][ib]

                        i1b = path[i - 1].b
                        jb = path[j].b
                        d_after_prev = self.distance_matrixes[letter_index][i1b][jb]


                    d1 = d_before_prev + d_before_post
                    d2 = d_after_prev + d_after_post

                    if (d1 - d2) > 0.0:
                        imp = 0

                        # We are reversing part of a list.
                        # e. g. a[n:m] = a[(m - 1):(n - 1)]
                        # Thus, we subtract ones from both indexes.
                        path[i:(j + 1)] = path[(j + 1 - 1):(i - 1):-1]
        
            imp = imp + 1

        return path

    @measure_time
    def nearest_neighbour_all_path_minimize(self):
        """
        Perform the nearest neigbour optimization on all paths in self.point_collection.
        """
        for letter_index in range(NUMBER_OF_LETTERS):
            path = self.point_collection[letter_index]
            self.point_collection[letter_index] = self.nearest_neighbour_single_path_optimize(path, letter_index)


    @measure_time
    def two_opt_all_path_minimize(self):

        for letter_index in range(NUMBER_OF_LETTERS):
            path = self.point_collection[letter_index]
            self.point_collection[letter_index] = self.two_opt_single_path_optimize(path, letter_index)

    """
    Stitch rules and defintions

    TL     TR
    +-----+
    |     |
    |     |
    |     |
    +-----+
    BL     BR

    TL - Top Left
    TR - Top Right
    BL - Bottom Left
    BR - Bottom Right

    Possible stitchs: 
    1. BL->TR->TL->BR
    2. BL->TR->BR->TL

    3. TR->BL->TL->BR
    4. TR->BL->BR->TL

    5. TL->BR->BL->TR
    6. TL->BR->TR->BL

    7. BR->TL->BL->TR
    8. BR->TL->TR->BL

    """    

    def x_y_to_corners(self, x, y):
        """
        x and y always mark the Top Left (TL) corner.
        """
        tr_x = x + 1
        tr_y = y

        bl_x = x
        bl_y = y + 1

        br_x = x + 1
        br_y = y + 1

        tl_x = x
        tl_y = y

        corners = (tl_x, tl_y, 
                   tr_x, tr_y,
                   bl_x, bl_y,
                   br_x, br_y,)

        return corners

    def tr_bl_br_tl(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (tr_x, bl_x, br_x, tl_x)
        stitch_y = (tr_y, bl_y, br_y, tl_y)

        return stitch_x, stitch_y


    def tr_bl_tl_br(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (tr_x, bl_x, tl_x, br_x)
        stitch_y = (tr_y, bl_y, tl_y, br_y)

        return stitch_x, stitch_y

    def tl_br_bl_tr(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (tl_x, br_x, bl_x, tr_x)
        stitch_y = (tl_y, br_y, bl_y, tr_y)

        return stitch_x, stitch_y

    def tl_br_tr_bl(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (tl_x, br_x, tr_x, bl_x)
        stitch_y = (tl_y, br_y, tr_y, bl_y)

        return stitch_x, stitch_y

    def br_tl_tr_bl(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (br_x, tl_x, tr_x, bl_x)
        stitch_y = (br_y, tl_y, tr_y, bl_y)

        return stitch_x, stitch_y

    def br_tl_bl_tr(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (br_x, tl_x, bl_x, tr_x)
        stitch_y = (br_y, tl_y, bl_y, tr_y)

        return stitch_x, stitch_y


    def bl_tr_tl_br(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (bl_x, tr_x, tl_x, br_x)
        stitch_y = (bl_y, tr_y, tl_y, br_y)

        return stitch_x, stitch_y

    def bl_tr_br_tl(self, x, y):

        tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y = self.x_y_to_corners(x, y)
        
        # Permute the output
        stitch_x = (bl_x, tr_x, br_x, tl_x)
        stitch_y = (bl_y, tr_y, br_y, tl_y)

        return stitch_x, stitch_y


    @measure_time
    def make_simple_stitchs(self):
        
        ret = []
        for letter_index in range(NUMBER_OF_LETTERS):

            path = self.point_collection[letter_index]
            n_points = len(path)

            if n_points == 0:
                continue

            letter = revert_get_letter_index(letter_index)
            ret.append(letter)

            # eprint("Number of points to append x 4: " + str(n_points*4))

            fp_x = -1
            fp_y = -1

            stitich_functions = [self.tr_bl_tl_br,
                                 self.tl_br_bl_tr,
                                 self.tl_br_tr_bl,
                                 self.br_tl_tr_bl,
                                 self.br_tl_bl_tr,
                                 self.bl_tr_tl_br,
                                 self.bl_tr_br_tl]

            # This is the simplest most possible stitching method.
            # Requires optimization.
            for p in path:

                stitch_x, stitch_y = self.tr_bl_br_tl(p.x, p.y)

                # The starting point of a new stitch must be different from
                # the ending point of the previous stitch.
                if (stitch_x[0] == fp_x) and (stitch_y[0] == fp_y):
                    index = 0
                    while True:
                        stitch_x, stitch_y = stitich_functions[index](p.x, p.y)
                        if (stitch_x[0] != fp_x) and (stitch_y[0] != fp_y):
                            break
                        index = index + 1

                ret.append(str(stitch_x[0]) + " " + str(stitch_y[0]))
                ret.append(str(stitch_x[1]) + " " + str(stitch_y[1]))
                ret.append(str(stitch_x[2]) + " " + str(stitch_y[2]))
                ret.append(str(stitch_x[3]) + " " + str(stitch_y[3]))
                
                fp_x = stitch_x[3]
                fp_y = stitch_y[3]

        return ret


class CrossStitch:
    
    @measure_time
    def embroider(self, pattern):

        N = len(pattern)
        # point_collection = [[] for i in range(NUMBER_OF_LETTERS)]
        point_collection = PointCollection()
        point_collection.fill_collection(pattern)
        
        eprint("\n Print initial paths \n")
        #point_collection.eprint()
        point_collection.eprint_collection_path_lengths()


        point_collection.nearest_neighbour_all_path_minimize()
        eprint("\n Print minimized paths (nearest neighbour) \n")
        #point_collection.eprint()
        point_collection.eprint_collection_path_lengths()

        #point_collection.two_opt_all_path_minimize()
        #eprint("\n Print minimized paths (2 - opt) \n")
        #point_collection.eprint()
        #point_collection.eprint_collection_path_lengths()
 
        ret = []
        flag = True

        if flag == True:        
            ret = point_collection.make_simple_stitchs()
        else:
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


        return ret

# -------8<------- end of solution submitted to the website -------8<-------

S = int(raw_input())
pattern = []
for i in range(S):
    pattern.append(raw_input().strip())

cs = CrossStitch()
ret = cs.embroider(pattern)
print(len(ret))
for st in ret:
    print(st)
sys.stdout.flush()
